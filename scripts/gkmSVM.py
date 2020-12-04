##!/usr/bin/env python
"""
    gkmSVM.py: Kernel manipulation, training and test with sklearn lib
    (internal algorithm is libSVM)

    Copyright (C) 2020 Seong Kyu Han

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import os, sys, pickle
import ctypes, random, logging
import numpy as np
from sklearn.svm import SVC
from sklearn.model_selection import KFold
from sklearn.metrics import roc_curve, auc
from sklearn.ensemble import GradientBoostingRegressor
from itertools import repeat


'''
struct svm_parameter
{
    int svm_type;
    int data_type;
    int kernel_type;
    int L;
    int k;
    int d;
    uint8_t M;
    double H;
    double gamma;

    /* these are for training only */
    double cache_size; /* in MB */
    double eps; /* stopping criteria */
    double C;   /* for C_SVC, EPSILON_SVR and NU_SVR */
    int nr_weight;      /* for C_SVC */
    int *weight_label;  /* for C_SVC */
    double* weight;     /* for C_SVC */
    double nu;  /* for NU_SVC, ONE_CLASS, and NU_SVR */
    double p;   /* for EPSILON_SVR */
    int shrinking;  /* use the shrinking heuristics */
    int probability; /* do probability estimates */
};
'''

##
# Options for GKM Kernel 
# Ctype-compatible python class
##
class svmParam(ctypes.Structure):
    _fields_ = (
        ('svm_type', ctypes.c_int),
        ('data_type', ctypes.c_int),
        ('kernel_type', ctypes.c_int),
        ('L', ctypes.c_int),     # L = 10 ## default values
        ('k', ctypes.c_int),     # k = 6
        ('d', ctypes.c_int),     # d = 3
        ('M', ctypes.c_uint8),
        ('H', ctypes.c_double),
        ('gamma', ctypes.c_double),
    )

##
# Compute Kernel Matrix
# (using optimized algorithm in gkmSVM-2.0)
##
def computeGkmKernel(opts_arg):
    svm_param = svmParam(*opts_arg[:len(svmParam._fields_)])
    nseq = svm_param.maxnumseq # max_n = 15,000 (= 7500 x 2)

    # create blank numpy obj with aligned memory addr
    # (to be compatible with double ** in ctype func)
    # enables call by reference for np.mat variable
    kmat = np.zeros(shape=(nseq, nseq))
    kmat_p = (kmat.ctypes.data + np.arange(kmat.shape[0]) * kmat.strides[0]).astype(np.uintp)
    array_2d_double = np.ctypeslib.ndpointer(dtype=np.uintp, ndim=1, flags='C')

    narr = np.array([0, 0], dtype=np.int)
    narr_p = ctypes.c_void_p(narr.ctypes.data)
    array_1d_int = np.ctypeslib.ndpointer(dtype=np.int, ndim=1, flags='C')

    # call ctype func in ../bin/GkmKernel.so
    libgkm = np.ctypeslib.load_library("libgkmqc.so", "../bin")
    libgkm.gkmKernelPyWrapper.restype = None
    libgkm.gkmKernelPyWrapper.argtypes = (ctypes.POINTER(svmParam), array_2d_double, array_1d_int,)
    libgkm.gkmKernelPyWrapper(svm_param, kmat_p, narr_p)

    n_pseqs = narr[0]
    n_nseqs = narr[1]
    n_seqs = n_pseqs + n_nseqs
    kmat = kmat[:n_seqs, :n_seqs] # shrink kernel matrix fit for pos/neg seqs

    return kmat, n_pseqs, n_nseqs

##
# cross-validate gkm-SVM models 
##
def crossValidate(kmat, n_pseqs, n_nseqs, ncv, nu_auc_regressor):
    # Answer set
    seqids = \
        list(map(lambda x: "p%4d" % x, range(n_pseqs))) +\
        list(map(lambda x: "p%4d" % x, range(n_nseqs)))
    y = np.concatenate(np.repeat(1, n_pseqs), np.repeat(0, n_nseqs))

    # Normal mode: 5-fold cross-validation
    if nu_auc_regressor == None:
        kf = KFold(n_splits=ncv)
        tprs = []
        aucs = []
        mean_fpr = np.linspace(0, 1, 100)

        for trainIdx, testIdx in kf.split(seqids):
            y_train, y_test = y[trainIdx], y[testIdx]
            k_train = kmat[trainIdx, :][:, trainIdx]
            k_test  = kmat[trainIdx, :][:, testIdx]
            sv = SVC(kernel="precomputed", C=1.0, tol=1e-3, shrinking=False, gamma=1.0, cache_size=256)
            y_prob_ = sv.fit(k_train, y_train).predict_proba(k_test)
            #y_pred = sv.predict(k_test)

            # interpolate tpr according to fpr
            fpr, tpr, _ = roc_curve(y_test, y_prob_[:, 1])
            interp_tpr = np.interp(mean_fpr, fpr, tpr)
            interp_tpr[0] = 0.0
            tprs.append(interp_tpr)
            aucs.append(auc(fpr, tpr))

        # calculate mean-AUC
        mean_tpr = np.mean(tprs, axis=0)
        mean_tpr[-1] = 1.0
        auc_score = auc(mean_fpr, mean_tpr)
        #auc_std = np.std(aucs)

    # Fast mode: AUC estimation based on nu values from single training
    # using Gradient-boost regressor
    else:
        sv = SVC(kernel="precomputed", C=1.0, tol=1e-3, shrinking=False, gamma=1.0, cache_size=256) # TODO: fit param
        sv.fit(kmat, y)
        auc_score = nu_auc_regressor.predict(np.atleast_2d([sv.nu]).T)[0]

    return auc_score # (auc_score, auc_std)

##
# init Function
##
def init(opts_arg, nu_auc_regressor=None):
    kmat, n_pseqs, n_nseqs = computeGkmKernel(opts_arg)
    ncv = opts_arg[15] # normal of folds
    auc_score = crossValidate(kmat, n_pseqs, n_nseqs, ncv, nu_auc_regressor)

    return auc_score
