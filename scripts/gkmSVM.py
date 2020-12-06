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

##
# Options for GKM Kernel 
# Ctype-compatible python class
##
class gkmOpt(ctypes.Structure):
    _fields_ = (
        ('kernel_type', ctypes.c_int),
        ('L', ctypes.c_int),     # L = 10 ## default values
        ('k', ctypes.c_int),     # k = 6
        ('d', ctypes.c_int),     # d = 3
        ('M', ctypes.c_uint8),
        ('H', ctypes.c_double),
        ('gamma', ctypes.c_double),
        ('posfile', ctypes.c_char_p),
        ('negfile', ctypes.c_char_p),
        ('nthreads', ctypes.c_int),
        ('verbosity', ctypes.c_int),
    )

##
# Compute Kernel Matrix
# (using optimized algorithm in LSGKM)
##
def computeGkmKernel(args_gkm):
    opts = gkmOpt(*args_gkm)

    # create blank numpy obj with aligned memory addr
    # (to be compatible with double ** in ctype func)
    # enables call by reference for np.mat variable
    kmat = np.zeros(shape=(15000, 15000))
    kmat_p = (kmat.ctypes.data + np.arange(kmat.shape[0]) * kmat.strides[0]).astype(np.uintp)
    array_2d_double = np.ctypeslib.ndpointer(dtype=np.uintp, ndim=1, flags='C')

    narr = np.ones(2, dtype=np.int32)
    c_int_p = ctypes.POINTER(ctypes.c_int)
    narr_p = narr.ctypes.data_as(c_int_p)

    # call ctype func in ../bin/GkmKernel.so
    _gkmkern_pylib = np.ctypeslib.load_library("gkmkern_pylib.so", "../bin")
    _gkmkern_pylib.gkm_main_pywrapper.restype = ctypes.c_int
    _gkmkern_pylib.gkm_main_pywrapper.argtypes = (ctypes.POINTER(gkmOpt), array_2d_double, c_int_p)
    ret = _gkmkern_pylib.gkm_main_pywrapper(opts, kmat_p, narr_p)

    if ret:
        sys.exit()

    n_pseqs, n_nseqs = narr
    n_seqs = n_pseqs + n_nseqs
    kmat = kmat[:n_seqs, :n_seqs] # shrink kernel matrix fit for pos/neg seqs

    return kmat, n_pseqs, n_nseqs

##
# cross-validate gkm-SVM models 
##
def crossValidate(args_svm, kmat, n_pseqs, n_nseqs):
    '''
        subparser_comm.add_argument("-C", "--regularization", type=float, default=0.1,
                help="regularization parameter C")
        subparser_comm.add_argument("-p", "--epsilon", type=float, default=0.1,
                help="epsilon in loss function of SVR")
        subparser_comm.add_argument("-e", "--precision", type=float, default=0.001,
                help="precision parameter")
        subparser_comm.add_argument("-M", "--cache-size", type=float, default=512,
                help="cache memory size in MB")
        subparser_comm.add_argument("-x", "--cv", type=int, default=5,
                help="x-fold cross validation for estimating effects of tags in training set")
        subparser_comm.add_argument("-s", "--random-seeds", type=int, default=1,
                help="random seed number for reproducibility of cross-validation")
        subparser_comm.add_argument("-r", "--repeats", type=int, default=1,
                help="number of repeats of CV training to reduce random variation")
        subparser_comm.add_argument("-T", "--threads", type=int, default=1,
                help="number of threads for SVR training; 1, 4, or 16")
    '''
    ncv = args_svm[0]
    nu_auc_regressor = args_svm[1]
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
def init(args):
    args_gkm = args[:11]
    kmat, n_pseqs, n_nseqs = computeGkmKernel(args_gkm)
    
    args_svm = args[11:]
    ncv = args[11] # normal of folds
    auc_score = crossValidate(args_svm, kmat, n_pseqs, n_nseqs)

    return auc_score
