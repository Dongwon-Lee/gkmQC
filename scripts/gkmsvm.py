#!/usr/bin/env python

"""
    gkmsvm.py: Kernel manipulation, training and test with sklearn lib
    (internal algorithm is libSVM)

    Copyright (C) 2020 Seong Kyu Han, Dongwon Lee

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
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from itertools import repeat
from multiprocessing import Pool

##
# nu-auc regressor
##
dir_this   = os.path.dirname(os.path.abspath(__file__))
dir_prnt   = os.path.dirname(dir_this)
base_data_dir = os.path.join(dir_prnt, "data")
bin_dir = os.path.join(dir_prnt, "bin")

#f = open("%s/nu_auc_gb_regressor.pkl" % base_data_dir, "rb")
#nu_auc_regressor = pickle.load(f)
#f.close()

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
    args_gkm[7] = ctypes.c_char_p(bytes(args_gkm[7], 'ascii'))
    args_gkm[8] = ctypes.c_char_p(bytes(args_gkm[8], 'ascii'))
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
    
    _gkmkern_pylib = np.ctypeslib.load_library("gkmkern_pylib.so", bin_dir)
    _gkmkern_pylib.gkm_main_pywrapper.restype = ctypes.c_int
    _gkmkern_pylib.gkm_main_pywrapper.argtypes = (ctypes.POINTER(gkmOpt), array_2d_double, c_int_p)
    ret = _gkmkern_pylib.gkm_main_pywrapper(opts, kmat_p, narr_p)

    if ret:
        print("error on kernel construction")
        sys.exit()

    n_pseqs, n_nseqs = narr
    n_seqs = n_pseqs + n_nseqs
    kmat = kmat[:n_seqs, :n_seqs] # shrink kernel matrix fit for pos/neg seqs
    kmat = np.maximum(kmat, kmat.T)

    return kmat, n_pseqs, n_nseqs

##
# cross-validate gkm-SVM models 
##
def _svm_train_proc(args_svm, y, trainIdx, testIdx):
    regularization, precision, shrinking, cache_size, _,\
    _, _, _, _ = args_svm
    y_train, y_test = y[trainIdx], y[testIdx]
    k_train = kmat[trainIdx, :][:, trainIdx]
    k_test  = kmat[testIdx, :][:, trainIdx]
    sv = SVC(
        kernel="precomputed",
        C=regularization,
        tol=precision,
        shrinking=bool(shrinking),
        gamma=1.0,
        cache_size=cache_size
    )
    y_score = sv.fit(k_train, y_train).decision_function(k_test)
    auc = roc_auc_score(y_test, y_score)
    return auc

def pool_wrapper_svm_train(args):
    return _svm_train_proc(*args)

def crossValidate(args_svm, _kmat, n_pseqs, n_nseqs):

    global kmat
    kmat = _kmat
    _, _, _, _, ncv,\
    repeats, fast_estimation, random_seeds, p = args_svm

    if random_seeds < 0:
        random_seeds = None
    
    seqids = \
        list(map(lambda x: "p%4d" % x, range(n_pseqs))) +\
        list(map(lambda x: "n%4d" % x, range(n_nseqs)))
 
    y = np.concatenate((np.repeat(1, n_pseqs), np.repeat(0, n_nseqs)))
    args_svm = args_svm

    # Normal mode: 5-fold cross-validation
    if fast_estimation == 0:
        args_l = []
        for _ in range(repeats):
            kf = StratifiedKFold(n_splits=ncv, shuffle=True, random_state=random_seeds)
            for trainIdx, testIdx in kf.split(seqids, y):
                args_l.append((args_svm, y, trainIdx, testIdx))

        pool = Pool(p)
        aucs = pool.map(pool_wrapper_svm_train, args_l)
        auc_score = np.mean(aucs)
        auc_std = np.std(aucs)

    # Fast mode: AUC estimation based on nu values from single training
    # using local regressor
    #else:
    #    sv = SVC(
    #        kernel="precomputed",
    #        C=regularization,
    #        tol=precision,
    #        shrinking=bool(shrinking),
    #        gamma=1.0,
    #        cache_size=cache_size
    #    )
    #    sv.fit(kmat, y) 
        #nu = np.sum(np.abs(sv.dual_coef_[0])) / len(y)
        #auc_score = nu_auc_regressor.predict(np.atleast_2d([nu]).T)[0]
            #auc_std = np.nan
    del kmat
    return (auc_score, auc_std)

##
# init Function
##
def init(pos_fa, neg_fa, args):

    args_gkm = [
        args.kernel_type,
        args.full_word_length, # l
        args.non_gap_length, # k
        args.max_num_gaps, # d
        args.init_decay, # M
        args.half_life_decay, # H
        args.rbf_gamma, # gamma
        pos_fa,
        neg_fa,
        args.n_processes,
        args.verbosity
    ]

    kmat, n_pseqs, n_nseqs = computeGkmKernel(args_gkm)

    args_svm = [
        args.regularization,
        args.precision,
        args.shrinking,
        args.cache_size,
        args.ncv,
        args.repeats,
        args.fast_estimation,
        args.random_seeds,
        args.n_processes
    ]
    
    auc_score, auc_std = crossValidate(args_svm, kmat, n_pseqs, n_nseqs)
    
    # write results to out file
    eval_out_file = args.name + ".gkmqc.eval.out"
    fa = open(eval_out_file, "a")
    fa.write('\t'.join(map(str, [pos_fa, neg_fa, n_pseqs, auc_score, auc_std])) + "\n")
    fa.close()

import argparse

def main():
    print("start gkm-SVM as main processor")

    desc_txt = "\n".join([
        "lsgkm-pywrapper with fast kernel-mat construction",
        "sequence-based predictive model with gapped-kmer kernel (Lee 2016).",
        "LIBSVM (Chang & Lin 2011) was used for implementing SVC.",
        "-- Seong Kyu Han (seongkyu.han@childrens.harvard.edu),",
        "-- Dongwon Lee (dongwon.lee@childrens.harvard.edu)"
    ])

    parser = argparse.ArgumentParser(description=desc_txt,
        formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("-p", "--pos-fa", type=str, required=True,
        help="positive fa file. REQUIRED")
    parser.add_argument("-n", "--neg-fa", type=str, required=True,
        help="negative fa file. REQUIRED")
    parser.add_argument("-w", "--name", type=str, required=True,
        help="prefix of output file to write AUC score. REQUIRED")
    parser.add_argument("-s", "--random-seeds", type=int, default=-1,
        help="random seed number\nfor reproducibility\n(default: no seed)")
    parser.add_argument("-@", "--n-processes", type=int, default=1,
        help="number of processes\n(default: 1)")
    parser.add_argument("-v", "--verbosity", type=int, default=1,
        help="verbosity\n(default: 1), 0: silent")

    group_gkm = parser.add_argument_group('gkm-kernel')
    group_gkm.add_argument("-t", "--kernel-type", type=int, default=4,
        help="\n".join([
            "kernel function\n(default: 4 wgkm)",
            "NOTE: RBF kernels (3 and 5) work best with -c 10 -g 2",
            "0 -- gapped-kmer",
            "1 -- estimated l-mer with full filter",
            "2 -- estimated l-mer with truncated filter (gkm)",
            "3 -- gkm + RBF (gkmrbf)",
            "4 -- gkm + center weighted (wgkm)",
            "     weight = max(M, floor(M*exp(-ln(2)*D/H)+1))",
            "5 -- gkm + center weighted + RBF (wgkmrbf)"
        ]))        
    group_gkm.add_argument("-L", "--full-word-length", type=int, default=10,
        help="full word length including gaps, 3<=L<=12\n(default: 10)")
    group_gkm.add_argument("-k", "--non-gap-length", type=int, default=6,
        help="number of non-gap positions, k<=L\n(default: 6)")
    group_gkm.add_argument("-d", "--max-num-gaps", type=int, default=3,
        help="maximum number of gaps allowed, d<=min(4, L-k)\n(default: 3)")
    group_gkm.add_argument("-M", "--init-decay", type=int, default=50,
        help="\n".join([
            "the initial value (M) of the exponential decay function",
            "for wgkm-kernels. max=255, -t 4 or 5 only\n(default: 50)"
        ]))
    group_gkm.add_argument("-H", "--half-life-decay", type=int, default=50,
        help="\n".join([
            "set the half-life parameter (H) that is the distance (D) required",
            "to fall to half of its initial value in the exponential decay",
            "function for wgkm-kernels. -t 4 or 5 only\n(default: 50)"
        ]))
    group_gkm.add_argument("-G", "--rbf-gamma", type=float, default=1.0,
        help="gamma for RBF kernel. -t 3 or 5 only\n(default: 1.0)")

    ## svm training options
    group_svm = parser.add_argument_group('SVM training')
    group_svm.add_argument("-C", "--regularization", type=float, default=1.0,
        help="regularization parameter C\n(default: 1.0)")
    group_svm.add_argument("-e", "--precision", type=float, default=0.001,
        help="set the precision parameter epsilon\n(default: 0.001)")
    group_svm.add_argument("-u", "--shrinking", type=int, default=0,
        help="if set, use the shrinking heuristics\n(default: 0)")
    group_svm.add_argument("-c", "--cache-size", type=int, default=512,
        help="cache memory size in MB\n(default: 512MB)")
    group_svm.add_argument("-x", "--ncv", type=int, default=5,
        help="x-fold cross validation\nfor estimating effects of tags in training set\n(default: 5)")
    group_svm.add_argument("-r", "--repeats", type=int, default=1,
        help="number of repeats of CV training\nto reduce random variation\n(default: 1)")
    group_svm.add_argument("-f", "--fast-estimation", type=int, default=0,
        help="fast estimation of AUC without nCV:\nusing nu score from trained SVM\n(default: 0)")

    # args processing codes
    args = parser.parse_args()
    init(args.pos_fa, args.neg_fa, args)

if __name__ == '__main__':
    main()
