#!/usr/bin/env python

"""
    gkmQC: gapped k-mer-SVM Quality Check

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

import sys
import argparse
import os.path
import math
import logging
from itertools import imap
import pickle
from sklearn.ensemble import GradientBoostingRegressor
from ctypes import *

__version__ = '1.0.0'

##
# nu-auc regressor - will move to gkmQC.py
##
base_data_dir = os.path.abspath('%s/../data' % os.path.dirname(__file__))
f = open("%s/nu_auc_gb_regressor.pkl" % base_data_dir, "rb")
nu_auc_regressor = pickle.load(f)
f.close()


HEADER  = "\n# ======================================="
HEADER += "\n#   gapped k-mer-SVM Quality Check (gkmQC)"
HEADER += "\n#   Version {0}".format(__version__)
HEADER += "\n#   (C) 2020 Seong Kyu Han, Dongwon Lee"
HEADER += "\n#   GNU General Public License v3"
HEADER += "\n# ======================================="

# Logging template
#logging.info("%d elements have <%d tags after filtering",
#            good_elems.count(False), args.min_num_tags)

def main():
    global HEADER

    desc_txt="perform quality evaluation of open-chromatin peaks\
    sequence-based predictive model  with gapped-kmer kernels (Ghandi et al. 2014; Lee 2016). \
    LIBSVM (Chang & Lin 2011) was used for implementing SVR. \
    -- by Seong Kyu Han (seongkyu.han@childrens.harvard.edu), Dongwon Lee (dongwon.lee@childrens.harvard.edu)"

    desc_nidx_txt = "Build random-sequence index"
    desc_eval_txt = "Evaluate open-chromatin peaks"
    desc_optz_txt = "Optimize open-chromatin peaks"

    # init parser
    parser = argparse.ArgumentParser(description=desc_txt,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparsers = parser.add_subparsers(title='commands', dest='commands')

    subparser_nidx = subparsers.add_parser('buildidx',
            help=desc_nidx_txt,
            description=desc_nidx_txt,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparser_eval = subparsers.add_parser('evaluate',
            help=desc_eval_txt,
            description=desc_eval_txt,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    subparser_optz = subparsers.add_parser('optimize',
            help=desc_optz_txt,
            description=desc_optz_txt,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # parser for the "evaluate" command
    subparser_eval.add_argument("narrowPeak", type=str,
            help="Open-chromatin peaks from calling software")
    subparser_eval.add_argument("nseqIdx", type=str,
            help="Null-seq idx pickle file")
    subparser_eval.add_argument("genome", type=str,
            help="Genome Assembly (e.g., hg38/19, mm10/9")
    subparser_eval.add_argument("-n", "--name", type=str, required=True,
            help="name of output prefix. REQUIRED")
    subparser_eval.add_argument("-t", "--min-num-tags", type=int, default=5,
            help="minimum number of tags per element for SVR training")
    subparser_eval.add_argument("-rs", "--rank-start", type=int, default=1,
            help="the rank number of peak subsets to start evaluation")
    subparser_eval.add_argument("-re", "--rank-end", type=int, default=20,
            help="the rank number of peak subsets to end evaluation")

    # parser for the "optimize" command

    # parser for common (evalute, optimize)
    for subparser_comm in [subparser_eval, subparser_optz]:
        subparser_comm.add_argument("-n", "--name", type=str, required=True,
                help="prefix used in 'build' command for training data set. this will also be used for output prefix. REQUIRED")
        subparser_comm.add_argument("-L", "--full-word-length", type=int, default=8,
                help="full word length including gaps, 3<=L<=12")
        subparser_comm.add_argument("-k", "--non-gap-length", type=int, default=4,
                help="number of non-gap positions, k<=L")
        subparser_comm.add_argument("-d", "--max-num-gaps", type=int, default=4,
                help="maximum number of gaps allowed, d<=min(4, L-k)")
        subparser_comm.add_argument("-R", "--use-revcomp", action='store_true',
                help="if set, reverse-complement tag sequences are also used")
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

    args = parser.parse_args()

    # compatible with clog format..
    logfmt_str = '%(levelname)s %(asctime)s: %(message)s'
    datefmt_str = '%Y-%m-%d %H:%M:%S'

    logging.basicConfig(stream=sys.stdout,
            format=logfmt_str, datefmt=datefmt_str,
            level=logging.INFO)

    HEADER += "\n# Command line:" + ' '.join(sys.argv)
    HEADER += "\n# Optional parameters:"
    if args.commands == "evaluate":
        HEADER += "\n#   EXP_NAME={0}".format(args.name)
        HEADER += "\n#   MIN_DNA_READ_CNTS={0}".format(args.min_dna_read_cnts)
        HEADER += "\n#   MIN_NUM_TAGS={0}".format(args.min_num_tags)
        HEADER += "\n#   LEFT_FLANKING_SEQ={0}".format(args.left_flanking_seq)
        HEADER += "\n#   RIGHT_FLANKING_SEQ={0}".format(args.right_flanking_seq)

    if args.commands == "optimize":
        HEADER += "\n#   EXP_NAME={0}".format(args.name)
        HEADER += "\n#   FULL_WORD_LENGTH={0}".format(args.full_word_length)
        HEADER += "\n#   NON_GAP_LENGTH={0}".format(args.non_gap_length)
        HEADER += "\n#   MAX_NUM_GAPS={0}".format(args.max_num_gaps)
        HEADER += "\n#   USE_REVCOMP={0}".format(True if args.use_revcomp else False)
        HEADER += "\n#   REGULARIZATION={0}".format(args.regularization)
        HEADER += "\n#   EPSILON={0}".format(args.epsilon)
        HEADER += "\n#   PRECISION={0}".format(args.precision)
        HEADER += "\n#   CACHE_SIZE={0}".format(args.cache_size)
        HEADER += "\n#   CV={0}".format(args.cv)
        HEADER += "\n#   RANDOM_SEEDS={0}".format(args.random_seeds)
        HEADER += "\n#   NUM_THREADS={0}".format(args.threads)
        HEADER += "\n#   CV_REPEATS={0}".format(args.repeats)
    
    logging.info(HEADER)

    # input/output files
    TTAGFILE = args.name + ".tag_rexpr.txt"
    TTAG_ADJ_FILE = args.name + ".tag_rexpr_adj.txt"
    RTAGFILE = args.name + ".excluded_tags.txt"
    MODELFILE = args.name + ".adj.model.txt" # model trained on adjusted data`
    RTAGSCOREFILE = args.name + ".excluded_tags.score.txt"
    CVFILE = args.name + ".cv.txt"

    if args.commands == "build":
        logging.info("### build training data")
        build_training_data(args, TTAGFILE, RTAGFILE)

    if args.commands == "train":
        # assuming that libmtsa.so is in the same directory as mtsa.py
        dll_name = os.path.join(os.path.dirname(__file__), "libmtsa.so")
        libmtsa=CDLL(dll_name)

        tdfile=c_char_p(TTAGFILE)
        outprefix=c_char_p(args.name)
        nthreads = c_int(args.threads)
        rseed = c_int(args.random_seeds)
        L = c_int(args.full_word_length)
        k = c_int(args.non_gap_length)
        d = c_int(args.max_num_gaps)
        norc = c_int(0) if args.use_revcomp else c_int(1)
        Cp = c_double(args.regularization)
        p = c_double(args.epsilon)
        eps = c_double(args.precision)
        ncv = c_int(args.cv)
        cache_size = c_double(args.cache_size)

        libmtsa.mtsa_init(2, nthreads)

        logging.info("### 1. perform initial cross-validation")
        libmtsa.mtsa_train_main(tdfile, outprefix, rseed,
                L, k, d, norc, Cp, p, eps, ncv, cache_size)

        logging.info("### 2. perform second cross-validation(s) with adjusted exprs")

        adjust_tag_group_effect(TTAGFILE, CVFILE, TTAG_ADJ_FILE)

        adj_tdfile=c_char_p(TTAG_ADJ_FILE)
        logging.info("repeat cross-valiations %d time(s)", args.repeats)
        #repeat cross validation with different random seeds
        for i in xrange(args.repeats):
            logging.info("repeated cross-valiation #%d", i+1)
            outprefix_adj = c_char_p(args.name + ".adj." + str(i))
            rseed = c_int(args.random_seeds+i)
            libmtsa.mtsa_train_main(adj_tdfile, outprefix_adj, rseed,
                    L, k, d, norc, Cp, p, eps, ncv, cache_size)

        logging.info("### 3. build a model using all data")
        outprefix_adj = c_char_p(args.name + ".adj")
        libmtsa.mtsa_train_main(adj_tdfile, outprefix_adj, rseed,
                L, k, d, norc, Cp, p, eps, 0, cache_size)

        logging.info("### 4. score excluded tags")
        testfile=c_char_p(RTAGFILE)
        modelfile=c_char_p(MODELFILE)
        outfile=c_char_p(RTAGSCOREFILE)
        libmtsa.mtsa_predict_main(testfile, modelfile, outfile)

        logging.info("### 5. calculate the sequence factor for normalization")
        cal_tag_sequence_factor(args)

    if args.commands == "predict":
        dll_name = os.path.join(os.path.dirname(__file__), "libmtsa.so")
        libmtsa=CDLL(dll_name)

        nthreads = c_int(args.threads)
        libmtsa.mtsa_init(2, nthreads)

        logging.info("### score sequences using the trained model")

        testfile=c_char_p(args.input_fn)
        modelfile=c_char_p(MODELFILE)
        outfile=c_char_p(args.output_fn)
        libmtsa.mtsa_predict_main(testfile, modelfile, outfile)

    if args.commands == "normalize":
        logging.info("### normalize mRNA count data")
        normalize_mrna_counts(args)

if __name__=='__main__':
    main()
