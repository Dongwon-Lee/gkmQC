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
import pickle
import logging

__version__ = '1.0.0'

# namespacing script dir
dir_this   = os.path.dirname(os.path.abspath(__file__))
dir_prnt   = os.path.dirname(dir_this)
dir_scripts = os.path.join(dir_prnt, "scripts")
dir_data = os.path.join(dir_prnt, "data")
sys.path.append(dir_scripts)

HEADER  = "\n# ==========================================="
HEADER += "\n#   gapped k-mer-SVM Quality Check (gkmQC)"
HEADER += "\n#   Version {0}".format(__version__)
HEADER += "\n#   (C) 2020 Seong Kyu Han, Dongwon Lee"
HEADER += "\n#   GNU General Public License v3"
HEADER += "\n# ==========================================="

# Logging template
#logging.info("%d elements have <%d tags after filtering",
#            good_elems.count(False), args.min_num_tags)

def main():
    global HEADER

    desc_txt = "\n".join([
        "Perform quality evaluation of open-chromatin peaks",
        "sequence-based predictive model with gapped-kmer kernel (Lee 2016).",
        "LIBSVM (Chang & Lin 2011) was used for implementing SVC.",
        "-- Seong Kyu Han (seongkyu.han@childrens.harvard.edu),",
        "-- Dongwon Lee (dongwon.lee@childrens.harvard.edu)"
    ])

    desc_nidx_txt = "Build genome index for generating null sequence"
    desc_eval_txt = "Evaluate peaks with AUCs of peak subsets"
    desc_optz_txt = "Optimize peaks with AUC-based thresholding"

    # init parser
    parser = argparse.ArgumentParser(description=desc_txt,
        formatter_class=argparse.RawTextHelpFormatter)

    parser._action_groups.pop()

    subparsers = parser.add_subparsers(title='commands', dest='commands')

    subparser_nidx = subparsers.add_parser('buildidx',
        help=desc_nidx_txt,
        description=desc_nidx_txt,
        formatter_class=argparse.RawTextHelpFormatter)

    subparser_eval = subparsers.add_parser('evaluate',
        help=desc_eval_txt,
        description=desc_eval_txt,
        formatter_class=argparse.RawTextHelpFormatter)

    subparser_optz = subparsers.add_parser('optimize',
        help=desc_optz_txt,
        description=desc_optz_txt,
        formatter_class=argparse.RawTextHelpFormatter)

    subparser_nidx._action_groups.pop()
    subparser_eval._action_groups.pop()
    subparser_optz._action_groups.pop()

    # parser for the "buildidx" command
    group_req_nidx = subparser_nidx.add_argument_group('required arguments')
    group_req_nidx.add_argument("-i", "--chrom-file", type=str, required=True,
       help="path of archived chromosome fa file (e.g., hg38.chromFa.tar.gz).")
    group_req_nidx.add_argument("-g", "--genome-assembly", type=str, required=True,
        help="genome assembly name for prefix annotation.")

    group_opt_nidx = subparser_nidx.add_argument_group('optional arguments')
    group_opt_nidx.add_argument("-w", "--window-bp", type=int, default=600,
        help="size of scanning window\n(default: 600bp)")
    group_opt_nidx.add_argument("-@", "--n-processes", type=int, default=10,
        help="number of processes\n(default: 10)")
    
    # parser for the "evaluate" command
    group_req_eval = subparser_eval.add_argument_group('required arguments')
    group_req_eval.add_argument("-i", "--peak-file", type=str, required=True,
       help="peak-calling file; *.narrowPeak.")
    group_req_eval.add_argument("-n", "--name", type=str, required=True,
        help="name of output prefix.")
    group_req_eval.add_argument("-g", "--genome-assembly", type=str, required=True,
        help="prefix used in 'buildidx' command for genome-indexing.")

    group_opt_eval = subparser_eval.add_argument_group('optional arguments\nbase')
    group_opt_eval.add_argument("-rs", "--rank-start", type=int, default=1,
        help="rank number of peak subsets to start evaluation\n(default: 1)")
    group_opt_eval.add_argument("-re", "--rank-end", type=int, default=20,
        help="rank number of peak subsets to end evaluation.\n(default: 20) set 0 not to limit.")

    # parser for the "optimize" command
    group_req_optz = subparser_optz.add_argument_group('required arguments')
    group_req_optz.add_argument("-i", "--gkmqc-file", type=str, required=True,
        help="peak-calling file or gkmQC AUC output file.")
    group_req_optz.add_argument("-n", "--name", type=str, required=True,
        help="name of output prefix.")
    group_req_optz.add_argument("-b", "--bam-file", type=str, required=True,
        help="processed bam file for peak-calling.")
    group_req_optz.add_argument("-g", "--genome-assembly", type=str, required=True,
        help="prefix used in 'buildidx' command for genome-indexing.")
    subparser_optz.add_argument("-ce", "--caller-executable", type=str, required=True,
        help="path of peak caller executable, e.g., /usr/bin/macs2")
    subparser_optz.add_argument("-cf", "--caller-flags", type=str, required=True,
        help="flags of peak caller")

    group_opt_optz = subparser_optz.add_argument_group('optional arguments\nbase')
    group_opt_optz.add_argument("-U", "--auc-coff", type=float, default=0.7,
        help="threshold of AUC to optimize the peaks\n(default: 0.7)")

    # parser for common (evalute, optimize)
    for subparser_comm, group_opt in [(subparser_eval, group_opt_eval), (subparser_optz, group_opt_optz)]:

        # optional arguments 
        group_opt.add_argument("-l", "--split-n", type=int, default=5000,
            help="number of peaks in a subset\n(default: 5000)")
        group_opt.add_argument("-o", "--score-col", type=int, default=8,
            help="col number by which be sorted\n(default: 8 [-log10P for MACS2])")
        group_opt.add_argument("-w", "--window-bp", type=int, default=600,
            help="size of scanning window\n(default: 600bp)")
        group_opt.add_argument("-mg", "--marginal-gc", type=float, default=0.02,
            help="GC errors allowed in generating null-seq\n(default=0.02)")
        group_opt.add_argument("-mr", "--marginal-rp", type=float, default=0.02,
            help="Repeat errors allowed in generating null-seq\n(default=0.02)")
        group_opt.add_argument("-s", "--random-seeds", type=int, default=-1,
            help="random seed number\nfor reproducibility\n(default: no seed)")
        group_opt.add_argument("-@", "--n-processes", type=int, default=10,
            help="number of processes\n(default: 10)")
        group_opt.add_argument("-v", "--verbosity", type=int, default=1,
            help="verbosity\n(default: 1), 0: silent")

        ## gkm-kernel options
        group_gkm = subparser_comm.add_argument_group('gkm-kernel')
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
        group_gkm.add_argument("-P", "--gkmsvm-mpi", type=str, default="none",
            help="Parallel computation of gkm-SVM using job schedular\n(default: none) options: slurm")

        ## svm training options
        group_svm = subparser_comm.add_argument_group('SVM training')
        group_svm.add_argument("-C", "--regularization", type=float, default=1.0,
            help="regularization parameter C\n(default: 1.0)")
        group_svm.add_argument("-e", "--precision", type=float, default=0.001,
            help="set the precision parameter epsilon\n(default: 0.001)")
        group_svm.add_argument("-u", "--shrinking", type=bool, default=False,
            help="if set, use the shrinking heuristics\n(default: False)")
        group_svm.add_argument("-c", "--cache-size", type=int, default=512,
            help="cache memory size in MB\n(default: 512MB)")
        group_svm.add_argument("-x", "--ncv", type=int, default=5,
            help="x-fold cross validation\nfor estimating effects of tags in training set\n(default: 5)")
        group_svm.add_argument("-r", "--repeats", type=int, default=1,
            help="number of repeats of CV training\nto reduce random variation\n(default: 1)")
        group_svm.add_argument("-f", "--fast-estimation", type=bool, default=False,
            help="fast estimation of AUC without nCV:\nusing nu score from trained SVM\n(default: False)")
        
    args = parser.parse_args()

    # formatting compatible with clog
    logfmt_str = '%(levelname)s %(asctime)s: %(message)s'
    datefmt_str = '%Y-%m-%d %H:%M:%S'

    logging.basicConfig(stream=sys.stdout,
        format=logfmt_str, datefmt=datefmt_str,
        level=logging.INFO)

    HEADER += "\n# Command line:" + ' '.join(sys.argv)
    HEADER += "\n# Parameters:"

    if args.commands == "buildidx":
        HEADER += "\n#   CHROM_FILE={0}".format(args.chrom_file)
        HEADER += "\n#   GENOME_ASSEMBLY={0}".format(args.genome_assembly)
        HEADER += "\n#   WINDOW_BP={0}".format(args.window_bp)
        HEADER += "\n#   N_PROCESSES={0}".format(args.n_processes)

    if args.commands == "evaluate":
        HEADER += "\n#   PEAK_FILE={0}".format(args.peak_file)
        HEADER += "\n#   EXP_NAME={0}".format(args.name)
        HEADER += "\n#   GENOME_ASSEMBLY={0}".format(args.genome_assembly)
        HEADER += "\n#   RANK_START={0}".format(args.rank_start)
        HEADER += "\n#   RANK_END={0}".format(args.rank_end)

    if args.commands == "optimize":
        HEADER += "\n#   GKMQC_FILE={0}".format(args.gkmqc_file)
        HEADER += "\n#   EXP_NAME={0}".format(args.name)
        HEADER += "\n#   BAM_FILE={0}".format(args.bam_file)
        HEADER += "\n#   GENOME_ASSEMBLY={0}".format(args.genome_assembly)
        HEADER += "\n#   CALLER_EXECUTABLE={0}".format(args.caller_executable)
        HEADER += "\n#   CALLER_FLAGS={0}".format(args.caller_flags)
        HEADER += "\n#   AUC_CUTOFF={0}".format(args.auc_coff)
    
    if args.commands == "evaluate" or args.commands == "optimize":
        # base
        HEADER += "\n#   N_PEAKS_PER_SUBSET={0}".format(args.split_n)
        HEADER += "\n#   SCORE_COL={0}".format(args.score_col)
        HEADER += "\n#   WINDOW_BP={0}".format(args.window_bp)
        HEADER += "\n#   MARGINAL_GC={0}".format(args.marginal_gc)
        HEADER += "\n#   MARGINAL_RP={0}".format(args.marginal_rp)
        HEADER += "\n#   RANDOM_SEEDS={0}".format(args.random_seeds)
        HEADER += "\n#   N_PROCESSES={0}".format(args.n_processes)
        HEADER += "\n#   VERBOSITY={0}".format(args.verbosity)

        # gkm-kernel
        HEADER += "\n#   KERNEL_TYPE={0}".format(args.kernel_type)
        HEADER += "\n#   FULL_WORD_LENGTH (L)={0}".format(args.full_word_length)
        HEADER += "\n#   NON_GAP_LENGTH   (k)={0}".format(args.non_gap_length)
        HEADER += "\n#   MAX_NUM_GAPS     (d)={0}".format(args.max_num_gaps)
        HEADER += "\n#   INIT_VAL_DECAY_F (M)={0}".format(args.init_decay)
        HEADER += "\n#   HALF_LIF_DECAY_F (H)={0}".format(args.half_life_decay)
        HEADER += "\n#   RBF_GAMMA        (G)={0}".format(args.rbf_gamma)
        HEADER += "\n#   GKMSVM_MPI={0}".format(args.gkmsvm_mpi)

        # svm training
        HEADER += "\n#   REGULARIZATION={0}".format(args.regularization)
        HEADER += "\n#   PRECISION={0}".format(args.precision)
        HEADER += "\n#   SHRINKING={0}".format(args.shrinking)
        HEADER += "\n#   CACHE_SIZE={0}".format(args.cache_size)
        HEADER += "\n#   CV={0}".format(args.ncv)
        HEADER += "\n#   CV_REPEATS={0}".format(args.repeats)
        HEADER += "\n#   FAST_ESTIMATE={0}".format(args.fast_estimation)

    logging.info(HEADER) 

    # output files
    import seqs_nullgen

    if args.commands == "buildidx":
        logging.info("### build null seq index")
        args_nidx = [args.chrom_file, args.genome_assembly, args.window_bp, args.n_processes]
        seqs_nullgen.build_nullseq_index(args_nidx)

    if args.commands == "evaluate":
        logging.info("### executing evaluate pipeline")
        print("## allocate output directory and make it if no exist.")
        gkmqc_out_dir = os.path.join(os.path.dirname(args.peak_file), args.name + ".gkmqc")

        if not os.path.isdir(gkmqc_out_dir):
            os.mkdir(gkmqc_out_dir)
        
        # Copy peak file
        os.system("cp %s %s" % (args.peak_file, gkmqc_out_dir))

        # Change Dir
        curdir = os.path.abspath('.')
        os.chdir(gkmqc_out_dir)
        
        # preprocess peaks
        import preprocess

        print("## QC and make a positive set")
        preprocess.make_qc_posset(gkmqc_out_dir, args)

        print("## split the positive set by p-value")
        ntests = preprocess.split_posset(gkmqc_out_dir, args)

        if args.rank_start > ntests:
            print("error: invalid range of ranks")
            sys.exit()
        
        if args.rank_end > ntests:
            args.rank_end = ntests

        print("## generate negative sets")
        pos_bed_files, neg_bed_files = preprocess.make_negset(gkmqc_out_dir, args)

        ## fa files
        pos_fa_files = list(map(lambda x: x.replace('.bed', '.fa'), pos_bed_files))
        neg_fa_files = list(map(lambda x: x.replace('.bed', '.fa'), neg_bed_files))

        # retrieve fa
        select_seqs_exe = os.path.join(dir_scripts, "seqs_select.py")
        for bed in pos_bed_files:
            fa = bed.replace('.bed', '.fa')
            if os.path.isfile(fa):
                continue
            os.system("python %s %s.e%d.qc.fa %s >%s" %\
                (select_seqs_exe, args.name, args.window_bp / 2, bed, fa))

        genome_fa_dir = os.path.join(dir_data, args.genome_assembly, 'fa')
        for bed in neg_bed_files:
            fa = bed.replace('.bed', '.fa')
            if os.path.isfile(fa):
                continue
            os.system("python %s/seqs_fetch.py -d %s %s %s" %\
                (dir_scripts, genome_fa_dir, bed, fa))
    
        ## cross-validate with gkmSVM
        # without job schedular
        print("## cross-validation with gkm-SVM:", end=' ')
        if args.gkmsvm_mpi == 'none':
            print("no job schedular mode")
            import gkmsvm
            for pos_fa, neg_fa in zip(pos_fa_files, neg_fa_files):
                print("cv: %s vs %s" % (pos_fa, neg_fa))
                gkmsvm.init(pos_fa, neg_fa, args)

        # job schedular: slurm
        elif args.gkmsvm_mpi == 'slurm':
            print("slurm")
            sbatch_exe = os.path.join(dir_scripts, "gkmsvm_slurm.sh")
            gkmsvm_py  = os.path.join(dir_scripts, "gkmsvm.py")
            gkmsvm_args = [
                "-w", "-s", "-@", "-v",
                "-t", "-L", "-k", "-d", "-M", "-H", "-G",
                "-C", "-e", "-u", "-c", "-x", "-r", "-f"
            ]
            gkmsvm_vals = [
                args.name, args.random_seeds, args.n_processes, args.verbosity,
                args.kernel_type, args.full_word_length, args.non_gap_length,
                args.max_num_gaps, args.init_decay, args.half_life_decay,
                args.rbf_gamma, args.regularization, args.precision, args.shrinking,
                args.cache_size, args.ncv, args.repeats, args.fast_estimation   
            ]
            gkmsvm_vals = list(map(str, gkmsvm_vals))
            args_vals_pairs = map(lambda x: ' '.join(x), list(zip(gkmsvm_args, gkmsvm_vals)))
            argc = ' '.join(list(args_vals_pairs))
            for pos_fa, neg_fa in zip(pos_fa_files, neg_fa_files):
                os.system("sbatch --export=NONE %s %s -p %s -n %s %s" %\
                (sbatch_exe, gkmsvm_py, pos_fa, neg_fa, argc))
        else:
            print("no available option for the job schedular")
            sys.exit()
        
        # return dir
        os.chdir(curdir)

    if args.commands == "optimize":
        # TODO: to be implemented
        logging.info("### score sequences using the trained model")
        print("TODO: under implementing")

if __name__=='__main__':
    main()