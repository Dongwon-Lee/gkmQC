import os, glob
import numpy as np
import subprocess
import logging

def optimize_peaks(args):

    prefix = args.gkmqc_prefix
    prefix_rc = args.gkmqc_rt_prefix
    base_dir = args.base_dir

    # inference ext-num
    file_prof = glob.glob("%s/%s.gkmqc/%s.e*.prof" % (base_dir, prefix, prefix))[0]
    ext = int(file_prof.split('.')[-2][1:])

    # files
    out_file = os.path.join(base_dir, "%s.gkmqc/%s.e%d.optz.bed" % (prefix, prefix, ext))
    file_gqc = os.path.join(base_dir, "%s.gkmqc/%s.gkmqc.eval.out" % (prefix, prefix))
    file_bed = os.path.join(base_dir, "%s.gkmqc/%s.e%d.bed" % (prefix, prefix, ext))
    
    # check last AUC
    f = open(file_gqc)
    l_auc = min(map(lambda x: float(x.split()[3]), f.readlines()))
    f.close()
    
    # if least AUC >0.75, then start optimization
    n_ori_flag = False
    if l_auc > args.auc_start_opt:
        logging.info("least AUC = %.3f > %.2f: start optimizing peaks from relaxed threshold",\
            l_auc, args.auc_start_opt)
        prefix = prefix_rc
        file_gqc = os.path.join(base_dir, "%s.gkmqc/%s.gkmqc.eval.out" % (prefix, prefix))
        file_bed = os.path.join(base_dir, "%s.gkmqc/%s.e%d.bed" % (prefix, prefix, ext))
        f = open(file_gqc)
        l_auc_opt = min(map(lambda x: float(x.split()[3]), f.readlines()))
        f.close()
        
        # if least AUC >0.7, then use optimized peak-calling file
        if l_auc_opt > args.auc_min_coff:
            logging.info("%.2f < least AUC from recalled peaks = %.3f < %.2f: use all peaks from relaxed threshold",\
                args.auc_min_coff, l_auc_opt, args.auc_start_opt)
            n_ori_flag = True
        else:
            logging.info("least AUC from recalled peaks = %.3f < %.2f: filtering peaks with gkmQC AUC",\
                l_auc_opt, args.auc_min_coff)
    
    # 0.7 < AUC < 0.75: use original peak-calling
    elif l_auc > args.auc_min_coff:
        logging.info("%.2f < least AUC = %.3f < %.2f: use all original peaks",\
            args.auc_min_coff, l_auc, args.auc_start_opt)
        n_ori_flag = True
    
    else:
        logging.info("least AUC = %.3f < %.2f: filtering peaks with gkmQC AUC",\
            l_auc, args.auc_min_coff)
    
    if n_ori_flag:
        f = open(file_bed)
        fo = open(out_file, "w")
        i = 0
        for line in f.readlines():
            line_tab = line.split()
            if int(line_tab[1]) > 0:
                fo.write(line)
                i += 1
        fo.close()
        f.close()

    else:
        f = open(file_gqc)
        ex_rank = np.inf
        i = 0
        for line in f.readlines():
            pf, _, _, auc_score, _ = line.split()
            rank = int(pf.split('.')[-2][3:])
            if float(auc_score) < args.auc_min_coff and rank < ex_rank:
                ex_rank = rank
            i += 1
        f.close()
        ex_rank -= 1
        
        # save file
        # get least score
        file_eps = os.path.join(base_dir, "%s.gkmqc/%s.e%d.qc.top%d.bed" % (prefix, prefix, ext, ex_rank))
        l_sig_val = np.inf
        f = open(file_eps)
        for line in f.readlines():
            curr_score = float(line.split()[-1])
            if curr_score < l_sig_val:
                l_sig_val = curr_score
        f.close()
        
        # filter peaks
        f = open(file_bed)
        fo = open(out_file, "w")
        i = 0
        for line in f.readlines():
            line_tab = line.split()
            if float(line_tab[4]) >= l_sig_val and int(line_tab[1]) > 0:
                fo.write(line)
                i += 1
        fo.close()
        f.close()
    
    logging.info("Done. Total optimized peaks = %d", i)
    logging.info("Optimized peaks have been saved to:")
    logging.info("%s", out_file)