import os
import numpy as np
import subprocess

def optimize_peaks(prefix, base_dir, rc_suffix="_rc", t=600):
    # files
    out_file = os.path.join(base_dir, "%s.gkmqc/%s.e%d.optz.bed" % (prefix, prefix, t/2))
    file_gqc = os.path.join(base_dir, "%s.gkmqc/%s.gkmqc.eval.out" % (prefix, prefix))
    file_bed = os.path.join(base_dir, "%s.gkmqc/%s.e%d.bed" % (prefix, prefix, t/2))
    
    # check last AUC
    f = open(file_gqc)
    l_auc = min(map(lambda x: float(x.split()[3]), f.readlines()))
    f.close()
    
    # if least AUC >0.75, then start optimization
    n_ori_flag = False
    if l_auc > 0.75:
        prefix += rc_suffix
        file_gqc = os.path.join(base_dir, "%s.gkmqc/%s.gkmqc.eval.out" % (prefix, prefix))
        file_bed = os.path.join(base_dir, "%s.gkmqc/%s.e%d.bed" % (prefix, prefix, t/2))
        f = open(file_gqc)
        l_auc_opt = min(map(lambda x: float(x.split()[3]), f.readlines()))
        f.close()
        
        # if least AUC >0.7, then use optimized peak-calling file
        if l_auc_opt > 0.7:
            n_ori_flag = True
    
    # 0.7 < AUC < 0.75: use original peak-calling
    elif l_auc > 0.7:
        n_ori_flag = True
    
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
        return i

    else:
        f = open(file_gqc)
        ex_rank = 100
        i = 0
        for line in f.readlines():
            pf, _, _, auc_score, _ = line.split()
            rank = int(pf.split('.')[-2][3:])
            if float(auc_score) < 0.7 and rank < ex_rank:
                ex_rank = rank
            i += 1
        f.close()
        ex_rank -= 1
        
        # save file
        # get least score
        file_eps = os.path.join(base_dir, "%s.gkmqc/%s.e%d.qc.top%d.bed" % (prefix, prefix, t/2, ex_rank))
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
        return i