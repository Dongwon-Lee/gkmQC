"""
    preprocess.py: treat and quality control of peak sequences

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

import os, random
from bitarray import bitarray
import seqs_nullgen
from seqs_nullgen import bitarray_fromfile
import logging

dir_this   = os.path.dirname(os.path.abspath(__file__))
dir_prnt   = os.path.dirname(dir_this)
base_data_dir = os.path.join(dir_prnt, "data")
dir_scripts = os.path.join(dir_prnt, "scripts")

def make_profile(bed_file, prof_file, genome_assembly):
    # load bit arraies
    arr_cg_dic = {}
    arr_rp_dic = {}
    arr_na_dic = {}

    def get_bit_array(genome_assembly, chr, pr):
        bit_dir = os.path.join(os.path.join(base_data_dir, genome_assembly, "bit"))
        return bitarray_fromfile(os.path.join(bit_dir, "%s.%s.bit" % (chr, pr)))[0]

    f = open(bed_file)
    fo = open(prof_file, "w")

    for line in f.readlines():
        line_tab = line.split()
        chr     = line_tab[0]
        start   = int(line_tab[1])
        end     = int(line_tab[2])
        seq_len = end - start
        seq_id  = '%s:%d-%d' % (chr, start+1, end)

        # gc
        if not chr in arr_cg_dic:
            arr_cg_dic[chr] = get_bit_array(genome_assembly, chr, "cg")
        cg = arr_cg_dic[chr][start:end].count(True) / seq_len
            
        # rp
        if not chr in arr_rp_dic:
            arr_rp_dic[chr] = get_bit_array(genome_assembly, chr, "rp")
        rp = arr_rp_dic[chr][start:end].count(True) / seq_len

        # na
        if not chr in arr_na_dic:
            arr_na_dic[chr] = get_bit_array(genome_assembly, chr, "na")
        na = arr_na_dic[chr][start:end].count(True) / seq_len
        
        # write
        fo.write('\t'.join(map(str, [seq_id, seq_len, cg, rp, na])) + '\n')

    fo.close()
    f.close()

##
# QC and make a positive set       
# score_col: 7 - hotspot2, 8 - macs2
##
def make_qc_posset(gkmqc_out_dir, args):

    peak_file = args.peak_file
    prefix = args.name
    window_bp = args.window_bp
    genome_assembly = args.genome_assembly
    score_col = args.score_col

    # FILES
    ext_len = window_bp / 2
    prefix = "%s.e%d" % (prefix, ext_len)
    posf0 = "%s.bed" % prefix
    posf0_prof = "%s.prof" % prefix
    posf = "%s.qc.bed" % prefix

    # 1. make fixed length peaks
    logging.info("make fixed length peaks")
    if os.path.isfile(posf0):
        logging.info("skip making %s", posf0)
    else:
        os.system("awk -v OFS='\t' -v SHFT=%d '$1 ~ /^chr[0-9XY]+$/ && $2+$10 > SHFT {\
            summit=$2+$10; print $1,summit-SHFT,summit+SHFT,$4,$%d}' %s >%s" %\
            (ext_len, score_col, peak_file, posf0))
    
    # 2. calculate gc/rp/na profiles of the fixed length peaks
    logging.info("calculate gc/rp/na profiles of the fixed length peaks")
    skip_flag = False
    if os.path.isfile(posf0_prof):
        nb = sum(1 for line in open(posf0))
        np = sum(1 for line in open(posf0_prof))
        if nb == np:
            logging.info("skip making %s", posf0_prof)
            skip_flag = True
    if not skip_flag:
        make_profile(posf0, posf0_prof, genome_assembly)
    

    logging.info("remove peaks with >1% of N bases & >70% of repeats")
    # 3. remove peaks with >1% of N bases & >70% of repeats
    if os.path.isfile(posf):
        logging.info("skip making %s", posf)
    else:
        os.system("paste %s %s | awk '$4<=0.7 && $5<=0.01' |cut -f 6- >%s" %\
            (posf0_prof, posf0, posf))

##
# split the positive set by p-value 
##
def split_posset(gkmqc_out_dir, args):

    prefix = args.name
    window_bp = args.window_bp
    split_n = args.split_n

    # FILES
    ext_len = window_bp / 2
    prefix = "%s.e%d" % (prefix, ext_len)
    posf   = "%s.qc.bed" % prefix

    # read
    ntot = 0
    posf_l = []
    f = open(posf)
    for line in f.readlines():
        ch, s, e, sid, score = line.split()
        posf_l.append((ch, int(s), int(e), sid, float(score)))
        ntot += 1
    f.close()

    ntests = int((ntot + int(split_n / 2)) / split_n)

    # Sort
    logging.info("sort peaks")
    posf_l.sort(key=lambda x: x[4], reverse=True)
    posf_lr = []
    prev_score = posf_l[0][4]
    prev_argi  = 0
    for i, posf_e in enumerate(posf_l):
        if posf_e[4] != prev_score or i == len(posf_l) - 1:
            sub = posf_l[prev_argi:i]
            if len(sub) > 1: random.shuffle(sub) # random shuffle peaks with same signal
            posf_lr += sub
            prev_score = posf_e[4]
            prev_argi = i
    del posf_l

    # Split
    logging.info("split processing")
    for i in range(ntests):
        s = split_n * i
        if i == ntests - 1: e = ntot
        else: e = split_n * (i + 1)
        
        fo = open("%s.top%d.bed" % (posf[:-4], i+1), "w")
        for line in sorted(posf_lr[s:e]):
            fo.write('\t'.join(map(str, line)) + '\n')
        fo.close()

    return ntests
##
# generate a negative set for each of the splitted positive sets (rank start and end)
##
def make_negset(gkmqc_out_dir, args):

    prefix = args.name
    window_bp = args.window_bp
    genome = args.genome_assembly
    rank_start = args.rank_start
    rank_end = args.rank_end
    rseed = args.random_seeds
    p = args.n_processes
    gc_margin = args.marginal_gc
    rp_margin = args.marginal_rp

    # FILES
    ext_len = window_bp / 2
    prefix = "%s.e%d" % (prefix, ext_len)

    pos_bed_files = list(map(lambda x: "%s.qc.top%d.bed" % (prefix, x), range(rank_start, rank_end + 1)))
    neg_bed_files = list(map(lambda x: "%s.qc.top%d.nr1.bed" % (prefix, x), range(rank_start, rank_end + 1)))

    pos_exist_l = len([f for f in pos_bed_files if os.path.isfile(f)])
    neg_exist_l = len([f for f in neg_bed_files if os.path.isfile(f)])

    if pos_exist_l == neg_exist_l:
        logging.info("skip making negative set")
    else:
        # generate null seq
        args_nseq_gen = [genome, window_bp, rseed, p, gc_margin, rp_margin]
        seqs_nullgen.fetch_nullseq_beds(pos_bed_files, neg_bed_files, args_nseq_gen)

    return (pos_bed_files, neg_bed_files)
    