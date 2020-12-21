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

import os, subprocess
from bitarray import bitarray
import seqs_nullgen
from seqs_nullgen import bitarray_fromfile
from subprocess import PIPE

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
    posf_fasta = "%s.qc.fa" % prefix

    # 1. make fixed length peaks
    print("make fixed length peaks")
    if os.path.isfile(posf0):
        print("skip making %s" % posf0)
    else:        
        os.system("awk -v OFS='\t' -v SHFT=%d '$1 ~ /^chr[0-9XY]+$/ {\
            summit=$2+$10; print $1,summit-SHFT,summit+SHFT,$4,$%d}' %s >%s" %\
            (ext_len, score_col, peak_file, posf0)) 
    
    # 2. calculate gc/rp/na profiles of the fixed length peaks
    print("calculate gc/rp/na profiles of the fixed length peaks")
    skip_flag = False
    if os.path.isfile(posf0_prof):
        nb = int(subprocess.getoutput("cat %s | wc -l" % posf0))
        np = int(subprocess.getoutput("cat %s | wc -l" % posf0_prof))
        if nb == np:
            print("skip making %s" % posf0_prof)
            skip_flag = True
    if not skip_flag:
        make_profile(posf0, posf0_prof, genome_assembly)
    
    print("remove peaks with >1% of N bases & >70% of repeats")
    # 3. remove peaks with >1% of N bases & >70% of repeats
    if os.path.isfile(posf):
        print("skip making %s" % posf)
    else:
        os.system("paste %s %s | awk '$4<=0.7 && $5<=0.01' |cut -f 6- >%s" %\
            (posf0_prof, posf0, posf))
    
    # 4. make fasta file
    print("make fasta file")
    skip_flag = False
    if os.path.isfile(posf_fasta):
        n2 = 2 * int(subprocess.getoutput("cat %s | wc -l" % posf))
        s2 = int(subprocess.getoutput("cat %s | wc -l" % posf_fasta))
        if n2 == s2:
            print("skip making %s" % posf_fasta)
            skip_flag = True

    if not skip_flag:
        genome_fa_dir = os.path.join(base_data_dir, genome_assembly, 'fa')
        os.system("python %s/seqs_fetch.py -d %s %s %s" %\
            (dir_scripts, genome_fa_dir, posf, posf_fasta))

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
    ntot   = int(subprocess.getoutput("cat %s | wc -l" % posf))
    ntests = int((ntot + int(split_n / 2)) / split_n)

    # Sort
    print("sort peaks")
    os.system("sort -gr -k5,5 %s >%s.tmp.sorted" %\
        (posf, posf))

    # Split
    print("split processing...")
    
    for i in range(1, ntests + 1):
        skipn = (i - 1) * split_n + 1
        if i == ntests:
            subprocess.Popen("tail -n +%d %s.tmp.sorted | sortBed >%s.top%d.bed" %\
                (skipn, posf, posf[:-4], i), shell=True, stdin=PIPE, stdout=PIPE)
        else:
            subprocess.Popen("tail -n +%d %s.tmp.sorted | head -n %d | sortBed >%s.top%d.bed" %\
                (skipn, posf, split_n, posf[:-4], i), shell=True, stdin=PIPE, stdout=PIPE)
    # Remove
    os.system("rm -f %s.tmp.sorted" % posf)

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

    pos_bed_files = list(map(lambda x: "%s.qc.top%d.bed" % (prefix, x), range(rank_start, rank_end)))
    neg_bed_files = list(map(lambda x: "%s.qc.top%d.nr1.bed" % (prefix, x), range(rank_start, rank_end)))

    pos_exist_l = len([f for f in pos_bed_files if os.path.isfile(f)])
    neg_exist_l = len([f for f in neg_bed_files if os.path.isfile(f)])

    if pos_exist_l == neg_exist_l:
        print("skip making negative set")
    else:
        # generate null seq
        args_nseq_gen = [genome, window_bp, rseed, p, gc_margin, rp_margin]
        seqs_nullgen.fetch_nullseq_beds(pos_bed_files, neg_bed_files, args_nseq_gen)

    return (pos_bed_files, neg_bed_files)
    