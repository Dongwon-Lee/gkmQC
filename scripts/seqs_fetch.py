"""
    seqs_fetch.py: fast fetch of positive and negative sequence
    of open-chromain peaks

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
import tarfile, zipfile, gzip
import optparse
import zipfile, tarfile
import random

from bitarray import bitarray
from multiprocessing import Process, Queue

base_data_dir = os.path.abspath('%s/../data' % os.path.dirname(__file__))

##
#  Sliding window with t=600 (default) and
#  calculate the count of repeats and CGs
##
def per_chrom_index(fn, t):
    try:
        # Open chromosome fa
        f = open(fn, "r")
    except IOError as err:
        print("I/O error: ", err)
        sys.exit(0)
    
    seq = ''.join(f.readlines()[1:]).strip()
    f.close()

    # 3D-list for build index
    idx_li = [[[] for col in range(t)] for row in range(t)] 

    # N/gc/repeat bit array
    na_arr = bitarray(map(lambda c: c in "nN",   seq))
    cg_arr = bitarray(map(lambda c: c in "cgCG", seq))
    rp_arr = bitarray(map(lambda c: c in "acgt", seq)) 

    del seq
    
    # scan with t-bp window
    # sliding and calculate the cnt of wildcard char
    na_cnt_c = na_arr[:t].count(True)
    cg_cnt_c = cg_arr[:t].count(True)
    rp_cnt_c = rp_arr[:t].count(True)

    for i in range(0, len(na_arr) - t):
        # indexing only for regions with N-count = 0
        if not na_cnt_c:
            idx_li[cg_cnt_c][rp_cnt_c].append(i) 

        # to minimize exhausitive count of trues
        na_cnt_c += (int(na_arr[i + t]) - int(na_arr[i]))
        cg_cnt_c += (int(cg_arr[i + t]) - int(cg_arr[i]))
        rp_cnt_c += (int(rp_arr[i + t]) - int(rp_arr[i])) 

    os.remove(fn)
    idx_bits = [na_arr, cg_arr, rp_arr]

    return idx_bits, idx_li

##
#  Save index list to pickle file
##
def save_idx_file(idx_bits, idx_li, prefix_dir, fn, t):

    # get prefix from chromosome filename
    prefix = '.'.join(fn.split('.')[:-1]) 
    if not os.path.isdir(prefix_dir):
        os.mkdir(prefix_dir)

    # Save pre-calculated idx list to pkl ext
    pkl_fn   = prefix + '.pkl'
    pkl_path = "%s/nidx_t%d/%s" % (prefix_dir, t, pkl_fn)
    fo = open(pkl_path, "wb")
    pickle.dump(idx_li, fo)
    fo.close()

    bit_pref = ['na', 'cg', 'rp']
    # Save bitarray into file
    for arr, pre in zip(idx_bits, bit_pref):
        bit_fn   = prefix + '.%s.bit' % pre
        bit_path = "%s/%s" % (prefix_dir, bit_fn)
        fo = open(bit_path, 'wb')
        arr.tofile(fo)
        fo.close()

##
# process function
##
def _proc_build_idx_func(fn, t, prefix_dir):
    idx_bits, idx_li = per_chrom_index(fn, t)
    save_idx_file(idx_bits, idx_li, prefix_dir, fn, t)
    return 0 # dummy output

##
# process functor: return worker function
# that receives args from task queue
# and does allocated function
##  
def worker_funct(_p_func):
    def worker(q_tasks, q_results):
        args = q_tasks.get()
        r = _p_func(*args)
        q_results.put(r)
    return worker

##
#  Build genome index of null sequences and
#  saves to pickle object
#  modified from Dongwon's code (nullseq_build_indices.py)
##
def build_nullseq_index(args_nidx):
    chrom_file = args_nidx[0]  # chromosome file (archived)
    prefix = args_nidx[1]      # name of genome assembly (e.g., hg38)
    t = args_nidx[2]           # size of window (e.g., 600bp)
    p = args_nidx[3]           # processes (< # of total chromosomes)
    prefix_dir = '%s/%s' % (base_data_dir, prefix)
    q_tasks = Queue()
    q_results = Queue() # dummy queue
    
    if zipfile.is_zipfile(chrom_file):
        zipfileobj = zipfile.ZipFile(chrom_file)
        for fn in zipfileobj.namelist():
            zipfileobj.extract(fn)
            q_tasks.put((fn, t, prefix_dir))      # put process into worker queue
        zipfileobj.close()

    elif tarfile.is_tarfile(chrom_file):
        tarfileobj = tarfile.open(chrom_file)
        for tarinfo in tarfileobj:
            if not tarinfo.isdir():
                fn = tarinfo.name
                tarfileobj.extract(fn)
                q_tasks.put((fn, t, prefix_dir))  # put process into worker queue
        tarfileobj.close()

    #elif chrom_file[-2:] == "gz":   # maybe not necessary in gkmQC
    #    fp = gzip.open(chrom_file, 'rt')
    #    q_tasks.put((fn, t, prefix_dir)))
    #    fp.close()
    else:
        return 1
    
    # multiprocessing: allocate workers
    workers = [Process(target=worker_funct(_proc_build_idx_func),\
         args=(q_tasks, q_results)) for _ in range(p)]
    for each in workers:
        each.start()

    # wait finish of worker processes    
    for each in workers:
        each.join()

    return 0

##
# read bed files
##
def read_bed_file(fn):
    try:
        f = open(fn, 'r')
        lines = f.readlines()
        f.close()
    except IOError as err:
        print("I/O error: ", err)
        sys.exit(0)
    
    positions = []
    for line in lines:
        if line[0] == '#':
            continue
        l = line.split()
        positions.append((l[0], int(l[1]), int(l[2])))
    return positions

##
# load .bit files (bit array - gc, rp, na)
##
def bitarray_fromfile(fn):
    try:
        fh = open(fn, 'rb')
        bits = bitarray()
        bits.fromfile(fh)
    except IOError as err:
        print("I/O error: ", err)
        sys.exit(0)
    return bits, fh

##
# per-chromosome retrive of null-seq location
# rewrote from Dongwon's code (nullseq_generate.py)
##
def _per_chrom_sample_nullseq_idx(position_l, genome, chrom, t, fold):
    
    # load bit array
    idxf_gc = os.path.join(base_data_dir, genome, '.'.join([chrom, 'gc', 'bit']))
    idxf_rp = os.path.join(base_data_dir, genome, '.'.join([chrom, 'rp', 'bit']))
    gc_arr, _ = bitarray_fromfile(idxf_gc)
    rp_arr, _ = bitarray_fromfile(idxf_rp)

    # load nullseq_idx pkl
    nidxf = os.path.join(base_data_dir, '%s/t%d/%s.pkl' % (genome, t, chrom))

    try:
        f = open(nidxf, 'rb')
        idx_li = pickle.load(f)
        f.close()
    except IOError as err:
        print("I/O error: ", err)
        sys.exit(0)

    sampled_posi_l = [[] for _ in range(len(position_l))]
    for i, positions in enumerate(position_l): # per peak subsets

        # make profile (pos, neg)
        #profiles_ = []
        sampled_posi = sampled_posi_l[i]
        for pos in positions:
            if pos[0] != chrom:
                continue
            #seqid = pos[0] + ':' + str(pos[1]+1) + '-' + str(pos[2]) // legacy code: seqid
            gc = gc_arr[pos[1]:pos[2]].count(True)
            rp = rp_arr[pos[1]:pos[2]].count(True)

            # random sampling of null seq with same gc/rp
            k = 0
            while k < fold:
                i = random.sample(idx_li[gc][rp], 1)[0]
                if not i in sampled_posi and i != pos[1]:
                    sampled_posi.append(i)
                    k += 1

        # this bit array is used to mark positions that are excluded from sampling
        # this will be updated as we sample more sequences in order to prevent sampled sequences from overlapping
        #idxf_na = os.path.join(base_data_dir, '.'.join([genome, chrom, 'na', 'bit']))
        #na_arr, naf = bitarray_fromfile(idxf_na)
        #naf.close()

    return (chrom, sampled_posi_l)

##
#  fetch genomic index of null (negative) sequences
#  using index file, rewrote from Dongwon's code (nullseq_generate.py)
##
def fetch_nullseq_beds(pos_bed_files, neg_bed_files, args_fetch_nb):
    
    genome = args_fetch_nb[0] # -b: genome prefix for index and bit array
    t      = args_fetch_nb[1] # -v: interval
    fold   = args_fetch_nb[2] # -x: number of sequence to sample, FOLD times of given dataset (default=1)
    rseed  = args_fetch_nb[3] # -r: random number seed (default=1)
    p      = args_fetch_nb[4] # -T: # of processes

    #gc_err = args_fetch_nb[5] # -g: GC errors allowed (default=0.02)   - no longer use in gkmQC
    #rp_err = args_fetch_nb[6] # -t: RPT errors allowed (default=0.02)  - no longer use in gkmQC
    #no_gc  = args_fetch_nb[7] # -G: do not match gc-contents  - no longer use in gkmQC
    #no_rp  = args_fetch_nb[8] # -R: do not match repeats - no longer use in gkmQC
    #fn_out = args_fetch_nb[9] # -o: set the name of output file (default=nullseq_output.bed) - processed from pos_bed_f
    
    # set random seed
    rseed = args_fetch_nb[4]
    random.seed(rseed)

    # args_fcts_null_bed
    # read positions of positive peaks
    
    position_l = []
    for pos_bed_file in pos_bed_files:
        position_l.append(read_bed_file(pos_bed_file))

    # chromosomes
    chrnames = sorted(set(map(lambda p: p[0], position_l[0])))

    q_tasks = Queue()
    q_results = Queue() # dummy queue
    for chrom in chrnames:
        position_l_by_chr = []
        for positions in position_l: 
            position_l_by_chr.append(list(filter(lambda x: x[0] == chrom, positions)))
        q_tasks.put((position_l_by_chr, genome, chrom, t, fold))

    # multiprocessing: allocate workers
    workers = [Process(target=worker_funct(_per_chrom_sample_nullseq_idx),\
         args=(q_tasks, q_results)) for _ in range(p)]
    for each in workers:
        each.start()

    # wait finish of worker processes    
    for each in workers:
        each.join()
    
    # 

        
    return

##
# Fetch sequence of positive and negative (null) sequences  
##
def fetch_seqs_fa(args_fcts_fa):
    # args_fcts_pos_neg_seq
    return

