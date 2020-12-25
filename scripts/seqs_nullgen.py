"""
    nullseq.py: fast sampling of positive and negative sequence
    of open-chromain peaks

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

import os, sys
import tarfile, zipfile, gzip
import optparse
import zipfile, tarfile
import random

from bitarray import bitarray
from pyfasta import Fasta
import numpy as np
from multiprocessing import Pool
from collections.abc import Iterable
import logging

dir_this   = os.path.dirname(os.path.abspath(__file__))
dir_prnt   = os.path.dirname(dir_this)
base_data_dir = os.path.join(dir_prnt, "data")

##
# load .bit files (bit array - gc, rp, na)
##
def bitarray_fromfile(fn):
    try:
        fh = open(fn, 'rb')
        bits = bitarray()
        bits.fromfile(fh)
    except IOError as err:
        logging.error("I/O error: %s", err)
        sys.exit(0)
    return bits, fh

def per_chrom_idx_bits(fn, prefix_dir, chr):

    # check preexisting file
    fa_fn = os.path.join(prefix_dir, 'fa', os.path.basename(fn))
    if not os.path.isfile(fa_fn):
        logging.info("move fa to idx dir: %s", fa_fn)
        os.rename(fn, fa_fn)
    try:
        # Open chromosome fa
        logging.info("read fa file: %s", chr)
        f = open(fa_fn, "r")
    except IOError as err:
        logging.error("I/O error: %s", err)
        sys.exit(0)
    
    seq = ''.join(f.readlines()[1:]).strip()
    seq = seq.replace('\n', '')
    f.close()

    # N/gc/repeat bit array
    wchar_l = ["nN", "cgCG", "acgt"]
    barr_n  = ["na", "cg", "rp"]
    
    arr_list = []
    for wchar, bn in zip(wchar_l, barr_n):
        bit_fn  = os.path.join(prefix_dir, 'bit/%s.%s.bit' % (chr, bn))
        if not os.path.isfile(bit_fn):
            logging.info("no wchar bit array for %s: making wcard bit index..", chr)
            arr = bitarray(map(lambda c: c in wchar, seq))

            # save bit file
            fo = open(bit_fn, 'wb')
            arr.tofile(fo)
            fo.close()
        else:
            logging.info("wchar bit array for %s already exists. loading ..", chr)
            arr, f = bitarray_fromfile(bit_fn)
            f.close()
        arr_list.append(arr)

    del seq
    return arr_list

##
#  Sliding window with t and
#  calculate the count of repeats and CGs
##

def flatten(l):
    for el in l:
        if isinstance(el, Iterable) and not isinstance(el, (str, bytes)):
            yield from flatten(el)
        else:
           yield el

def per_chrom_nidx_l(fn, prefix_dir, chr, t, arr_list):

    nidx_pos_fn = os.path.join(prefix_dir, 'nidx_t%d/%s_pos.npy' % (t, chr))
    nidx_ptr_fn = os.path.join(prefix_dir, 'nidx_t%d/%s_ptr.npz' % (t, chr))

    if not (os.path.isfile(nidx_pos_fn) and os.path.isfile(nidx_ptr_fn)):
    
        # 3D-list for build index
        nidx_l = [[[] for col in range(t+1)] for row in range(t+1)] 

        # N/gc/repeat bit array
        na_arr, cg_arr, rp_arr = arr_list
        
        # scan with t-bp window
        # sliding and calculate the cnt of wildcard char
        na_cnt_c = na_arr[:t].count(True)
        cg_cnt_c = cg_arr[:t].count(True)
        rp_cnt_c = rp_arr[:t].count(True)

        logging.info("making nulls index for %s", fn)
        for i in range(0, len(na_arr) - t):
            # indexing only for regions with N-count = 0
            if not na_cnt_c:
                nidx_l[cg_cnt_c][rp_cnt_c].append(i) # start from 0

            # to minimize exhausitive count of trues
            na_cnt_c += (int(na_arr[i + t]) - int(na_arr[i]))
            cg_cnt_c += (int(cg_arr[i + t]) - int(cg_arr[i]))
            rp_cnt_c += (int(rp_arr[i + t]) - int(rp_arr[i]))
        
        n = 0
        nidx_ptr = np.ones(shape=(t+1, t+1), dtype=np.int32)
        for i, nl in enumerate(nidx_l):
            for j, el in enumerate(nl):
                nidx_ptr[i][j] = n
                n += len(el)

        np.savez_compressed(nidx_ptr_fn, ptr=nidx_ptr, len=n)
        del nidx_ptr
        
        nidx_pos = np.fromiter(flatten(nidx_l), count=-1, dtype=np.int32)
        del nidx_l
        
        # save nidx mat
        np.save(nidx_pos_fn, nidx_pos)
        del nidx_pos_fn

    else:
        logging.info("already have nidx_pos/ptr matrices for %s, skip.", fn)

##
# process function
##
def _proc_build_idx_func(fn, t, prefix_dir):

    chr = '.'.join(os.path.basename(fn).split('.')[:-1])
    barr_list = per_chrom_idx_bits(fn, prefix_dir, chr)
    per_chrom_nidx_l(fn, prefix_dir, chr, t, barr_list)
    
    return 0 # dummy output

def pool_wrapper_nidx_build(args):
    return _proc_build_idx_func(*args)

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

    if not os.path.isdir(prefix_dir):
        os.mkdir(prefix_dir)

    fseq_dir = prefix_dir + '/fa'
    if not os.path.isdir(fseq_dir):
        os.mkdir(fseq_dir)

    nidx_dir = prefix_dir + '/nidx_t%d' % t
    if not os.path.isdir(nidx_dir):
        os.mkdir(nidx_dir)

    barr_dir = prefix_dir + '/bit'
    if not os.path.isdir(barr_dir):
        os.mkdir(barr_dir)
    
    args_l = []
    logging.info("extract genome fa files...: %s", chrom_file)
    if zipfile.is_zipfile(chrom_file):
        zipfileobj = zipfile.ZipFile(chrom_file)
        for fn in zipfileobj.namelist():
            fn_dest = os.path.join(fseq_dir, os.path.basename(fn))
            if not os.path.isfile(fn_dest): # check preexisting fa file
                logging.info("no file in idx dir: extracting file : %s", fn) 
                zipfileobj.extract(fn)
            args_l.append((fn, t, prefix_dir))
        zipfileobj.close()

    elif tarfile.is_tarfile(chrom_file):
        tarfileobj = tarfile.open(chrom_file)
        for tarinfo in tarfileobj:
            if not tarinfo.isdir():
                fn = tarinfo.name
                fn_dest = os.path.join(fseq_dir, os.path.basename(fn))
                if not os.path.isfile(fn_dest): # check preexisting fa file
                    logging.info("no file in idx dir: extracting file : %s", fn) 
                    tarfileobj.extract(fn)
                args_l.append((fn, t, prefix_dir))
        tarfileobj.close()

    #elif chrom_file[-2:] == "gz":   # maybe not necessary in gkmQC
    #    fp = gzip.open(chrom_file, 'rt')
    #    q_tasks.put((fn, t, prefix_dir)))
    #    fp.close()
    else:
        logging.error("error: needs .zip or .tar.gz file")
        return 1
    
    # multiprocessing
    pool = Pool(p)
    pool.map(pool_wrapper_nidx_build, args_l)

    return 0

## TEST CODE
# build_nullseq_index(["hg38.chromFa.tar.gz", "hg38", 600, 10])

##
# read bed files
##
def read_bed_file(fn):
    try:
        f = open(fn, 'r')
        lines = f.readlines()
        f.close()
    except IOError as err:
        logging.error("I/O error: %s", err)
        sys.exit(0)
    
    posi_dic = {}
    for line in lines:
        if line[0] == '#':
            continue
        l = line.split()
        if not l[0] in posi_dic:
            posi_dic[l[0]] = []
        posi_dic[l[0]].append(int(l[1]))

    return posi_dic

##
# per-chromosome retrive of null-seq location
##
def _per_chrom_sample_nullseq_idx(pos_posi_l, genome, chrom, t, p, fold, gc_margin, rp_margin):
    
    # load bit array
    #print("loading bit array: gc, np, na")
    idxf_gc = os.path.join(base_data_dir, '%s/bit/%s' % (genome, '.'.join([chrom, 'cg', 'bit'])))
    idxf_rp = os.path.join(base_data_dir, '%s/bit/%s' % (genome, '.'.join([chrom, 'rp', 'bit'])))
    idxf_na = os.path.join(base_data_dir, '%s/bit/%s' % (genome, '.'.join([chrom, 'na', 'bit'])))
    gc_arr, _ = bitarray_fromfile(idxf_gc)
    rp_arr, _ = bitarray_fromfile(idxf_rp)
    na_arr, _ = bitarray_fromfile(idxf_na)

    # load nullseq_idx pos and ptr mat
    try:
        nidx_ptr_fn = os.path.join(base_data_dir, '%s/nidx_t%d/%s_ptr.npz' % (genome, t, chrom))
        nidx_ptr_d = np.load(nidx_ptr_fn)
        nidx_ptr = nidx_ptr_d['ptr']
        n = nidx_ptr_d['len']

        nidx_pos_fn = os.path.join(base_data_dir, '%s/nidx_t%d/%s_pos.npy' % (genome, t, chrom))
        nidx_pos = np.memmap(nidx_pos_fn, dtype="int32", mode="r", shape=(n,))

    except IOError as err:
        logging.error("I/O error: %s", err)
        sys.exit(0)

    sampled_posi_l = []
    for i, pos_posi in enumerate(pos_posi_l): # per peak subsets

        # check counted num of a set
        nidx_l_incr_ptr = np.full((t+1, t+1), 0, dtype=np.uint8)

        # mark-up NA for positive region
        na_arr_sub = na_arr.copy()

        for pos in pos_posi:
            na_arr_sub[pos:pos + t] = True

        sampled_posi = []

        # start search
        l_pos_posi = len(pos_posi)
        pos_i = 0
        eol_flag = False

        while len(sampled_posi) < l_pos_posi:
            if eol_flag:
                pos = random.choice(pos_posi)
            else:
                pos = pos_posi[pos_i]

            gc = gc_arr[pos:pos + t].count(True)
            rp = rp_arr[pos:pos + t].count(True)
            #print("#### seq %d: " % pos, end=' ')

            n_start = nidx_ptr[gc][rp]
            n_end = nidx_ptr[gc+int((rp+1)/(t+1))][(rp+1)%(t+1)]
            target_nidx = nidx_pos[n_start:n_end]
            target_ptr = nidx_l_incr_ptr[gc][rp]

            # random sampling of null seq with same gc/rp
            k = 0
            gc_d = 1 # margin
            rp_d = 1
            gc_i = 1 if random.random() < 0.5 else -1 # direction (plus, minus)
            rp_i = 1 if random.random() < 0.5 else -1
            ex_t = 1 if random.random() < 0.5 else -1  # 1: gc, -1: rp
            end_flag = False

            while k < fold:
                # nidx_length == 0 or scan finished -> change gc/rp
                while target_ptr == len(target_nidx):
                    if ex_t > 0: # gc
                        gc += gc_d * gc_i
                        gc_d += 1
                        gc_i *= -1
                        if gc_d > gc_margin:
                            end_flag = True
                            break

                    else: # rp
                        rp += rp_d * rp_i
                        rp_d += 1
                        rp_i *= -1
                        if rp_d > rp_margin:
                            end_flag = True
                            break

                    n_start = nidx_ptr[gc][rp]
                    n_end = nidx_ptr[gc+int((rp+1)/(t+1))][(rp+1)%(t+1)]
                    target_nidx = nidx_pos[n_start:n_end]
                    target_ptr = nidx_l_incr_ptr[gc][rp] # update
                    ex_t *= -1
                    #print ("ext,", end = ' ')

                if end_flag:
                    #print("none..")
                    break # cannot find no correspond neg sequence

                s = random.choice(target_nidx) # random choice
                if not na_arr_sub[s:s + t].any():
                    sampled_posi.append(s)
                    na_arr_sub[s:s + t] = True
                    k += 1
                    #print("found!")

                # update target ptr
                target_ptr += 1

            if not eol_flag:
                pos_i += 1

            if pos_i == l_pos_posi:
                eol_flag = True

        del na_arr_sub
        del nidx_l_incr_ptr

        sampled_posi_l.append((i, sampled_posi))
        logging.info("%s: finished %d-set!", chrom, i)
    
    del gc_arr, rp_arr, na_arr
    del nidx_ptr, nidx_pos
    del nidx_ptr_d

    return (chrom, sampled_posi_l)


def pool_wrapper_nidx_sample(args):
    return _per_chrom_sample_nullseq_idx(*args)

##
#  fetch genomic index of null (negative) sequences
##
def fetch_nullseq_beds(pos_bed_files, neg_bed_files, args_fetch_nb):
    
    genome = args_fetch_nb[0] # -b: genome prefix for index and bit array
    t      = args_fetch_nb[1] # -v: interval
    rseed  = args_fetch_nb[2] # -r: random number seed (default=1)
    p      = args_fetch_nb[3] # -T: # of processes
    fold   = 1 # -x: number of sequence to sample, FOLD times of given dataset (default=1)

    gc_margin = int(args_fetch_nb[4] * t) # -g: GC errors allowed (default=0.02) 
    rp_margin = int(args_fetch_nb[5] * t) # -t: RPT errors allowed (default=0.02) 
    
    # set random seed
    if rseed >= 0:
        random.seed(rseed)

    # read positions of positive peaks
    pos_posi_l = []
    for pos_bed_file in pos_bed_files:
        pos_posi_l.append(read_bed_file(pos_bed_file))

    # chromosomes
    chrnames = set()
    for pos_posi in pos_posi_l:
        chrnames.update(set(pos_posi.keys()))
    chrnames = sorted(list(chrnames))

    positive_l = [] # by chr

    args_l = []
    for chrom in chrnames:
        pos_posi_l_by_chr = []
        for pos_posi_dic in pos_posi_l: 
            pos_posi_l_by_chr.append(pos_posi_dic.get(chrom, []))
        args_l.append((pos_posi_l_by_chr, genome, chrom, t, p, fold, gc_margin, rp_margin))
        positive_l.append(pos_posi_l_by_chr)

    # multiprocessing
    pool = Pool(p)
    results_l = pool.map(pool_wrapper_nidx_sample, args_l)

    #n = len(chrnames)
    results_l.sort(key=lambda x: x[0])

    # write bed files
    fo_l = list(map(lambda x: open(x, "w"), neg_bed_files))
    for chrom, neg_posi_l in results_l:
        for i, neg_posi in neg_posi_l:
            #print(ineg_posi[1])
            outstr_l = map(lambda x: "%s\t%d\t%d" % (chrom, x, x + t), sorted(neg_posi))
            fo_l[i].write('\n'.join(list(outstr_l)) + '\n')

    for fo in fo_l:
        fo.close()
    
    logging.info("fetch fasta seq")
    # write fa files (pos, neg)
    fa_files = list(map(lambda x: x.replace('.bed', '.fa'), pos_bed_files + neg_bed_files))
    fo_l = list(map(lambda x: open(x, "w"), fa_files))
    for pos_posi_l, (chrom, neg_posi_l) in zip(positive_l, results_l):
        logging.info(chrom)
        # load fasta object
        chr_fa = os.path.join(base_data_dir, '%s/fa/%s.fa' % (genome, chrom))
        f = Fasta(chr_fa)[chrom]

        # write fa to files
        for pos_posi, (i, neg_posi) in zip(pos_posi_l, neg_posi_l):
            for x in pos_posi:
                outstr = ">%s:%d-%d\n%s\n\n" % (chrom, x+1, x+t, f[x:x+t].upper())
                fo_l[i].write(outstr)

            for x in sorted(neg_posi):
                outstr = ">%s:%d-%d\n%s\n\n" % (chrom, x+1, x+t, f[x:x+t].upper())
                fo_l[i+len(pos_bed_files)].write(outstr)

    for fo in fo_l:
        fo.close()
