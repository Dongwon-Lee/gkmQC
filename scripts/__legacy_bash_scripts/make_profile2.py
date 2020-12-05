#!/usr/bin/env python
"""
    make_profile.py: generate sequence profile of the given data set

    Copyright (C) 2011 Dongwon Lee

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

import os
import os.path
import sys
import random
import optparse

from bitarray import bitarray


def read_pos_file(filename):
    """
    """
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    positions = []

    for line in lines:
        if line[0] == '#': 
            continue

        l = line.split()
        g = l[0].split(':')
        h = g[1].split('-')

        positions.append((g[0], int(h[0])-1, int(h[1]))) #correct to 0-based half open coordinate system

    return positions

def read_bed_file(filename):
    """
    """
    f = open(filename, 'r')
    lines = f.readlines()
    f.close()

    positions = []

    for line in lines:
        if line[0] == '#': 
            continue

        l = line.split()

        positions.append((l[0], int(l[1]), int(l[2])))

    return positions


def bitarray_fromfile(filename):
    """
    """
    fh = open(filename, 'rb')
    bits = bitarray()
    bits.fromfile(fh)

    return bits, fh


def get_seqid(buildname, pos):
    return '_'.join( [buildname, pos[0], str(pos[1]), str(pos[2]), '+'] )

def get_seqid_old(buildname, pos):
    return pos[0] + ':' + str(pos[1]+1) + '-' + str(pos[2])

def make_profile(positions, buildname, basedir):
    """
    """
    chrnames = sorted(set(map(lambda p: p[0], positions)))

    profiles = {}
    for chrom in chrnames:
        idxf_gc = os.path.join(basedir, '.'.join([buildname, chrom, 'gc', 'out']))
        idxf_rpt = os.path.join(basedir, '.'.join([buildname, chrom, 'rpt', 'out']))
        idxf_na = os.path.join(basedir, '.'.join([buildname, chrom, 'na', 'out']))

        #if os.path.exists(idxf_gc) == False or os.path.exists(idxf_rpt) == False:
        #   continue
        bits_gc, gcf = bitarray_fromfile(idxf_gc)
        bits_rpt, rptf = bitarray_fromfile(idxf_rpt)
        bits_na, naf = bitarray_fromfile(idxf_na)

        for pos in positions:
            if pos[0] != chrom:
                continue

            seqid = get_seqid_old(buildname, pos)
            slen = pos[2]-pos[1]
            gc = bits_gc[pos[1]:pos[2]].count(True)
            rpt = bits_rpt[pos[1]:pos[2]].count(True)
            na = bits_na[pos[1]:pos[2]].count(True)

            profiles[seqid] = (slen, gc, rpt, na)

        gcf.close()
        rptf.close()
        naf.close()

    return profiles

def main():
    if len(sys.argv) != 5:
        print("Usage:" + sys.argv[0] + " BEDFILE(or POSFILE) BUILDNAME BASE_DIR OUT_FILE")
        sys.exit()

    bedfile = sys.argv[1]
    buildname = sys.argv[2]
    basedir = sys.argv[3]
    output = sys.argv[4]

    if bedfile[-3:] == 'bed':
        positions = read_bed_file(bedfile)
    else:
        positions = read_pos_file(bedfile)

    seqids = []
    for pos in positions:
        seqids.append(get_seqid_old(buildname, pos))

    profiles = make_profile(positions, buildname, basedir)

    f = open(output, 'w')
    for seqid in seqids:
        prof = profiles[seqid]
        l = float(prof[0])
        f.write('\t'.join( map(str, [seqid, prof[0], prof[1]/l, prof[2]/l, prof[3]/l]) ) + '\n')
    f.close()

if __name__ == "__main__": main()
