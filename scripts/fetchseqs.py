#!/usr/bin/env python
"""
    fetchseqs.py: fetch sequences of the given positions from the genome

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

import sys
import os
import os.path
import optparse
import gzip
import pickle
import zlib


def read_bedfile(fn, seqid_col=0):
    f = open(fn, 'r')
    positions = []

    for line in f:
        if line[0] == "#": 
            continue

        if line.find('\t') > 0:
            l = line.strip().split('\t')
        else:
            l = line.strip().split(' ')

        chrom = l[0]
        pos_start = int(l[1])
        pos_end = int(l[2])

        seqid = ""
        if seqid_col > 0:
            seqid = l[seqid_col-1]

        positions.append((chrom, pos_start, pos_end, seqid))

    f.close()
    return positions


def read_posfile(fn, seqid_col=0):
    f = open(fn, 'r')
    positions = []

    for line in f:
        if line[0] == "#": 
            continue

        if line.find('\t') > 0:
            l = line.strip().split('\t')
        else:
            l = line.strip().split(' ')

        m = l[0].split(':')
        n = m[1].split('-')

        chrom = m[0]
        pos_start = int(n[0])-1 #1-based to 0-based
        pos_end = int(n[1])

        seqid = ""
        if seqid_col > 0:
            seqid = l[seqid_col-1]

        positions.append((chrom, pos_start, pos_end, seqid)) 
    
    f.close()
    return positions


def read_snpfile(fn, seqid_col=0):
    if fn[-2:] == 'gz':
        f = gzip.open(fn, 'r')
    else:
        f = open(fn, 'r')

    positions = []

    for line in f:
        if line[0] == "#": 
            continue

        if line.find('\t') > 0:
            l = line.strip().split('\t')
        else:
            l = line.strip().split(' ')

        if l[0] == "MT":
            chrom = "chrM"
        else:
            chrom = "chr" + l[0]

        snppos = int(l[1]) # 1-based
        len_refallele = len(l[3])

        seqid = line.strip() #entire line by default
        if seqid_col > 0:
            seqid = l[seqid_col-1]
        
        positions.append((chrom, snppos-1, snppos-1+len_refallele, seqid))

    f.close()
    return positions


def read_chromfile(filename):
    try:
        f = open(filename, 'r')
        seq = ''.join(map(lambda s: s.rstrip('\n').upper(), f.readlines()[1:]))
        f.close()

    except IOError as err:
        print("I/O error: ", err)
        sys.exit(0)

    return seq


def read_chromfile_gz(filename):
    try:
        f = gzip.open(filename, 'rt')
        seq = ''.join(map(lambda s: s.rstrip('\n').upper(), f.readlines()[1:]))
        f.close()

    except IOError as err:
        print("I/O error: ", err)
        sys.exit(0)

    return seq


def read_fasta_gz_file(filename):
    sequences = dict()
    chrname = None

    try:
        fp = gzip.open(filename, 'rt')
        for line in fp:
            if line[0] == '>':
                # compress the previous one
                if chrname:
                    sequences[chrname] = zlib.compress(pickle.dumps(sequences[chrname]))

                chrname = line[1:].rstrip('\n')
                sequences[chrname] = []
                print(chrname)
            else:
                sequences[chrname].append(line.rstrip('\n').upper())

        # process the last one
        sequences[chrname] = zlib.compress(pickle.dumps(sequences[chrname]))
        fp.close()
    except IOError as err:
        print("I/O error: ", err)
        sys.exit(0)
    
    return sequences

def main(argv = sys.argv):
    usage = "Usage: %prog [options] POSITION_FILE OUTFILE"
    desc  = "fetch genomic DNA sequences"

    parser = optparse.OptionParser(usage=usage, description=desc)                                                                              
    parser.add_option("-f", "--format", dest="format", default="bed", \
            help="set the input file format. supported types are (bed), (pos)ition, and (snp). (default=bed)")

    parser.add_option("-l", "--left", dest="lextend", type="int", default=0, \
            help="set the number of bases to extend in 5'(-) direction (default=0)")

    parser.add_option("-r", "--right", dest="rextend", type="int", default=0, \
            help="set the number of bases to extend in 3'(+) direction (default=0)")

    parser.add_option("-d", "--dir", dest="basedir", default='.', \
            help="set the base directory of reference sequences (default='.')")

    parser.add_option("-s", "--seqid", dest="seqid_column", type="int", default=0, \
            help="set a specific column (1-based) which will be used as sequence id (default=coordinate)")

    (options, args) = parser.parse_args()

    if len(args) == 0:
        parser.print_help()
        sys.exit(0)

    if len(args) != 2:
        parser.error("incorrect number of arguments")
        parser.print_help()
        sys.exit(0)

    position_fn = args[0]
    outfile = args[1]

    if options.format == "bed":
        positions = read_bedfile(position_fn, options.seqid_column)
    elif options.format == "pos":
        positions = read_posfile(position_fn, options.seqid_column)
    elif options.format == "snp":
        positions = read_snpfile(position_fn, options.seqid_column)
    else:
        sys.stderr.write('unknown file type: ' + options.format + '\n')
        sys.exit(0)

    chrnames = sorted(set(map(lambda p: p[0], positions)))
    
    fout = open(outfile, 'w')

    # read the genome sequence file if exists
    genomeseq = None
    genome_gz_fn = os.path.join(options.basedir, '.'.join(("chromosomes", 'fa', 'gz')))
    if os.path.isfile(genome_gz_fn):
        genomeseq = read_fasta_gz_file(genome_gz_fn)

    for chrom in chrnames:
        print(chrom)
        refseq_fn = os.path.join(options.basedir, '.'.join((chrom, 'fa')))
        refseq_gz_fn = os.path.join(options.basedir, '.'.join((chrom, 'fa', 'gz')))
        if genomeseq:
            refseq = ''.join(pickle.loads(zlib.decompress(genomeseq[chrom])))
        elif os.path.isfile(refseq_fn):
            refseq = read_chromfile(refseq_fn)
        elif os.path.isfile(refseq_gz_fn):
            refseq = read_chromfile_gz(refseq_gz_fn)
        else:
            sys.stderr.write('no sequence file in the given directory: ' + options.basedir + '\n')
            sys.exit(0)

        for pos in sorted(filter(lambda p: p[0] == chrom, positions), key=lambda pos: pos[1]):
            pos_start = max(0, pos[1]-options.lextend)
            pos_end = (pos[2]+options.rextend)
            if options.seqid_column > 0:
                fout.write(">%s\n" % (pos[3]))
            else:
                if options.format == "snp":
                    fout.write(">%s:%d-%d\t%s\n" % (chrom, pos_start+1, pos_end, pos[3])) #convert 0-based to 1-based and add SNP name
                else:
                    fout.write(">%s:%d-%d\n" % (chrom, pos_start+1, pos_end)) #convert 0-based to 1-based
            fout.write("%s\n" % refseq[pos_start:pos_end])
    fout.close()

if __name__ == '__main__': main()
