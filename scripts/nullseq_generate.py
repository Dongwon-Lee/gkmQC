#!/usr/bin/env python
"""
	nullseq_generate.py: sample random genomic sequences by matching the 
	given profile 

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


def bitarray_fromfile(filename):
	"""
	"""
	fh = open(filename, 'rb')
	bits = bitarray()
	bits.fromfile(fh)

	return bits, fh


def make_profile(positions, buildname, basedir):
	"""
	"""
	chrnames = sorted(set(map(lambda p: p[0], positions)))

	profiles = []
	for chrom in chrnames:
		idxf_gc = os.path.join(basedir, '.'.join([buildname, chrom, 'gc', 'out']))
		idxf_rpt = os.path.join(basedir, '.'.join([buildname, chrom, 'rpt', 'out']))

		#if os.path.exists(idxf_gc) == False or os.path.exists(idxf_rpt) == False:
		#	continue
		bits_gc, gcf = bitarray_fromfile(idxf_gc)
		bits_rpt, rptf = bitarray_fromfile(idxf_rpt)

		for pos in positions:
			if pos[0] != chrom:
				continue

			seqid = pos[0] + ':' + str(pos[1]+1) + '-' + str(pos[2])
			slen = pos[2]-pos[1]
			gc = bits_gc[pos[1]:pos[2]].count(True)
			rpt = bits_rpt[pos[1]:pos[2]].count(True)

			profiles.append((seqid, slen, gc, rpt))

		gcf.close()
		rptf.close()

	return profiles


def sample_sequences(positions, buildname, basedir, options):
	"""
	"""
	rpt_err = options.rpt_err
	gc_err = options.gc_err
	max_trys = options.max_trys
	norpt = options.norpt
	nogc = options.nogc

	chrnames = sorted(set(map(lambda p: p[0], positions)))
	profiles = make_profile(positions, buildname, basedir)

	excluded = []
	if options.excludefile:
		excluded = read_bed_file(options.excludefile)

	#truncate (clear) file
	f = open(options.output, 'w')
	f.close()
	
	for chrom in chrnames:
		if options.quiet == False:
			sys.stderr.write("sampling from " + chrom + "\n")

		idxf_na = os.path.join(basedir, '.'.join([buildname, chrom, 'na', 'out']))
		idxf_gc = os.path.join(basedir, '.'.join([buildname, chrom, 'gc', 'out']))
		idxf_rpt = os.path.join(basedir, '.'.join([buildname, chrom, 'rpt', 'out']))

		#this bit array is used to mark positions that are excluded from sampling
		#this will be updated as we sample more sequences in order to prevent sampled sequences from overlapping
		bits_na, naf = bitarray_fromfile(idxf_na)
		bits_gc, gcf = bitarray_fromfile(idxf_gc)
		bits_rpt, rptf = bitarray_fromfile(idxf_rpt)

		#mark excluded regions
		for pos in excluded:
			if pos[0] != chrom:
				continue
			bits_na[pos[1]:pos[2]] = True

		npos = 0
		#mark positive regions
		for pos in positions:
			if pos[0] != chrom:
				continue
			bits_na[pos[1]:pos[2]] = True
			npos+=1

		if options.count == 0:
			count = options.fold*npos
		else:
			count = options.count

		sampled_positions = []
		while len(sampled_positions) < count:
			sampled_prof = random.choice(profiles) # positive region
			sampled_len = sampled_prof[1] # To be controlled 
			sampled_gc = sampled_prof[2] # ""
			sampled_rpt = sampled_prof[3] 

			rpt_err_allowed = int(rpt_err*sampled_len)
			gc_err_allowed = int(gc_err*sampled_len)
			trys = 0
			while trys < max_trys:
				trys += 1

				pos = random.randint(1, bits_na.length() - sampled_len) # random choose
				pos_e = pos+sampled_len

				if bits_na[pos:pos_e].any():
					continue

				if not norpt:
					pos_rpt = bits_rpt[pos:pos_e].count(True)
					if abs(sampled_rpt - pos_rpt) > rpt_err_allowed:
						continue

				if not nogc:
					pos_gc = bits_gc[pos:pos_e].count(True)
					if abs(sampled_gc - pos_gc) > gc_err_allowed:
						continue

				#accept the sampled position
				#mark the sampled regions
				bits_na[pos:pos_e] = True

				sampled_positions.append((chrom, pos, pos_e))

				#print trys, chrom, pos, pos_e, sampled_len, pos_rpt, sampled_rpt, pos_gc, sampled_gc
				break
			else:
				if options.quiet == False:
					sys.stderr.write(' '.join(["fail to sample from", \
							"len=", str(sampled_len), \
							"rpt=", str(sampled_rpt), \
							"gc=", str(sampled_gc)]) + '\n')

		naf.close()
		gcf.close()
		rptf.close()

		f = open(options.output, 'a')
		for spos in sorted(sampled_positions, key=lambda s: s[1]):
			f.write('\t'.join([spos[0], str(spos[1]), str(spos[2])]) + '\n')
		f.close()


def main(argv=sys.argv):
	usage = "usage: %prog [options] <Input File (bed or pos)> <Genome Build Name> <Base Directory of Index Files>"

	desc  = "generate null sequences"
	parser = optparse.OptionParser(usage=usage, description=desc)

	parser.add_option("-i", dest="input_type", \
		default = "bed", help="set the type of input file. bed and pos are available (default=bed)")

	parser.add_option("-x", dest="fold", type="int", \
		default = 1, help="number of sequence to sample, FOLD times of given dataset (default=1)")

	parser.add_option("-c", dest="count", type="int", \
		default=0, help="number of sequences to sample, override -x option (default=NA, obsolete)")

	parser.add_option("-r", dest="rseed", type="int", \
		default=1, help="random number seed (default=1)")

	parser.add_option("-g", dest="gc_err", type="float", \
		default=0.02, help="GC errors allowed (default=0.02)")

	parser.add_option("-t", dest="rpt_err", type="float", \
		default=0.02, help="RPT errors allowed (default=0.02)")

	parser.add_option("-e", dest="excludefile", \
		default="", help="filename that contains regions to be excluded (default=NA)")

	parser.add_option("-G", dest="nogc", action="store_true", \
		default=False, help="do not match gc-contents")

	parser.add_option("-R", dest="norpt", action="store_true", \
		default=False, help="do not match repeats")

	parser.add_option("-m", dest="max_trys", type="int", \
		default=10000, help="number of maximum trys to sample of one sequence (default=10000)")

	parser.add_option("-o", dest="output", default="nullseq_output.bed", \
			help="set the name of output file (default=nullseq_output.bed)")

	parser.add_option("-q", dest="quiet", default=False, action="store_true", \
  			help="supress messages (default=false)")

	(options, args) = parser.parse_args()

	if len(args) == 0:
		parser.print_help()
		sys.exit(0)

	if len(args) != 3:
		parser.error("incorrect number of arguments")
		parser.print_help()
		sys.exit(0)


	posfile = args[0]
	buildname = args[1]
	basedir = args[2]

	random.seed(options.rseed)

	if options.quiet == False:
		sys.stderr.write('Options:\n')
		sys.stderr.write('  input type: ' + options.input_type + '\n')
		sys.stderr.write('  N fold: ' + str(options.fold) + '\n')
		sys.stderr.write('  GC match: ' + str(not options.nogc) + '\n')
		sys.stderr.write('  repeat match: ' + str(not options.norpt) + '\n')
		sys.stderr.write('  GC error allowed: ' + str(options.gc_err) + '\n')
		sys.stderr.write('  repeat error allowed: ' + str(options.rpt_err) + '\n')
		sys.stderr.write('  random seed: ' + str(options.rseed) + '\n')
		sys.stderr.write('  max trys: ' + str(options.max_trys) + '\n')
		sys.stderr.write('  excluded regions: ' + options.excludefile+ '\n')
		sys.stderr.write('  output: ' + options.output + '\n')
		sys.stderr.write('\n')

		sys.stderr.write('Input args:\n')
		sys.stderr.write('  input file: ' + posfile+ '\n')
		sys.stderr.write('  genome build name: ' + buildname+ '\n')
		sys.stderr.write('  basedir: ' + basedir + '\n')
		sys.stderr.write('\n')

	if options.input_type == 'bed':
		positions = read_bed_file(posfile)
	elif options.input_type == 'pos':
		positions = read_pos_file(posfile)
	else:
		parser.error(options.input_type + " is not supported")
		parser.print_help()
		sys.exit(0)
		
	sample_sequences(positions, buildname, basedir, options)
		
if __name__ == "__main__": main()
