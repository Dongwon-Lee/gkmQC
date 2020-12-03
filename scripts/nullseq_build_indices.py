import os
import sys
import tarfile
import zipfile
import gzip
import optparse

from bitarray import bitarray

def clear_indexes(sid, buildname):
    na = '.'.join([buildname, sid, "na", "out"])
    gc = '.'.join([buildname, sid, "gc", "out"])
    rpt = '.'.join([buildname, sid, "rpt", "out"])

    #truncate files
    for fn in (na, gc, rpt):
        f = open(fn, 'wb')
        f.close()


def append_indexes(seq, sid, buildname):
    na = '.'.join([buildname, sid, "na", "out"])
    gc = '.'.join([buildname, sid, "gc", "out"])
    rpt = '.'.join([buildname, sid, "rpt", "out"])

    f = open(na, 'ab')
    bitarray(map(lambda c: c in 'nN', seq)).tofile(f)
    f.close()

    f = open(gc, 'ab')
    bitarray(map(lambda c: c in 'cgCG', seq)).tofile(f)
    f.close()

    f = open(rpt, 'ab')
    bitarray(map(lambda c: c in 'acgt', seq)).tofile(f)
    f.close()


def build_indexes(fn, buildname):
    save_interval = 8*32*1024

    try:
        f = open(fn, 'r')
        
        seq = [] 
        sid = ''
        nlines = 0
        for line in f:
            if line[0] == '>':
                if sid: 
                    append_indexes("".join(seq), sid, buildname)
                    seq = []

                sid = line[1:].rstrip('\n').split()[0]
                clear_indexes(sid, buildname)
            else:
                nlines += 1
                seq.append(line.rstrip('\n'))

                if nlines % save_interval == 0:
                    append_indexes("".join(seq), sid, buildname)
                    seq = []

        #the last remaining sequence 
        append_indexes("".join(seq), sid, buildname)

    except IOError as err:
        print("I/O error: ", err)
        sys.exit(0)

def build_indexes_fp(fp, buildname):
    save_interval = 8*32*1024

    try:
        seq = [] 
        sid = ''
        nlines = 0
        for line in fp:
            if line[0] == '>':
                if sid: 
                    append_indexes("".join(seq), sid, buildname)
                    seq = []

                sid = line[1:].rstrip('\n').split()[0]
                clear_indexes(sid, buildname)
            else:
                nlines += 1
                seq.append(line.rstrip('\n'))

                if nlines % save_interval == 0:
                    append_indexes("".join(seq), sid, buildname)
                    seq = []

        #the last remaining sequence 
        append_indexes("".join(seq), sid, buildname)

    except IOError as err:
        print("I/O error: ", err)
        sys.exit(0)


def main(argv=sys.argv):
    usage = "usage: %prog [options] <Chromosome File(TARBALL gzip (tar.gz) or zip)> <Genome Build Name>"
    desc  = "generate bit index files for generating null sequences"

    parser = optparse.OptionParser(usage=usage, description=desc)

    parser.add_option("-q", dest="quiet", default=False, action="store_true", \
            help="supress messages (default=false)")

    (options, args) = parser.parse_args()

    if len(args) == 0:
        parser.print_help()
        sys.exit(0)

    if len(args) != 2:
        parser.error("incorrect number of arguments")
        parser.print_help()
        sys.exit(0)

    chrom_file = args[0]
    genome = args[1]

    if zipfile.is_zipfile(chrom_file):
        if options.quiet == False:
            sys.stderr.write("detected file type is zip.\n")

        zipfileobj = zipfile.ZipFile(chrom_file)

        for fn in zipfileobj.namelist():
            if options.quiet == False:
                sys.stderr.write(' '.join(["processing", fn, "\n"]))

            zipfileobj.extract(fn)
            build_indexes(fn, genome)
            os.remove(fn)

        zipfileobj.close()

    elif tarfile.is_tarfile(chrom_file):
        if options.quiet == False:
            sys.stderr.write("detected file type is tar.\n")

        tarfileobj = tarfile.open(chrom_file)

        for tarinfo in tarfileobj:
            if not tarinfo.isdir():
                fn = tarinfo.name
                if options.quiet == False:
                    sys.stderr.write(' '.join(["processing", fn, "\n"]))

                tarfileobj.extract(fn)
                build_indexes(fn, genome)
                os.remove(fn)

        tarfileobj.close()
    elif chrom_file[-2:] == "gz":
        fp = gzip.open(chrom_file, 'rt')
        build_indexes_fp(fp, genome)
        fp.close()
    else:
        sys.stderr.write(' '.join(["unknown input file:", chrom_file, "\n"]))

if __name__ == "__main__": main()
