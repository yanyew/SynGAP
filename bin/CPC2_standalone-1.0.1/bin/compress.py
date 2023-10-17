'''
This module deals with compressed file (.gz or .bz2)
'''

import gzip
import bz2
import sys
def gz_file(fq_file,mode,level=6):
	try:
		if fq_file.endswith("gz"):
			fq_fp = gzip.open(fq_file,mode+"b",level)
		else:
			sys.stderr.write("[INFO] read file '%s'\n"%fq_file)
			fq_fp = open(fq_file,mode)
	except:
		sys.stderr.write("Error: Fail to IO file: %s\n"%(fq_file))
		sys.exit(1)
	return fq_fp


def bz2file(f):
	fz = None
	if f.endswith("bz2"):
		fz = bz2.BZ2File(f)
	else:
		sys.stderr.write("Error: Fail to IO file: %s\n"%(f))
		sys.exit(1)
	return fz

