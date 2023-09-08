#!/usr/bin/env python

"""
/*
 * Copyright (C) 2010,2011,2012 Scott Clark. All rights reserved.
 *
 * Developed by:
 * Scott Clark
 * Cornell University Center for Applied Mathematics
 * http://cam.cornell.edu
 * AND
 * Rob Egan
 * Department of Energy Joint Genome Institute
 * http://jgi.doe.gov
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a 
 * copy of this software and associated documentation files (the "Software"), 
 * to deal with the Software without restriction, including without limitation 
 * the rights to use, copy, modify, merge, publish, distribute, sublicense, 
 * and/or sell copies of the Software, and to permit persons to whom the 
 * Software is furnished to do so, subject to the following conditions:
 *
 *   1. Redistributions of source code must retain the above copyright notice, 
 *      this list of conditions and the following disclaimers.
 *   2. Redistributions in binary form must reproduce the above copyright 
 *      notice, this list of conditions and the following disclaimers in the 
 *      documentation and/or other materials provided with the distribution.
 *   3. Neither the names of Cornell University, The Joint Genome Institute, 
 *      nor the names of its contributors may be used to endorse or promote 
 *      products derived from this Software without specific prior written 
 *      permission.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
 * CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER 
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING 
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS WITH THE SOFTWARE.
 */

// For more information on the license please see 
// The University of Illinois/NCSA Open Source License
// http://www.opensource.org/licenses/UoI-NCSA.php
"""

def convertToWiggle(inFile):
	kwig = open(inFile + "-kmer.wig", 'w')
	pwig = open(inFile + "-place.wig", 'w')
	iwig = open(inFile + "-insert.wig", 'w')
	dwig = open(inFile + "-depth.wig", 'w')
	wig = open(inFile + ".wig", 'w')

	pwig.write("track name=ALE-place color=0,0,255 group=ALE priority=1\n")
	iwig.write("track name=ALE-insert color=255,0,255 group=ALE priority=1\n")
	dwig.write("track name=ALE-depth color=255,0,0 group=ALE priority=2\n")
	kwig.write("track name=ALE-kmer color=0,255,0 group=ALE priority=3\n")
	wig.write("track name=depth color=0,0,0 group=ALE priotity=4\n")

	for line in file(inFile):
		line = line.rstrip()
		if line[0] == "#":
			sp = line.split()
			if sp[1] == "Reference:":
				pwig.write("fixedStep chrom=" + sp[2] + " start=1 step=1\n")
				iwig.write("fixedStep chrom=" + sp[2] + " start=1 step=1\n")
				dwig.write("fixedStep chrom=" + sp[2] + " start=1 step=1\n")
				kwig.write("fixedStep chrom=" + sp[2] + " start=1 step=1\n")
				wig.write("fixedStep chrom=" + sp[2] + " start=1 step=1\n")
			continue
		contig,position,depth,depthLike,placeLike,insertLike,kmerLike = line.split()
		kwig.write(str(kmerLike) + "\n")	
		pwig.write(str(placeLike) + "\n")
		iwig.write(str(insertLike) + "\n")	
		dwig.write(str(depthLike) + "\n")	
		wig.write(str(depth) + "\n")

def main(argv):
        """ """
        convertToWiggle(argv[1])

def usage(exit_value):
    print >> sys.stderr, "ale2wiggle.py ale_report.ale"
    sys.exit(exit_value)

if __name__ == "__main__":
    import sys
    import os
    if len(sys.argv) == 1: usage(1)
    if sys.argv[1] == '-h' or sys.argv[1] == '--help': usage(0)
    if not os.path.isfile(sys.argv[1]): 
       print >> sys.stderr, 'No suchfile or directory:', sys.argv[1]
       sys.exit(0)
    sys.exit(main(sys.argv))
    

