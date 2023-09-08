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

// ALE.h

#ifndef _ALE_H_
#define _ALE_H_

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "ALEhelpers.h"
#include "ALElike.h"
#include "geneTree.h"

static const char WELCOME_MSG[] = "Welcome to the Assembly Likelihood Estimator!\n(C) 2010 Scott Clark\n\n";
static const char USAGE[] = "Usage: ALE [-options] alignments.[s|b]am assembly.fasta[.gz] ALEoutput.txt\n";
#define DEFAULT_REALIGN_OPTIONS "1,3,11,4,8"
static const char LONG_OPTIONS[] =
		"Options: <i>nt <f>loat <s>tring [default]\n"
		"-h or --help    : print out this help\n"
		"--kmer <f>      : Kmer depth for kmer stats [4]\n"
		"--qOff <i>      : Quality ascii offset (illumina) [33] or 64 (or 0)\n"
		"--pl <s>        : placementOutputBAM\n"
		"--pm <s>        : library parameter file (auto outputs .param)\n"
		"--nout          : only output meta information (no per base) [off]\n"
		"--minLL         : the minimum log Likelihood (-120)\n"
		"--metagenome    : Evaluate each contig independently for depth & kmer metrics\n"
		"--realign[=matchScore,misMatchPenalty,gapOpenPenalty,gapExtPenalty,minimumSoftClip (default: " DEFAULT_REALIGN_OPTIONS ") ]\n"
		"                   Realign reads with Striped-Smith-Waterman honoring ambiguous reference bases\n"
		"                   and stacking homo-polymer indels\n"
		"                   for PacBio, try --realign=1,5,2,1,20 (similar to BWA-SW recommendations)\n"
		"--SNPreport <s> : Creates a new text file reporting all SNP phasing \n"
		"                   observed by a read against ambiguous bases in the reference\n"
		"--minQual <i>   : Minimum quality score to use in Z-normalization (default 3).\n"
		"                   Illumina quality scores can be unreliable below this threshold\n"
		"\n\n";

#endif
