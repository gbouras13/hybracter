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

// ALE.c

#include "ALE.h"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <getopt.h>
#include "ALEhelpers.h"
#include "ALElike.h"

static struct option long_options[] = {
		{"help", 0, 0, 0},
		{"kmer", 1, 0, 0},
		{"qOff", 1, 0, 0},
		{"pl", 1, 0, 0},
		{"pm", 1, 0, 0},
		{"nout", 0, 0, 0},
		{"minLL", 0, 0, 0},
		{"metagenome", 0, 0, 0},
		{"realign", 2, 0, 0},
		{"SNPreport", 1, 0, 0},
        {0, 0, 0, 0}
};

void usage() {
    printf("%s", WELCOME_MSG);
    printf("%s", USAGE);
    printf("\n%s", LONG_OPTIONS);
}
int main(int argc, char **argv){
	int c = 0, digit_optind = 0;
    int kmerLen = 4;
    char placementOut[256];
    samfile_t *placementBam = NULL;
    *placementOut = '\0';
    libraryParametersT *libParams = NULL;
    double outlierFraction = 0.02;
    int qOff = -1;
    int printAllAleOutput = 1;

    int doRealign = 0;
    // use BWA bwasw default matrix
    char defaultRealignOptions[] = DEFAULT_REALIGN_OPTIONS;
    int matchScore, mismatchPenalty, gapOpenPenalty, gapExtendPenalty, minSoftClip;
    FILE *snpPhaseFile = NULL;


    while(1) {
        int this_option_optind = optind ? optind : 1;
        int option_index = 0;
        c = getopt_long(argc, argv, "h",
                 long_options, &option_index);
        if (c == -1)
        	break;
        switch(c) {
        case 0:
        	// long options
        	if (strcmp(long_options[option_index].name, "kmer") == 0) {
                kmerLen = atoi(optarg);
                if(kmerLen > 20){
                    printf("--kmer option of %i not in range [2,20], set to default [4].\n", kmerLen);
                    kmerLen = 4;
                }
        	} else if (strcmp(long_options[option_index].name, "qOff") == 0) {
                qOff = atoi(optarg);
                if(qOff != 64 && qOff != 33 && qOff != 0){
                    printf("--qOff option of %i not in set [33,64], will be set to 33.\n", qOff);
                    qOff = 33;  // SAM/BAM specification is for ascii - 33.
                }
        	} else if (strcmp(long_options[option_index].name, "pl") == 0) {
                strcpy(placementOut, optarg);
        	} else if (strcmp(long_options[option_index].name, "pm") == 0) {
                libParams = malloc(sizeof(libraryParametersT));
                importLibraryParameters(libParams, optarg);
        	} else if (strcmp(long_options[option_index].name, "nout") == 0) {
                printAllAleOutput = 0;
                printf("Turned off per-base output\n");
        	} else if (strcmp(long_options[option_index].name, "minLL") == 0) {
            	setMinLogLike(atof(optarg));
            	printf("Set minLogLike to: %0.1f\n", getMinLogLike());
        	} else if (strcmp(long_options[option_index].name, "metagenome") == 0) {
            	setMetagenome();
            	printf("This dataset will be evaluated as a Metagenome\n");
        	} else if (strcmp(long_options[option_index].name, "realign") == 0) {
            	doRealign = 1;
            	char *swopts = defaultRealignOptions;
            	if (optarg) {
            		swopts = optarg;
            	}
            	// parse SW options...
            	if (5 != sscanf(swopts, "%d,%d,%d,%d,%d", &matchScore, &mismatchPenalty, &gapOpenPenalty, &gapExtendPenalty, &minSoftClip)) {
            		usage();
            		fprintf(stderr, "Could not parse %s as --realign parameter.  Expecting 5 integers.\n", optarg);
            		exit(1);
            	}
            	printf("This dataset will be realigned with Smith-Waterman letting ambiguous bases match %s\n", swopts);
        	} else if (strcmp(long_options[option_index].name, "SNPreport") == 0) {
        		snpPhaseFile = fopen(optarg, "w");
        		if (snpPhaseFile == NULL) {
        			fprintf(stderr, "Could not open '%s' for writing the SNPreport!\n", optarg);
        			exit(1);
        		}
        		printSNPhaseHeader(snpPhaseFile);
        		printf("Reporting SNP phasing into %s\n", optarg);
        	} else if (strcmp(long_options[option_index].name, "minQual") == 0) {
        		fprintf(stderr, "Using minimum quality score of %s\n", optarg);
        		setMinimumQuality( atoi(optarg) );
        		exit(0);
         	} else if (strcmp(long_options[option_index].name, "help") == 0) {
        		usage();
        		exit(0);
         	} else {
         		fprintf(stderr, "unknown option: %s\n", long_options[option_index].name);
         		usage();
         		exit(1);
         	}
        	break;
        case 'h':
        default:
        	usage();
        	exit(0);
        }
    }

    if (argc - optind < 3) {
        // the input was shorter than anticipated, print long options if -h specified
    	usage();
        return 0;
    }
    
    // input and output files
    char *bamFile = argv[optind], *assemblyFile = argv[optind+1], *aleFile = argv[optind+2];
    printf("BAM file: %s\n", bamFile);
    printf("Assembly fasta file: %s\n", assemblyFile);
    printf("ALE Output file: %s\n", aleFile);

    // attempt to open the bam input file
    samfile_t *ins = openSamOrBam(bamFile);
    
    printf("Reading in assembly...\n");
    assemblyT *theAssembly = loadAssembly(assemblyFile);
    if (!validateAssemblyIsSameAsAlignment(ins->header, theAssembly)) {
        printf("Error! Assembly fasta %s does not match alignment file %s!\n", assemblyFile, bamFile);
        exit(1);
    }
    
    if (*placementOut != '\0') {
        printf("Placement file: %s\n", placementOut);
        char mode[] = "wbu";
        if (placementOut[0] != '-' )
            mode[2] = '\0'; // compress bam if not part of a pipe
        placementBam = samopen(placementOut, mode, ins->header);
        if(placementBam == 0){
            printf("Error! Could not write to the placement BAM file: %s\n", placementOut);
            exit(1);
        }
    }

    printf("Reading in the map and computing statistics...\n");

    if(libParams != NULL){
      saveLibraryParameters(libParams, aleFile); // the ALE output file name
    }

    // calculate the insert mean/std if not given
    if(libParams == NULL){
        printf("Insert length and std not given, will be calculated from input map.\n");

        libParams = computeLibraryParameters(ins, outlierFraction, qOff);
        saveLibraryParameters(libParams, aleFile); // the ALE output file name

        // close and re-open bam file
        samclose(ins);
        ins = openSamOrBam(bamFile);
    }

    // set assembly avgReadSize to that of library (rounded to nearest int)
    theAssembly->avgReadSize = (double)libParams->avgReadSize;

    // place reads and compute statistics on the assembly
    printf("Computing read placements and depths\n");
    if (doRealign)
    	initRefNumForRealign(theAssembly, matchScore, mismatchPenalty, gapOpenPenalty, gapExtendPenalty, minSoftClip);
    computeReadPlacements(ins, theAssembly, libParams, placementBam, snpPhaseFile);
    if (doRealign)
    	destroyRefNumForRealign(theAssembly);
    if (snpPhaseFile != NULL) {
    	if (ftell(snpPhaseFile) > 0) {
    		fprintf(snpPhaseFile, "\n");
    	}
    	fclose(snpPhaseFile);
    	snpPhaseFile = NULL;
    }
    
    // compute statistics on assembly
    printf("Computing k-mer statistics...\n");
    computeKmerStats(theAssembly, kmerLen);
    printf("Done computing k-mer statistics.\n");
    
    printf("Computing depth statistics...\n");
    computeDepthStats(theAssembly, libParams);
    printf("Done computing depth statistics.\n");

    printf("Computing statistics on expected missing...\n");
    applyExpectedMissingLength(theAssembly);
    printf("Done computing statistics on expected missing.\n");
    
    FILE *out = fopen(aleFile, "w");
    if(out == NULL){
        printf("Error! Could not open output file: %s\n", aleFile);
    }
    writeToOutput(theAssembly, printAllAleOutput, out);
    fclose(out);
    
    printf("Output is in file: %s\n", aleFile);

    if (placementBam != NULL) {
        printf("Closing placement file\n");
        samclose(placementBam);
    }

    printf("Closing input file\n");
    samclose(ins);
    
    if (libParams != NULL)
        free(libParams);

    freeAssembly(theAssembly);
    return 0;
}
