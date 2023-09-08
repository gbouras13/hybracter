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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int usage(char *);

int usage(char *prg) {
    printf("Usage: %s [options] reads\nOutputs two files part1_reads and part2_reads for bowtie consumption\nOptions: -ro for reads only (no quality info)\n", prg);
}

int main(int argc, char **argv){
    static char *prog;
    char *p;
    prog = argv[0];
    if((p = strrchr(prog, '/')) != NULL)
        prog = ++p;
    if (argc < 2) { usage(prog); return 1; }
    if (strcmp(argv[1], "-h") ==  0) { usage(prog); return 0; }
    if (strcmp(argv[1], "--helph")  == 0) { usage(prog); return 0; }
    printf("Input file to split: %s\n", argv[argc - 1]);
    
    int hasQualityInfo = 1;
    
    if (argc == 3) {
        if(strcmp(argv[1], "-ro")==0){
            hasQualityInfo = 0;
        }else{
            printf("Could not find option %s\n", argv[1]);
            return 0;
        }
    }
    
    // attempt to open the input file
    FILE *ins = fopen(argv[argc - 1], "r");
    if(ins == NULL){
        printf("Error! Could not open input file: %s\n", argv[argc - 1]);
        return 1;
    }
    
    
    
    // open up output file
    FILE *fo1, *fo2;
    char fileName1[100] = "part1_", fileName2[100] = "part2_";
    strcat(fileName1, argv[argc - 1]);
    strcat(fileName2, argv[argc - 1]);
    fo1 = fopen(fileName1, "w");
    fo2 = fopen(fileName2, "w");
    if(fo1 == NULL || fo2 == NULL){
        printf("Error! Could not open output files: %s or %s\n", fileName1, fileName2);
    }
    
    // read in and output the files
    int keepGoing = 1;
    char seqName[256];
    char seq[256];
    char qual[256];
    char temp[5];
    if(hasQualityInfo == 1){
        while(keepGoing > 0){
            // first read
            keepGoing = fscanf( ins, "%255s%255s%5s%255s", seqName, seq, temp, qual);
            if(keepGoing > 0){
                fprintf(fo1, "%s\n%s\n+\n%s\n", seqName, seq, qual);
                keepGoing = fscanf( ins, "%255s%255s%5s%255s", seqName, seq, temp, qual);
                fprintf(fo2, "%s\n%s\n+\n%s\n", seqName, seq, qual);
            }
        }
        fclose(fo1);
        fclose(fo2);
        fclose(ins);
    }else{
        while(keepGoing > 0){
            // first read
            keepGoing = fscanf( ins, "%255s%255s", seqName, seq);
            if(keepGoing > 0){
                fprintf(fo1, "%s\n%s\n", seqName, seq);
                keepGoing = fscanf( ins, "%255s%255s", seqName, seq);
                fprintf(fo2, "%s\n%s\n", seqName, seq);
            }
        }
        fclose(fo1);
        fclose(fo2);
        fclose(ins);
    }
    return 1;
}
