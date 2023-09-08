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

// cc -g -O2 synthReadGen2.c -o synthReadGen2 -lz -lm
// bowtie -t -I 0 -X 1000 --rf -a -l 10 -v 3 -S --sam-nohead --phred64-quals E_coli_first10k -1 part1_inwardBow100k.fna -2 part2_inwardBow100k.fna e_coli_bow.map

// $ ./synthReadGen2 -h
// Welcome to the Synthetic Read Generator of ALE!                                                                                                                                 
// (C) 2010 Scott Clark                                                                                                                                                            
// 
// Usage: ./synthReadGen [options] <inputFile> <outputFile>
// 
// Options: <i>nt <f>loat [default]
//   -h      : print out this help
//   -id <i> : set distribution used for insert length
//             [1 = normal], 2 = poisson
//   -ld <i> : set distribution used for read length
//             [1 = normal], 2 = poisson
//   -im <f> : inward insert length mean [200.0]
//   -om <f> : outward insert length mean [500.0]
//   -is <f> : inward insert length std dev [10.0]
//   -os <f> : outward insert length std dev [15.0]
//   -ip <f> : probability for an inward read [0.5]
//   -er <c> : illumina error char [^]
//   -nr <i> : number of reads to make [1000]
//   -rl <x> : read length mean [85.0]
//   -rs <x> : read length sigma [7.0]
//   -ps <x> : no error for first x bases in a read [0]
//   -b      : outputs two fastq files for bowtie mapping [off]

#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
#include <stdlib.h>
#include <time.h>
#include <math.h>
KSEQ_INIT(gzFile, gzread)

// get the compliment of any base pair
char getComplimentRes(const char res){
  if(res == 'A'){return 'T';}
  if(res == 'T'){return 'A';}
  if(res == 'G'){return 'C';}
  return 'G';
}

// translates an illumina char error rate to a numerical one
static double QtoP[63] = {0.0,0.205671765276,0.36904265552,0.498812766373,0.601892829447,0.683772233983,0.748811356849,0.800473768503,0.841510680754,0.874107458821,0.9,0.920567176528,0.936904265552,0.949881276637,0.960189282945,0.968377223398,0.974881135685,0.98004737685,0.984151068075,0.987410745882,0.99,0.992056717653,0.993690426555,0.994988127664,0.996018928294,0.99683772234,0.997488113568,0.998004737685,0.998415106808,0.998741074588,0.999,0.999205671765,0.999369042656,0.999498812766,0.999601892829,0.999683772234,0.999748811357,0.999800473769,0.999841510681,0.999874107459,0.9999,0.999920567177,0.999936904266,0.999949881277,0.999960189283,0.999968377223,0.999974881136,0.999980047377,0.999984151068,0.999987410746,0.99999,0.999992056718,0.999993690427,0.999994988128,0.999996018928,0.999996837722,0.999997488114,0.999998004738,0.999998415107,0.999998741075,0.999999,0.999999205672,0.999999369043};

int intFact(int temp){
  if(temp > 0){
    return intFact(temp - 1);
  }else{
    return 1;
  }
}

// returns the base, randomly substituted for another with error rate 1 - errorRate
char getBaseWithError(const char res, const double errorRate){
  if(errorRate > rand() / ( RAND_MAX + 1.0 )){
    //printf("True.\n");
      if(res == 'A' || res == 'T' || res == 'G' || res == 'C'){
          return res;
      }else{
          // ambiguous
          if(0.25 < rand() / ( RAND_MAX + 1.0 )){
              return 'A';
          }else{
              if(0.33333333 < rand() / ( RAND_MAX + 1.0 )){
                  return 'T';
              }else{
                  if(0.5 < rand() / ( RAND_MAX + 1.0 )){
                      return 'G';
                  }else{
                      return 'C';
                  }
              }
          }
      }
  }else{ // error!
    //printf("False.\n");
    if(res == 'A'){
      if(0.33333333 < rand() / ( RAND_MAX + 1.0 )){
    return 'T';
      }else{
    if(0.5 < rand() / ( RAND_MAX + 1.0 )){
      return 'G';
    }else{
      return 'C';
    }
      }
    }
    if(res == 'T'){
      if(0.33333333 < rand() / ( RAND_MAX + 1.0 )){
    return 'A';
      }else{
    if(0.5 < rand() / ( RAND_MAX + 1.0 )){
      return 'G';
    }else{
      return 'C';
    }
      }
    }
    if(res == 'G'){
      if(0.33333333 < rand() / ( RAND_MAX + 1.0 )){
    return 'T';
      }else{
    if(0.5 < rand() / ( RAND_MAX + 1.0 )){
      return 'A';
    }else{
      return 'C';
    }
      }
    }
    if(res == 'C'){
      if(0.33333333 < rand() / ( RAND_MAX + 1.0 )){
    return 'T';
      }else{
    if(0.5 < rand() / ( RAND_MAX + 1.0 )){
      return 'G';
    }else{
      return 'A';
    }
      }
    }
  }
  // ambiguous return a random base
  if(0.25 < rand() / ( RAND_MAX + 1.0 )){
      return 'A';
  }else{
      if(0.33333333 < rand() / ( RAND_MAX + 1.0 )){
          return 'T';
      }else{
          if(0.5 < rand() / ( RAND_MAX + 1.0 )){
              return 'G';
          }else{
              return 'C';
          }
      }
  }
}

// this is taken from pg 369 of numerical recipes
double NormalVal(float mu, float sigma){
  double u, v, x, y, q;
  int stop = 0;
  while(stop == 0){
    u = rand()/( RAND_MAX + 1.0 );
    v = 1.7156*(rand()/( RAND_MAX + 1.0 ) - 0.5);
    x = u - 0.449871;
    y = abs(v) + 0.386595;
    q = sqrt(x) + y*(0.19600*y-0.25472*x);
    if( q > 0.27597 && ( q > 0.27846 || sqrt(v) > -4.0*log(u)*sqrt(u))){
      stop = 1;
    }
  }
  return mu + sigma*v/u;
}

double TruncatedNormal(float mu, float sigma, float min, float max){
  double val = NormalVal(mu, sigma);
  if(val < min){
    return min;
  }else if(val > max){
    return max;
  }else{
    return val;
  }
}

// get the numerical value for the quality of a base call (illumina)
double getQualityP(const char quality){
  return QtoP[quality - 64];
}

static const char WELCOME_MSG[80] = "Welcome to the Synthetic Read Generator of ALE!\n(C) 2010 Scott Clark\n\n";
static const char SHORT_OPTIONS[80] = "    Options:\n    -h : print out help\n";
static const char LONG_OPTIONS[1024] = "Options: <i>nt <f>loat [default]\n  -h      : print out this help\n  -id <i> : set distribution used for insert length\n            [1 = normal], 2 = poisson\n  -ld <i> : set distribution used for read length\n            [1 = normal], 2 = poisson\n  -im <f> : inward insert length mean [200.0]\n  -om <f> : outward insert length mean [500.0]\n  -is <f> : inward insert length std dev [10.0]\n  -os <f> : outward insert length std dev [15.0]\n  -ip <f> : probability for an inward read [0.5]\n  -er <c> : illumina error char [^]\n  -nr <i> : number of reads to make [1000]\n  -rl <x> : read length mean [85.0]\n  -rs <x> : read length sigma [7.0]\n  -ps <x> : no error for first x bases in a read [0]\n  -b      : outputs two fastq files for bowtie mapping [off]\n";

int main(int argc, char **argv){
  if (argc == 1) {
    printf("%s", WELCOME_MSG);
    printf("Usage: %s [options] <inputFile> <outputFile>\n", argv[0]);
    printf("%s", SHORT_OPTIONS);
    return 1;
  }
  if (argc < 3){
    if(argv[1][0] == '-' && argv[1][1] == 'h'){
      printf("%s", WELCOME_MSG);
      printf("Usage: %s [options] <inputFile> <outputFile>\n\n", argv[0]);
      printf("%s", LONG_OPTIONS);
      return 0;
    }
    printf("%s", WELCOME_MSG);
    printf("Usage: %s [options] <inputFile> <outputFile>\n", argv[0]);
    printf("%s", SHORT_OPTIONS);
    return 1;
  }
  int i;
  int insDistToUse = 1; // distribution to use for insert lengths 1 = normal, 2 = poisson
  int lenDistToUse = 1; // distribution to use for read lengths 1 = normal, 2 = poisson
  int perfectStart = 0; // no error for first x bases in a read
  float inwardInsLenMean = 200.0; // inward insert length mean
  float inwardInsLenSigma = 10.0; // inward insert length sigma
  float inwardProb = 0.5; // probability of an inward insert
  float outwardInsLenMean = 500.0; // outward insert length mean
  float outwardInsLenSigma = 15.0; // inward insert length sigma
  float readLengthMean = 85.0; // read length mean
  float readLengthSigma = 7.0; // read length sigma
  int bowtieOutput = 0;
  char illuminaChar = '^';
  float errorRate = getQualityP(illuminaChar); // 1 - prob of substitution error
  int numReads = 1000;
  if(argv[1][0] == '-' && argv[1][1] == 'h'){
      printf("%s", WELCOME_MSG);
      printf("Usage: %s [options] <inputFile> <outputFile>\n\n", argv[0]);
      printf("%s", LONG_OPTIONS);
      return 0;
  }
  if (argc > 4){
    for(i = 1; i < argc - 2; i++){
      if(argv[i][0] == '-' && argv[i][1] == 'i' && argv[i][2] == 'd'){ // insert distance change
    insDistToUse = atoi(argv[i+1]);
      }else if(argv[i][0] == '-' && argv[i][1] == 'l' && argv[i][2] == 'd'){
    lenDistToUse = atoi(argv[i+1]);
      }else if(argv[i][0] == '-' && argv[i][1] == 'i' && argv[i][2] == 'm'){
    inwardInsLenMean = atof(argv[i+1]);
      }else if(argv[i][0] == '-' && argv[i][1] == 'i' && argv[i][2] == 'p'){
    inwardProb = atof(argv[i+1]);
      }else if(argv[i][0] == '-' && argv[i][1] == 'o' && argv[i][2] == 'm'){
    outwardInsLenMean = atof(argv[i+1]);
      }else if(argv[i][0] == '-' && argv[i][1] == 'e' && argv[i][2] == 'r'){
    illuminaChar = argv[i+1][0];
    errorRate = getQualityP(illuminaChar);
    //errorRate = atof(argv[i+1]);
      }else if(argv[i][0] == '-' && argv[i][1] == 'n' && argv[i][2] == 'r'){
    numReads = atof(argv[i+1]);
      }else if(argv[i][0] == '-' && argv[i][1] == 'i' && argv[i][2] == 's'){
    inwardInsLenSigma = atof(argv[i+1]);
      }else if(argv[i][0] == '-' && argv[i][1] == 'o' && argv[i][2] == 's'){
    outwardInsLenSigma = atof(argv[i+1]);
      }else if(argv[i][0] == '-' && argv[i][1] == 'r' && argv[i][2] == 'l'){
    readLengthMean = atof(argv[i+1]);
      }else if(argv[i][0] == '-' && argv[i][1] == 'r' && argv[i][2] == 's'){
    readLengthSigma = atof(argv[i+1]);
      }else if(argv[i][0] == '-' && argv[i][1] == 'p' && argv[i][2] == 's'){
    perfectStart = atoi(argv[i+1]);
      }else if(argv[i][0] == '-' && argv[i][1] == 'b'){
          bowtieOutput = 1;
          i--;
      }else{
    printf("Error! Could not recognize option: %s with value %s\n", argv[i], argv[i+1]);
    printf("Usage: %s [options] <inputFile> <outputFile>\n", argv[0]);
    printf("%s", SHORT_OPTIONS);
    return 0;
      }
      i++;
    }
  }
  printf("Starting with parameters: id = %i; ld = %i; im = %f; om = %f; ip = %f; er = %f (%c); is = %f; os = %f; rm = %f; rs = %f; ps = %i\n", insDistToUse, lenDistToUse, inwardInsLenMean, outwardInsLenMean, inwardProb,  errorRate, illuminaChar, inwardInsLenSigma, outwardInsLenSigma, readLengthMean, readLengthSigma, perfectStart);
  
  printf("Input file: %s\n", argv[argc - 2]);
  printf("Output file: %s\n", argv[argc - 1]);
  
  // open up input file (with the assembly to parse into reads)
  gzFile fp;
  kseq_t *seq;
  int l;
  fp = gzopen(argv[argc - 2], "r");
  seq = kseq_init(fp);
  printf("Read in file, now making sequences...\n");
  int seqLen;
  
  // open up output file
  FILE *fo1, *fo2, *foD;
  char fileName1[100] = "part1_", fileName2[100] = "part2_";
  char debugName[100] = "synthDebug_";
  if(bowtieOutput == 0){
      fo1 = fopen(argv[argc - 1], "w");
      fo2 = fo1;
  }else{
      strcat(fileName1, argv[argc - 1]);
      strcat(fileName2, argv[argc - 1]);
      strcat(debugName, argv[argc - 1]);
      fo1 = fopen(fileName1, "w");
      fo2 = fopen(fileName2, "w");
      foD = fopen(debugName, "w");
  }
  
  // seed the random number generator
  unsigned int iseed = (unsigned int)time(NULL);
  srand (iseed);
  
  int insertLength;
  int placement;
  int leftReadLength;
  int rightReadLength;
  int readOn = 0;

  while ((l = kseq_read(seq)) >= 0) {
    seqLen = (unsigned int)(seq->seq.l);
    printf("Name of assembly seq to parse: %s\n", seq->name.s);
    if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
    printf("length: %d\n", seqLen);

    int *placementDepth = (int*) calloc(seqLen, sizeof(int));
    for(i=0; i<seqLen; i++){
        placementDepth[i] = 0;
    }
    
    for(readOn = 0; readOn < numReads; readOn++){
      // decide if this is going to be an inward or outward facing read
      
      // INWARD CODE FIRST
      
      if(rand()/( RAND_MAX + 1.0 ) < inwardProb){ // we are making an inward facing read!
        if(insDistToUse == 1){ // use normal distribution
          insertLength = (int)TruncatedNormal(inwardInsLenMean, inwardInsLenSigma, 10, seqLen);
        }
        // now we want to place it uniformly
        placement = (int)(rand()/( RAND_MAX + 1.0 )*(seqLen - insertLength));
        if(lenDistToUse == 1){ // use normal distribution
          leftReadLength = (int)TruncatedNormal(readLengthMean, readLengthSigma, 10, insertLength);
          rightReadLength = (int)TruncatedNormal(readLengthMean, readLengthSigma, 10, insertLength);
        }

        for(i = 0; i < leftReadLength; i++){
            placementDepth[placement + i]++;
        }
        for(i = 0; i < rightReadLength; i++){
            placementDepth[placement + insertLength - 1 - i]++;
        }
    
        if(rand()/( RAND_MAX + 1.0 ) < 0.5){ // forward first
          fprintf(fo1, "@synthR%iFi1,p:%i,i:%i\n", readOn, placement, insertLength);
          for(i = 0; i < leftReadLength; i++){
            if(i < perfectStart){
              fprintf(fo1, "%c",getBaseWithError(seq->seq.s[placement + i], 1.0));
            }else{
              fprintf(fo1, "%c",getBaseWithError(seq->seq.s[placement + i], errorRate));
            }
          }
          fprintf(fo1, "\n+\n");
          for(i = 0; i < leftReadLength; i++){fprintf(fo1, "%c",illuminaChar);}
          fprintf(fo1, "\n");
          fprintf(fo2, "@synthR%iFi1,p:%i,i:%i\n", readOn, placement, insertLength);
          for(i = 0; i < rightReadLength; i++){
            if(i < perfectStart){
              fprintf(fo2, "%c",getBaseWithError(getComplimentRes(seq->seq.s[placement + insertLength - 1 - i]), 1.0));
            }else{
              fprintf(fo2, "%c",getBaseWithError(getComplimentRes(seq->seq.s[placement + insertLength - 1 - i]), errorRate));
            }
          }
          fprintf(fo2, "\n+\n");
          for(i = 0; i < rightReadLength; i++){fprintf(fo2, "%c",illuminaChar);}
          fprintf(fo2, "\n");
        }else{ // backward first
          fprintf(fo1, "@synthR%iFi1,p:%i,i:%i\n", readOn, placement, insertLength);
          for(i = 0; i < rightReadLength; i++){
            if(i < perfectStart){
              fprintf(fo1, "%c",getBaseWithError(getComplimentRes(seq->seq.s[placement + insertLength - 1 - i]), 1.0));
            }else{
              fprintf(fo1, "%c",getBaseWithError(getComplimentRes(seq->seq.s[placement + insertLength - 1 - i]), errorRate));
            }
          }
          fprintf(fo1, "\n+\n");
          for(i = 0; i < rightReadLength; i++){fprintf(fo1, "%c",illuminaChar);}
          fprintf(fo1, "\n");
          fprintf(fo2, "@synthR%iFi1,p:%i,i:%i\n", readOn, placement, insertLength);
          for(i = 0; i < leftReadLength; i++){
            if(i < perfectStart){
              fprintf(fo2, "%c",getBaseWithError(seq->seq.s[placement + i], 1.0));
            }else{
              fprintf(fo2, "%c",getBaseWithError(seq->seq.s[placement + i], errorRate));
            }
          }
          fprintf(fo2, "\n+\n");
          for(i = 0; i < leftReadLength; i++){fprintf(fo2, "%c",illuminaChar);}
          fprintf(fo2, "\n");
        }
      }else{ // make an outward facing read

        // OUTWARD FACING READ

        if(insDistToUse == 1){ // use normal distribution
          insertLength = (int)TruncatedNormal(outwardInsLenMean, outwardInsLenSigma, 10, seqLen);
        }
        // now we want to place it uniformly
        placement = (int)(rand()/( RAND_MAX + 1.0 )*(seqLen - insertLength));
        if(lenDistToUse == 1){ // use normal distribution
          leftReadLength = (int)TruncatedNormal(readLengthMean, readLengthSigma, 10, insertLength);
          rightReadLength = (int)TruncatedNormal(readLengthMean, readLengthSigma, 10, insertLength);
        }

        for(i = 0; i < leftReadLength; i++){
            placementDepth[placement + leftReadLength - i]++;
        }
        for(i = 0; i < rightReadLength; i++){
            placementDepth[placement + insertLength - 1 - rightReadLength + i]++;
        }

        if(rand()/( RAND_MAX + 1.0 ) < 0.5){ // the left read first
          fprintf(fo1, "@synthR%iL,p:%i,i:%i\n", readOn, placement, insertLength);
          for(i = 0; i < leftReadLength; i++){
            if(i < perfectStart){
              fprintf(fo1, "%c",getBaseWithError(getComplimentRes(seq->seq.s[placement + leftReadLength - i]), 1.0));
            }else{
              fprintf(fo1, "%c",getBaseWithError(getComplimentRes(seq->seq.s[placement + leftReadLength - i]), errorRate));
            }
          }
          fprintf(fo1, "\n+\n");
          for(i = 0; i < leftReadLength; i++){fprintf(fo1, "%c",illuminaChar);}
          fprintf(fo1, "\n");
          fprintf(fo2, "@synthR%iL,p:%i,i:%i\n", readOn, placement, insertLength);
          for(i = 0; i < rightReadLength; i++){
            if(i < perfectStart){
              fprintf(fo2, "%c",getBaseWithError(seq->seq.s[placement + insertLength - 1 - rightReadLength + i], 1.0));
            }else{
              fprintf(fo2, "%c",getBaseWithError(seq->seq.s[placement + insertLength - 1 - rightReadLength + i], errorRate));
            }
          }
          fprintf(fo2, "\n+\n");
          for(i = 0; i < rightReadLength; i++){fprintf(fo2, "%c",illuminaChar);}
          fprintf(fo2, "\n");
        }else{ // the right read first
          fprintf(fo1, "@synthR%iL,p:%i,i:%i\n", readOn, placement, insertLength);
          for(i = 0; i < rightReadLength; i++){
            if(i < perfectStart){
              fprintf(fo1, "%c",getBaseWithError(seq->seq.s[placement + insertLength - 1 - rightReadLength + i], 1.0));
            }else{
              fprintf(fo1, "%c",getBaseWithError(seq->seq.s[placement + insertLength - 1 - rightReadLength + i], errorRate));
            }
          }
          fprintf(fo1, "\n+\n");
          for(i = 0; i < rightReadLength; i++){fprintf(fo1, "%c",illuminaChar);}
          fprintf(fo1, "\n");
          fprintf(fo2, "@synthR%il,p:%i,i:%i\n", readOn, placement, insertLength);
          for(i = 0; i < leftReadLength; i++){
            if(i < perfectStart){
              fprintf(fo2, "%c",getBaseWithError(getComplimentRes(seq->seq.s[placement + leftReadLength - i]), 1.0));
            }else{
              fprintf(fo2, "%c",getBaseWithError(getComplimentRes(seq->seq.s[placement + leftReadLength - i]), errorRate));
            }
          }
          fprintf(fo2, "\n+\n");
          for(i = 0; i < leftReadLength; i++){fprintf(fo2, "%c",illuminaChar);}
          fprintf(fo2, "\n");
        }
      }
    }

    for(i=0; i<seqLen; i++){
        fprintf(foD, "%d\n",placementDepth[i]);
    }
    free(placementDepth);
  }
  fclose(fo1);
  if (bowtieOutput != 0) {
	  fclose(fo2);
	  fclose(foD);
  }
  return 0;
}
