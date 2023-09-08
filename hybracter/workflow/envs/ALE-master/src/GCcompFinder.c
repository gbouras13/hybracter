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
#include <math.h>

int main(int argc, char **argv){
  
  if (argc < 4) {
    printf("Usage: %s assemblyFile.fasta runningWindow output\n(C) 2011 Scott Clark\nFinds the running window GC content in a fasta file\nJust sequence! No gene name/info!\n", argv[0]);
    return 0;
  }
  
  // attempt to open the first input file
  FILE *ins = fopen(argv[argc - 3], "r");
  if(ins == NULL){
      printf("Error! Could not open the first input file: %s\n", argv[argc - 3]);
      return 1;
  }
  
  // attempt to open the output file
  FILE *fout = fopen(argv[argc - 1], "w");
  if(fout == NULL){
      printf("Error! Could not open the output file: %s\n", argv[argc - 1]);
      return 1;
  }
  
  int runningWindow = atoi(argv[argc - 2]);
  int i;
  char base;
  char *backTrace = malloc(sizeof(char)*runningWindow);
  int gcTot = 0;
  int keepGoing;
  
  for(i = 0; i < runningWindow; i++){
    keepGoing = fscanf(ins, "%c", &base);
//     if(base == '\0' || base == '\n'){
//       i--;
//     }
    backTrace[i] = base;
    if(base == 'G' || base == 'g' || base == 'C' || base == 'c'){
      gcTot++;
    }
  }
  
  printf("%f\n", (float)gcTot/((float)runningWindow));
  
  while(keepGoing > 0){
    keepGoing = fscanf(ins, "%c", &base);
    if(base == '\0' || base == '\n'){
      i--;
    }
    if(base == 'G' || base == 'g' || base == 'C' || base == 'c'){
      gcTot++;
    }
  }
}
