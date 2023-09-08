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

// geneTree.c

#include "geneTree.h"

int OutputIndicies(treeBranch_t *pRoot, const char sequence[], const int klen, int index[][2]){
  int i;
  treeBranch_t *currentBranch = pRoot;
  treeBranch_t *nextBranch;
  for(i = 0; i < klen; i++){
    if(sequence[i] == 'A'){
      if(currentBranch->subBranches[0] == NULL){
	return -1;
      }else{
	nextBranch = currentBranch->subBranches[0];
      }
    }
    if(sequence[i] == 'T'){
      if(currentBranch->subBranches[1] == NULL){
	return -1;
      }else{
	nextBranch = currentBranch->subBranches[1];
      }
    }
    if(sequence[i] == 'C'){
      if(currentBranch->subBranches[2] == NULL){
	return -1;
      }else{
	nextBranch = currentBranch->subBranches[2];
      }
    }
    if(sequence[i] == 'G'){
      if(currentBranch->subBranches[3] == NULL){
	return -1;
      }else{
	nextBranch = currentBranch->subBranches[3];
      }
    }
    currentBranch = nextBranch;
  }
  //printf("Number indicies: %i.\n", currentBranch->leaf->current);
  for(i = 0; i < currentBranch->leaf->current; i++){
    index[i][0] = currentBranch->leaf->indicies[i][0];
    index[i][1] = currentBranch->leaf->indicies[i][1];
  }
  return currentBranch->leaf->current;
}

treeBranch_t MakeTree(const assembly_t theAssem, const int klen, const int NUM_ASSEMBLY_PARTS){
  
  int i, j;
  treeBranch_t theRoot, *pRoot = &theRoot;
  theRoot.subBranches[0] = NULL;
  theRoot.subBranches[1] = NULL;
  theRoot.subBranches[2] = NULL;
  theRoot.subBranches[3] = NULL;
  theRoot.leaf = NULL;
  for(i = 0; i < NUM_ASSEMBLY_PARTS; i++){
    for(j = 0; j < theAssem.assemblyParts[i].seqLen - klen; j++){
      //printf("Adding sequence %i, %i to tree.\n", i, j);
      //printf("Attempting to add to tree: %i.\n", j);
      AddSeqToTree(theAssem.assemblyParts[i].sequence, j, klen, pRoot, i);
    }
  }
  return theRoot;
}

int AddSeqToTree(const char sequence[], const int offset, const int klen, treeBranch_t *pRoot, const int assemPart){
  //printf("Offset: %i, klen: %i, assemPart: %i\n", offset, klen, assemPart);
  int i;
  char currentRes;
  treeBranch_t *currentBranch = pRoot;
  treeBranch_t *nextBranch;
  leafNode_t *theLeaf, *oldLeaf;
  for(i = 0; i < klen; i++){
    //currentRes = getCharFromSeqByLoc(sequence, offset + i);
    currentRes = sequence[offset + i];
    //printf("Adding %c (res: %i).\n", currentRes, offset + i);
    if(currentRes == 'A'){
      if(currentBranch->subBranches[0] == NULL){
	//malloc another treebranch
	nextBranch = malloc(sizeof (treeBranch_t));
	nextBranch->subBranches[0] = NULL;
	nextBranch->subBranches[1] = NULL;
	nextBranch->subBranches[2] = NULL;
	nextBranch->subBranches[3] = NULL;
	nextBranch->leaf = NULL;
	currentBranch->subBranches[0] = nextBranch;
      }
      currentBranch = currentBranch->subBranches[0];
    }
    if(currentRes == 'T'){
      if(currentBranch->subBranches[1] == NULL){
	//malloc another treebranch
	nextBranch = malloc(sizeof (treeBranch_t));
	nextBranch->subBranches[0] = NULL;
	nextBranch->subBranches[1] = NULL;
	nextBranch->subBranches[2] = NULL;
	nextBranch->subBranches[3] = NULL;
	nextBranch->leaf = NULL;
	currentBranch->subBranches[1] = nextBranch;
      }
      currentBranch = currentBranch->subBranches[1];
    }
    if(currentRes == 'C'){
      if(currentBranch->subBranches[2] == NULL){
	//malloc another treebranch
	nextBranch = malloc(sizeof (treeBranch_t));
	nextBranch->subBranches[0] = NULL;
	nextBranch->subBranches[1] = NULL;
	nextBranch->subBranches[2] = NULL;
	nextBranch->subBranches[3] = NULL;
	nextBranch->leaf = NULL;
	currentBranch->subBranches[2] = nextBranch;
      }
      currentBranch = currentBranch->subBranches[2];
    }
    if(currentRes == 'G'){
      if(currentBranch->subBranches[3] == NULL){
	//malloc another treebranch
	nextBranch = malloc(sizeof (treeBranch_t));
	nextBranch->subBranches[0] = NULL;
	nextBranch->subBranches[1] = NULL;
	nextBranch->subBranches[2] = NULL;
	nextBranch->subBranches[3] = NULL;
	nextBranch->leaf = NULL;
	currentBranch->subBranches[3] = nextBranch;
      }
      currentBranch = currentBranch->subBranches[3];
    }
    if(currentRes == 'N'){
      return 0;
    }
  }

  //currentRes = getCharFromSeqByLoc(sequence, offset + klen - 1);
  //printf("Making a leaf.\n");
    if(currentBranch->leaf == NULL){
      //printf("malloc the leaf indicies.\n");
      theLeaf = malloc(2*sizeof(int) + 2*sizeof(int)*START_LENGTH);
      //printf("Initialize the leaf.\n");
      currentBranch->leaf = theLeaf;
      currentBranch->leaf->length = START_LENGTH;
      currentBranch->leaf->current = 0;
    }else{
      if(currentBranch->leaf->current > currentBranch->leaf->length - 1){
	// malloc a bigger leaf
	// printf("Making a bigger leaf.");
	oldLeaf = currentBranch->leaf;
	theLeaf = malloc(2*sizeof(int) + 2*sizeof(int)*(currentBranch->leaf->length + START_LENGTH));
	theLeaf->current = oldLeaf->current;
	theLeaf->length = currentBranch->leaf->length + START_LENGTH;
	for(i = 0; i < theLeaf->current; i++){
	  theLeaf->indicies[i][0] = oldLeaf->indicies[i][0];
	  theLeaf->indicies[i][1] = oldLeaf->indicies[i][1];
	}
	currentBranch->leaf = theLeaf;
	free(oldLeaf);
      }
    }
    //printf("Set indicies.\n");
    currentBranch->leaf->indicies[currentBranch->leaf->current][0] = assemPart;
    currentBranch->leaf->indicies[currentBranch->leaf->current][1] = offset;
    currentBranch->leaf->current += 1;
    return 1;
}

