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

// ALEhelpers.c

#include "ALEhelpers.h"
#include "ALElike.h"

// initialize ALEhelpers static variables
static double _minLogLike = -120.0;
static int _metagenome = 0;
static double _leastIdentity = 0.95; // 95% identity is floor, unless there are worse reads that map
static long _multiple_primary_alignments = 0;
static int _minimumQuality = 3;

void setMinLogLike(double min) {
	_minLogLike = min;
}
double getMinLogLike() {
	return _minLogLike;
}
double validateLogLikelihood(const double logLikelihood) {
	if (logLikelihood < getMinLogLike() || isnan(logLikelihood) || isinf(logLikelihood))
		return getMinLogLike();
	else
		return logLikelihood;
}
void setMetagenome() {
	_metagenome = 1;
}
int isMetagenome() {
	return _metagenome != 0;
}
void setLeastIdentity(double identity) {
	if (identity < _leastIdentity)
		_leastIdentity = identity;
}
double getLeastIdentity() {
	return _leastIdentity;
}

void setMultiplePrimaryAlignments() {
	_multiple_primary_alignments++;
}

long getMultiplePrimaryAlignments() {
	return _multiple_primary_alignments;
}

void setMinimumQuality(int min) {
	_minimumQuality = min;
}
int getMinimumQuality() {
	return _minimumQuality;
}

double getQtoP(char qualChar, int qOff) {
    int idx = qualChar - qOff;
    if (idx < 0 || idx >= 63 ) {
        printf("WARNING: getQtoP called out of range: %c %d %d\n", qualChar, qOff, idx);
        if (idx < 0) idx = 0;
        else idx = 62;
    }
    return QtoP[idx];
}
double getQtoLogP(char qualChar, int qOff) {
    int idx = qualChar - qOff;
    if (idx < 0 || idx >= 63 ) {
        printf("WARNING: getQtoLogP called out of range: %c %d %d\n", qualChar, qOff, idx);
        if (idx < 0) idx = 0;
        else idx=62;
    }
    return QtoLogP[idx];
}
double getQtoLogPMiss(char qualChar, int qOff) {
    int idx = qualChar - qOff;
    if (idx < 0 || idx >= 63 ) {
        printf("WARNING: getQtoLogPMiss called out of range: %c %d %d\n", qualChar, qOff, idx);
        if (idx < 0) idx = 0;
        else idx = 62;
    }
    return QtoLogPMiss[idx]; // TODO switch to (1-Q)*Q
}

void IncreaseAssemblyPartsByOne(assembly_t *theAssembly, int numParts){
  assemblyPart_t *tempPartPointer = malloc(numParts* sizeof(assemblyPart_t));
  int i;
  for(i = 0; i < numParts - 1; i++){
    tempPartPointer[i] = theAssembly->assemblyParts[i];
  }
  //free(theAssembly->assemblyParts);
  theAssembly->assemblyParts = tempPartPointer;
  //return theAssembly;
}

double poissonInt(int k, double lambda){
  int i;
  double answer = 1.0;
  for(i = 1; i < k+1; i++){
    answer = answer*lambda*exp(-lambda/(float)k)/(float)i;
  }
  return answer;
}

//uses Stirlings approximation to high precision
double lnfact(double input){
  return (input - 0.5)*log(input) - input + lnfactconst - 1.0/(12.0*input) - 1.0/(360.0*input*input*input) - 1.0/(1260.0*input*input*input*input*input);
}

// convert a 4-mer into a byte
unsigned char seqToChar(const char pos1, const char pos2, const char pos3, const char pos4){
    unsigned char seqChar = 0;
    seqChar += (pos1 == 'A' || pos1 == 'G');
    seqChar += (pos1 == 'A' || pos1 == 'C')*2;
    seqChar += (pos2 == 'A' || pos2 == 'G')*4;
    seqChar += (pos2 == 'A' || pos2 == 'C')*8;
    seqChar += (pos3 == 'A' || pos3 == 'G')*16;
    seqChar += (pos3 == 'A' || pos3 == 'C')*32;
    seqChar += (pos4 == 'A' || pos4 == 'G')*64;
    seqChar += (pos4 == 'A' || pos4 == 'C')*128;
    return seqChar;
}

// saves the quality information of a sequence
void makeQual(const char qual[], const unsigned int seqLen, char quality[], int qualityHistogram[]){
  int i;
  for(i = 0; i < seqLen; i++){
    quality[i] = qual[i];
    qualityHistogram[qual[i] - 64]++;
  }
}

// saves the quality information of a sequence pretending it is as high as it can be
void makeQualPerfect(const unsigned int seqLen, char quality[]){
  int i;
  for(i = 0; i < seqLen; i++){
    quality[i] = 127;
  }
}

// saves the quality information of a sequence
void makeAssemblySeq(const char seq[], const unsigned int seqLen, unsigned char sequence[]){
  int i;
  for(i = 0; i < seqLen; i++){
    sequence[i] = toupper(seq[i]);
  }
}

// converts a sequence to the condensed char equivalent
void makeSeq(const char seq[], const unsigned int seqLen, unsigned char sequence[]){
    int i;
    // make all multiples of 4 possible
    for(i = 0; i < seqLen/4; i++){
        sequence[i] = seqToChar(seq[i*4], seq[i*4 + 1], seq[i*4 + 2], seq[i*4 + 3]);
    }
    if(seqLen%4 == 3){
        sequence[seqLen/4] = seqToChar(seq[i*4], seq[i*4 + 1], seq[i*4 + 2], 'T');
    }else if(seqLen%4 == 2){
        sequence[seqLen/4] = seqToChar(seq[i*4], seq[i*4 + 1], 'T', 'T');
    }else if(seqLen%4 == 1){
        sequence[seqLen/4] = seqToChar(seq[i*4], 'T', 'T', 'T');
    }
}

// gives you the nucleotide at position loc from the sequence <seq>, not efficient for many close
// calls, use charToSeq to get 4 at a time.
char getCharFromSeqByLoc(const unsigned char seq[], const unsigned int loc){
    char subSeq[4];
    charToSeqFour(seq[loc/4], subSeq); // grab the 4-mer that we want
    return subSeq[loc%4]; // return the actual residue we want
}



void charToSeqFour(unsigned char num, char seq[]){
  seq[0] = theFourConverter[num][0];
  seq[1] = theFourConverter[num][1];
  seq[2] = theFourConverter[num][2];
  seq[3] = theFourConverter[num][3];
}

// unwraps a byte <num> into it's <len> residues and stores them in <seq>
// useful for output or direct manipulation of the sequence
void charToSeq(unsigned char num, char seq[], const int len){
    int i;
    for(i = 0; i < len; i++){
        if(num%2){ // A or G
            num = num >> 1;
            if(num%2){ // A or C match?
                seq[i] = 'A';
            }else{
                seq[i] = 'G';
            }
        }else{
            num = num >> 1;
            if(num%2){ // A or C match?
                seq[i] = 'C';
            }else{
                seq[i] = 'T';
            }
        }
        num = num >> 1;
    }
    for(i = len; i < 4; i++){ // fill in the rest with spaces
        seq[i] = ' ';
    }
}

// prints out the sequence in it's numeric (byte-wrapped) form
void PrintSequenceNumeric(const unsigned char sequence[], const unsigned int seqLen){
    int j;
    for(j = 0; j < seqLen/4; j++){
        //printf("%i ", sequence[j]);
    }
    if(seqLen%4 != 0){ //print the end
        //printf("%i", sequence[seqLen/4]);
    }
    //printf("\n");
}

void PrintSequenceB(const unsigned char sequence[], const unsigned int seqLen){
    int j;
    for(j = seqLen-1; j > 0; j--){
        //printf("%c", getCharFromSeqByLoc(sequence, j));
    }
    //printf("\n");
}

void PrintSequence(const unsigned char sequence[], const unsigned int seqLen){
    char seqer[4];
    int j;
    for(j = 0; j < seqLen/4; j++){
        charToSeq(sequence[j], seqer, 4);
        //printf("%.4s", seqer);
    }
    if(seqLen%4 != 0){
        charToSeq(sequence[seqLen/4], seqer, seqLen%4);
        //printf("%.4s", seqer);
    }
    //printf("\n");
}

void PrintQuality(const char quality[], const unsigned int seqLen){
  int i;
  for(i = 0; i < seqLen; i++){
    //printf("%c", quality[i]);
  }
  //printf("\n");
}

void PrintAssembly(const char sequence[], const unsigned int seqLen){
  int i;
  for(i = 0; i < seqLen; i++){
    //printf("%c", sequence[i]);
  }
  //printf("\n");
}

// get the numerical value for the quality of a base call
double getQualityP(const char quality[], const unsigned int i){
  return QtoP[quality[i] - 64];
}

/* This is a secret function, its magics are UNKNOWN */
int intMax(int a, int b){
    if(a > b){
        return a;
    }
    return b;
}

/* This is a secret function, its magics are UNKNOWN */
int intMin(int a, int b){
    if(a < b){
        return a;
    }
    return b;
}

char getComplimentRes(const char res){
  if(res == 'A'){return 'T';}
  if(res == 'T'){return 'A';}
  if(res == 'G'){return 'C';}
  return 'G';
}

void PrintPlacements(pairedRead_t theRead){
  int i;
  for(i = 0; i < theRead.numPlacements; i++){
    //printf("L: %f, o1: %i, o2: %i, i: %i an: %i\n", theRead.placements[i].likelihood, theRead.placements[i].offset1, theRead.placements[i].offset2, theRead.placements[i].placeInfo, theRead.placements[i].assemPart);
  }
}

void initAlignment(alignSet_t *dst) {
    dst->likelihood = 0.0;
    dst->likelihoodInsert = 0.0;
    dst->start1 = -1;
    dst->start2 = -1;
    dst->end1 = -1;
    dst->end2 = -1;
    dst->contigId1 = -1;
    dst->contigId2 = -1;
    dst->name = NULL;
    dst->orientation = NO_READS;
    dst->nextAlignment = NULL;
}

void destroyAlignment(alignSet_t *dst) {
    assert(dst != NULL);
    if (dst->name != NULL)
        free(dst->name);
}

void copyAlignment(alignSet_t *dst, const alignSet_t *src) {
    assert(dst != NULL);
    assert(src != NULL);
    destroyAlignment(dst);

    dst->likelihood = src->likelihood;
    dst->likelihoodInsert = src->likelihoodInsert;
    dst->start1 = src->start1;
    dst->start2 = src->start2;
    dst->end1 = src->end1;
    dst->end2 = src->end2;
    dst->contigId1 = src->contigId1;
    dst->contigId2 = src->contigId2;
    if (src->name != NULL)
        dst->name = strdup(src->name);
    else
        dst->name = NULL;
    dst->orientation = src->orientation;
    dst->nextAlignment = src->nextAlignment;
    dst->bamOfAlignment1 = src->bamOfAlignment1;
    dst->bamOfAlignment2 = src->bamOfAlignment2;
}

void printAlignment(const alignSet_t *src) {
    fprintf(stderr, "l: %f li: %f, s1: %d, s2: %d, e1: %d, e2: %d, c1: %d, c2: %d, %s, %s, %p, b1: %p, b2: %p\n", src->likelihood, src->likelihoodInsert, src->start1, src->start2, src->end1, src->end2, src->contigId1, src->contigId2, src->name, MATE_ORIENTATION_LABELS[src->orientation], src->nextAlignment, src->bamOfAlignment1, src->bamOfAlignment2);
}

void swap(void **x, void **y) {
    void *t = *x;
    *x = *y;
    *y = t;
}

int isGC(char seq) {
    return (seq == 'G' || seq == 'C' || seq == 'g' || seq == 'c');
}

int getGCtotal(char seq1[], int seq1len){
    int GCtot = 0, i;
    for(i = 0; i < seq1len; i++){
        if(isGC(seq1[i])){
            GCtot++;
        }
    }
    return GCtot;
}

int findNumAssemPieces(kseq_t *ins){
    int l, count = 0;
    while ((l = kseq_read(ins)) >= 0) {
        count++;
    }
    return count;
}

void readAssembly(kseq_t *ins, assemblyT *theAssembly){
    int contigLen, i, l, j = 0;

    // read contigs into a linked list.  Consolidate to an array later
    contig_ll *tmp, *head = malloc(sizeof(contig_ll));
    head->contig = NULL;
    head->next = NULL;
    tmp = head;
    theAssembly->totalAssemLen = 0;

    while ((l = kseq_read(ins)) >= 0) {
        contigLen = (int)(ins->seq.l);
        theAssembly->totalAssemLen += (long)(ins->seq.l);
        //printf("Found contig %d, contigLen = %i, name=%s\n", j, contigLen, ins->name.s);

        contig_t *contig = tmp->contig = (contig_t*) malloc(sizeof(contig_t));
        contig->seqLen = contigLen;
        contig->isCircular = 0; // assume nothing is circular until there reasonable evidence
        contig->seq = malloc(contigLen*sizeof(char));
        contig->seqNum = NULL; // only use if realigning sequences
        contig->ambiguousBaseCount = 0;
        contig->ambiguousBasePositions = NULL;
        contig->depth = malloc(contigLen*sizeof(float));
        contig->matchLogLikelihood = malloc(contigLen*sizeof(float));
        contig->insertLogLikelihood = malloc(contigLen*sizeof(float));
        contig->depthLogLikelihood = malloc(contigLen*sizeof(float));
        contig->kmerLogLikelihood = malloc(contigLen*sizeof(float));
        contig->name = strdup(ins->name.s);
        for(i = 0; i < contigLen; i++){
            contig->seq[i] = toupper(ins->seq.s[i]);
            contig->depth[i] = 0.0;
            contig->matchLogLikelihood[i] = 0.0;
            contig->insertLogLikelihood[i] = 0.0;
            contig->depthLogLikelihood[i] = 0.0;
            contig->kmerLogLikelihood[i] = 0.0;
        }
        j++;
        tmp->next = malloc(sizeof(contig_ll));
        tmp = tmp->next;
        tmp->contig = NULL;
        tmp->next = NULL;

    }
    int numberAssemblyPieces = j;
    //printf("Found %d contigs\n", numberAssemblyPieces);

    theAssembly->contigs = (contig_t**)malloc((numberAssemblyPieces)*sizeof(contig_t*));
    theAssembly->numContigs = numberAssemblyPieces;
    
    theAssembly->totalScore = 0.0;
    theAssembly->kmerAvgSum = 0.0;
    theAssembly->kmerAvgNorm = 0.0;
    theAssembly->placeAvgSum = 0.0;
    theAssembly->placeAvgNorm = 0.0;
    theAssembly->insertAvgSum = 0.0;
    theAssembly->insertAvgNorm = 0.0;
    theAssembly->depthScoreAvgSum = 0.0;
    theAssembly->depthScoreAvgNorm = 0.0;
    theAssembly->depthAvgSum = 0.0;
    theAssembly->depthAvgNorm = 0.0;
    theAssembly->overlapAvgSum = 0.0;
    theAssembly->overlapAvgNorm = 0.0;
    theAssembly->totalReads = 0;
    theAssembly->totalUnmappedReads = 0;
    theAssembly->totalMappedReads = 0;
    theAssembly->totalPlacedReads = 0;
    theAssembly->avgReadSize = 0.0;

    // consolidate linked list into array, free linked list
    tmp = head;
    i = 0;
    while (tmp != NULL) {
        if (tmp->contig != NULL)
            theAssembly->contigs[i++] = tmp->contig;
        head = tmp;
        tmp = head->next;
        free(head);
    }

}

int validateAssemblyIsSameAsAlignment(bam_header_t *header, assemblyT *theAssembly) {
    int i;
    //printf("Validating assembly and alignment files consisting of %d contigs\n", header->n_targets);
    if (header->n_targets != theAssembly->numContigs) {
        //printf("Different number of contigs in assembly (%d) and alignmentfile (%d)\n", theAssembly->numContigs, header->n_targets);
        return 0;
    }
    for (i = 0; i < header->n_targets; i++) {
        if (header->target_len[i] != theAssembly->contigs[i]->seqLen) {
        	//printf("Different contig length in contig %d: %d vs %d\n", i, header->target_len[i], theAssembly->contigs[i]->seqLen);
        	return 0;
        }
        if (strcmp(header->target_name[i], theAssembly->contigs[i]->name) != 0) {
        	//printf("Warning assembly and alignment files disagree on the name of contig %d: %s vs %s\n", i, header->target_name[i], theAssembly->contigs[i]->name);
        }
    }
    return 1;
}

// returns a new array of GCcont for this contig.  User must free after using.
// below is my attempt at a hanning window convolution, I coded it from scratch so watch for bugs!
unsigned char *calculateContigGCcont(contig_t *contig, int windowSize) {
	int j, baseGC;
	int *GCpast = malloc(sizeof(int) * windowSize);
	unsigned char *GCcont = (unsigned char*) malloc(sizeof(unsigned char) * contig->seqLen);
	if (contig->seqLen < 2 * windowSize) {
		// contig is too small to estimate per-base windowing
		baseGC = getGCtotal(contig->seq, contig->seqLen);
		baseGC = floor(100.0*(double)baseGC / (double) contig->seqLen);
		for(j=0; j < contig->seqLen ; j++) {
			GCcont[j] = baseGC;
		}
	} else {
		baseGC = getGCtotal(contig->seq, windowSize);
		GCpast[0] = baseGC;
		for(j = 0; j < windowSize - 1; j++){
			GCcont[j] = floor(100.0*(double)baseGC/(double)((j+1)*windowSize));
			if (GCcont[j] > 100)
				printf("GC correction out of range %d %s %d len %d '%c'\n", baseGC, contig->name, j, contig->seqLen, contig->seq[j]);
			GCpast[(j+1)%windowSize] = GCpast[j%windowSize];
			if(isGC(contig->seq[j])){
				GCpast[(j+1)%windowSize]--;
			}
			if(isGC(contig->seq[j+windowSize])){
				GCpast[(j+1)%windowSize]++;
			}
			baseGC += GCpast[(j+1)%windowSize];
		}
		for(j = windowSize - 1; j < contig->seqLen - windowSize; j++){
			GCcont[j] = floor(100.0*(double)baseGC/(double)(windowSize*windowSize));
			if (GCcont[j] > 100)
				printf("GC correction out of range %d %s %d len %d '%c'\n", baseGC, contig->name, j, contig->seqLen, contig->seq[j]);
			baseGC -= GCpast[(j+1)%windowSize];
			GCpast[(j+1)%windowSize] = GCpast[j%windowSize];
			if(isGC(contig->seq[j])){
				GCpast[(j+1)%windowSize]--;
			}
			if(isGC(contig->seq[j+windowSize])) {
				GCpast[(j+1)%windowSize]++;
			}
			baseGC += GCpast[(j+1)%windowSize];
		}
		for(j = contig->seqLen - windowSize; j < contig->seqLen; j++){
			GCcont[j] = floor(100.0*(double)baseGC/(double)((contig->seqLen - j)*windowSize));
			if (GCcont[j] > 100)
				printf("GC correction out of range %d %s %d len %d '%c'\n", baseGC, contig->name, j, contig->seqLen, contig->seq[j]);
			baseGC -= GCpast[(j+1)%windowSize];
		}
	}
	free(GCpast);
	return GCcont;
}

double getContigAvgDepth(contig_t *contig) {
	double depth = 0.0;
	int j;
	for(j=0; j < contig->seqLen; j++)
		depth += contig->depth[j];
	return depth / contig->seqLen;
}

int getSeqMapLenBAM(bam1_t *read) {
    assert(read != NULL);
    return bam_cigar2qlen(&read->core, bam1_cigar(read));
}

int getFragmentMapLenBAM(bam1_t *read1) {
    assert(read1 != NULL);

    int left = read1->core.pos < read1->core.mpos ? read1->core.pos : read1->core.mpos;
    int readLength = getSeqMapLenBAM(read1);
    int right1 = read1->core.pos + readLength;
    int right2 = read1->core.mpos + readLength;
    int right = right1 > right2 ? right1 : right2;
    assert(right >= left);
    return right - left;
}

enum MATE_ORIENTATION getPairedMateOrientation(bam1_t *read1) {
    char paired = (read1->core.flag & BAM_FPAIRED) == BAM_FPAIRED;
    char read1unmap = ((read1->core.flag & BAM_FUNMAP)  == BAM_FUNMAP)  | (read1->core.tid < 0);
    char read2unmap = ((read1->core.flag & BAM_FMUNMAP) == BAM_FMUNMAP) | ( paired && read1->core.mtid < 0);

    if (read1unmap | read2unmap) {
        // read or mate is not mapped
        if (paired) {
        	// paired
        	if (read1unmap && read2unmap) {
        		// neither read is mapped
        		return UNMAPPED_PAIR;
        	} else {
        		// only one read in the pair is mapped
        		if ((read1->core.flag & (BAM_FREAD1|BAM_FREAD2)) != 0) {
        			return read1unmap ? UNMAPPED_MATE : SINGLE_UNMAPPED_MATE;
        		} else {
        			// not simply a pair of two reads, treat like single (think PacBio...)
        			return read1unmap ? UNMAPPED_SINGLE : SINGLE_READ;
        		}
        	}
        } else {
        	// not paired, so single.
        	return read1unmap ? UNMAPPED_SINGLE : SINGLE_READ;
        }
    }

    // this read (and potential mate) is mapped
    assert((read1->core.flag & BAM_FUNMAP) != BAM_FUNMAP);

    if (!paired) {
    	//printf("%s is single\n", bam1_qname(read1));
        return SINGLE_READ;
    }
    if (((read1->core.flag & (BAM_FREAD1 | BAM_FREAD2)) == 0) || (read1->core.isize == 0 && read1->core.tid >= 0 && read1->core.tid == read1->core.mtid)) {
        // required for PacBio reads that claim they are pairs but are not in a strict sense
    	//printf("%s is single (forced)\n", bam1_qname(read1));
    	return SINGLE_READ;
    }
    assert((read1->core.flag & BAM_FPAIRED) == BAM_FPAIRED);

    char isProper = ((read1->core.flag & BAM_FPROPER_PAIR) == BAM_FPROPER_PAIR);
    if (read1->core.tid == read1->core.mtid) {
        // reads map to same contig
        int read1Dir = (read1->core.flag & BAM_FREVERSE) == BAM_FREVERSE ? 1 : 0;
        int read2Dir = (read1->core.flag & BAM_FMREVERSE) == BAM_FMREVERSE ? 1 : 0;
        if (read1Dir == read2Dir)
        	return (isProper ? VALID_FF : NOT_PROPER_FF);
        else {
        	// TODO rethink this if read sizes are different could use read1->core.isize instead
        	int readLength = getSeqMapLenBAM(read1);
        	if (read1Dir == 0) {
        		if (read1->core.pos <= read1->core.mpos + readLength)
        			return (isProper ? VALID_FR : NOT_PROPER_FR);
                else
        		    return (isProper ? VALID_RF : NOT_PROPER_RF);
        	} else {
        		if (read1->core.mpos <= read1->core.pos + readLength)
        			return (isProper ? VALID_FR : NOT_PROPER_FR);
                else
        		    return (isProper ? VALID_RF : NOT_PROPER_RF);
        	}
        }
    } else {
        // reads map to different contig
        return CHIMER;
    }

}

#ifndef BAM_FSUPPLEMENTARY
#define BAM_FSUPPLEMENTARY 2048
#endif

enum MATE_ORIENTATION readNextBAM(samfile_t *ins, bam1_t *read1) {
    assert(ins != NULL);
    assert(read1 != NULL);

    while (1) {
      int bytesRead = samread(ins, read1);
      if (bytesRead <= 0) {
        return NO_READS;
      } else if ( ( read1->core.flag & ( BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP | BAM_FSUPPLEMENTARY)) == 0) { // only consider the read if it is primary, passedqc and not duplicated
        enum MATE_ORIENTATION o = getPairedMateOrientation(read1);
        //printf("readNextBam: %s %d %d %s\n", bam1_qname(read1), read1->core.flag, o, MATE_ORIENTATION_LABELS[o]);
        return o;
      }
    }
}


// prints out all of the alignments in the linked list
void printAlignments(alignSet_t *head){
    // print out the head
    //printf("Alignment 1 for read %s: %f at %i-%i and %i-%i.\n", head->name, head->likelihood, head->start1, head->end1, head->start2, head->end2);
    alignSet_t *current = head;
    int i = 1;
    while(current->nextAlignment != NULL){
        current = current->nextAlignment;
        i++;
        //printf("Alignment %i for read %s: %f at %i-%i and %i-%i.\n", i, current->name, current->likelihood, current->start1, current->end1, current->start2, current->end2);
    }
}

void writeToOutput(assemblyT *theAssembly, int fullOut, FILE *out){
    int i, j;
    printf("Set minLogLike to: %0.1f\n", getMinLogLike());
    //printf("Writing statistics to output file.\n");
    fprintf(out, "# ALE_score: %lf\n", theAssembly->totalScore);
    fprintf(out, "# numContigs: %d\n", theAssembly->numContigs);
    fprintf(out, "# totalAssemLen: %ld\n", theAssembly->totalAssemLen);
    fprintf(out, "# placeAvg: %lf\n", theAssembly->placeAvgSum/theAssembly->placeAvgNorm);
    fprintf(out, "# insertAvg: %lf\n", theAssembly->insertAvgSum/theAssembly->insertAvgNorm);
    fprintf(out, "# kmerAvg: %lf\n", theAssembly->kmerAvgSum/theAssembly->kmerAvgNorm);
    fprintf(out, "# depthScoreAvg: %lf\n", theAssembly->depthScoreAvgSum/theAssembly->depthScoreAvgNorm);
    fprintf(out, "# depthAvg: %lf\n", theAssembly->depthAvgSum/theAssembly->depthAvgNorm);
    fprintf(out, "# totalReads: %ld\n", theAssembly->totalReads);
    fprintf(out, "# totalMappedReads: %ld\n", theAssembly->totalMappedReads);
    fprintf(out, "# totalUnmappedReads: %ld\n", theAssembly->totalUnmappedReads);
    fprintf(out, "# totalPlacedReads: %ld\n", theAssembly->totalPlacedReads);
    fprintf(out, "# readAvgLen: %lf\n", theAssembly->avgReadSize);
    fprintf(out, "# avgReadOverlap: %lf\n", theAssembly->overlapAvgSum/theAssembly->overlapAvgNorm);
    
    if(fullOut == 1){
        for(i = 0; i < theAssembly->numContigs; i++){
            contig_t *contig = theAssembly->contigs[i];
            fprintf(out, "# Reference: %s %i %lf\n# contig position depth ln(depthLike) ln(placeLike) ln(insertLike) ln(kmerLike)\n", contig->name, contig->seqLen, getContigAvgDepth(contig));
            for(j = 0; j < contig->seqLen; j++){
                fprintf(out, "%d %d %0.3f %0.3f %0.3f %0.3f %0.3f\n", i, j, contig->depth[j], contig->depthLogLikelihood[j], contig->matchLogLikelihood[j], contig->insertLogLikelihood[j], contig->kmerLogLikelihood[j]);
            }
        }
    }
    printf("Total ALE Score: %lf\n", theAssembly->totalScore);
}

int assemblySanityCheck(assemblyT *theAssembly){
    int i, j, num = theAssembly->numContigs;
    int error = 1;
    int countAmbiguous = 0;
    for(j=0; j < num ; j++){
        contig_t *contig = theAssembly->contigs[j];
        for(i = 0; i < contig->seqLen; i++){
            if(contig->seq[i] != 'A' && contig->seq[i] != 'T' && contig->seq[i] != 'C' && contig->seq[i] != 'G' && contig->seq[i] != 'N'){
                //printf("Found an ambiguous base in the assembly, contig %d, position %d = %c\n", j, i, contig->seq[i]);
                //contig->seq[i] = 'N';
                error = 0;
                countAmbiguous++;
                int c = contig->ambiguousBaseCount;
                contig->ambiguousBasePositions = (int32_t*) realloc(contig->ambiguousBasePositions, sizeof(*contig->ambiguousBasePositions) * kroundup32(c));
                contig->ambiguousBasePositions[ contig->ambiguousBaseCount++ ] = i;
            }
        }
    }
    if(error == 0){
      printf("Found %d ambiguous bases (excluding N) in the assembly.\n", countAmbiguous);
    }
    return error;
}

assemblyT *loadAssembly(char *filename) {

    // attempt to open the input file
    gzFile assemblyFile = gzopen(filename, "r");
    kseq_t *Aseq;
    if(assemblyFile == NULL){
        printf("Error! Could not open assembly file: %s\n", filename);
        exit(1);
    }

    assemblyT *theAssembly = malloc(sizeof(assemblyT));
    if (theAssembly == NULL)
        exit(1);

    Aseq = kseq_init(assemblyFile);

    readAssembly(Aseq, theAssembly);

    kseq_destroy(Aseq);

    //printf("Done reading in assembly.\n");

    //printAssembly(theAssembly);
    assemblySanityCheck(theAssembly);
    gzclose(assemblyFile);

    return theAssembly;
}

void freeContig(contig_t *contig) {
    if (contig == NULL)
        return;
    free(contig->name);
    free(contig->seq);
    free(contig->depth);
    free(contig->seqNum);
    free(contig->ambiguousBasePositions);
    free(contig->matchLogLikelihood);
    free(contig->insertLogLikelihood);
    free(contig->kmerLogLikelihood);
    free(contig->depthLogLikelihood);
    free(contig);
}

void freeAssembly(assemblyT *theAssembly) {
    if (theAssembly != NULL) {
        if (theAssembly->contigs != NULL) {
        	int i;
        	for (i = 0; i < theAssembly->numContigs; i++)
        		freeContig(theAssembly->contigs[i]);
        	free(theAssembly->contigs);
        }
        free(theAssembly);
    }
}

samfile_t *openSamOrBam(const char *fileName) {
    samfile_t *in = samopen(fileName, "rb", 0);
    if (in == NULL || in->header == NULL) {
        printf("Checking if %s is a SAM formatted file, instead of BAM\n", fileName);
        in = samopen(fileName, "r", 0);
        if (in == NULL || in->header == NULL) {
            printf("Error! Failed to open BAM/SAM file %s\n", fileName);
            exit(1);
        }
    }
    return in;
}

double zNormalizationInsertStd(double std) {
  // int((pmf of normal(0,sigma))^2, 0, infty)
  double expIns = 1.0;
  if (std > 0.0){
      expIns = 1.0/(2.0*sqrt(3.14159265)*std);
  }
  return expIns;
}

int importLibraryParameters(libraryParametersT *libParams, char paramFile[256]){
  FILE *in = fopen(paramFile, "r");
  int temp;
  if(in == NULL){
    printf("Could not open parameter file: %s\n", paramFile);
    libParams = NULL;
    return 0;
  }else{
    
    // read in the file the same way we outputed it
    char comments[256];
    int j;
    temp = fscanf(in, "%s", comments);
    for(j=0; j < MATE_ORIENTATION_MAX; j++) {
        temp = fscanf(in, "%lf", &libParams->mateParameters[j].insertLength);
        temp = fscanf(in, "%lf", &libParams->mateParameters[j].insertStd);
        temp = fscanf(in, "%lf", &libParams->mateParameters[j].libraryFraction);
        temp = fscanf(in, "%ld", &libParams->mateParameters[j].count);
        temp = fscanf(in, "%ld", &libParams->mateParameters[j].mapped);
        // placed is always calculated
        temp = fscanf(in, "%d", &libParams->mateParameters[j].isValid);

        // These are not stored, but calculated
        libParams->mateParameters[j].placed = 0;
        if (libParams->mateParameters[j].isValid && j <= MAPPED_PAIRED_ORIENTATION) {
        	libParams->mateParameters[j].zNormalizationInsert = zNormalizationInsertStd(libParams->mateParameters[j].insertStd);
        } else {
        	libParams->mateParameters[j].zNormalizationInsert = 1.0;
        }
    }
    temp = fscanf(in, "%ld", &libParams->avgReadSize);
    temp = fscanf(in, "%ld", &libParams->numReads);
    temp = fscanf(in, "%lf", &libParams->totalValidSingleFraction);
    temp = fscanf(in, "%lf", &libParams->totalValidMateFraction);
    temp = fscanf(in, "%lf", &libParams->totalChimerMateFraction);
    temp = fscanf(in, "%lf", &libParams->totalUnmappedFraction);
    temp = fscanf(in, "%d", &libParams->qOff);
    temp = fscanf(in, "%d", &libParams->isSortedByName);
    temp = fscanf(in, "%u", &libParams->primaryOrientation);
    fclose(in);
    printf("Read in library parameters.\n");
    return 1;
  }
}

void saveLibraryParameters(libraryParametersT *libParams, char aleFile[256]){
  char paramFile[256];
  int j;
  strcpy(paramFile, aleFile);
  strcat(paramFile, ".param");
  FILE *out = fopen(paramFile, "w");
  fprintf(out, "#ALE_library_parameter_file,_see_doc\n");
  for(j=0; j < MATE_ORIENTATION_MAX; j++) {
      fprintf(out, "%lf\n", libParams->mateParameters[j].insertLength);
      fprintf(out, "%lf\n", libParams->mateParameters[j].insertStd);
      fprintf(out, "%lf\n", libParams->mateParameters[j].libraryFraction);
      fprintf(out, "%ld\n", libParams->mateParameters[j].count);
      fprintf(out, "%ld\n", libParams->mateParameters[j].mapped);
      // placed is always calculated
      fprintf(out, "%d\n", libParams->mateParameters[j].isValid);
  }
  fprintf(out, "%ld\n", libParams->avgReadSize);
  fprintf(out, "%ld\n", libParams->numReads);
  fprintf(out, "%lf\n", libParams->totalValidSingleFraction);
  fprintf(out, "%lf\n", libParams->totalValidMateFraction);
  fprintf(out, "%lf\n", libParams->totalChimerMateFraction);
  fprintf(out, "%lf\n", libParams->totalUnmappedFraction);
  fprintf(out, "%d\n", libParams->qOff);
  fprintf(out, "%d\n", libParams->isSortedByName);
  fprintf(out, "%d\n", libParams->primaryOrientation);
  fclose(out);
  printf("Saved library parameters to %s\n", paramFile);
}

libraryParametersT *computeLibraryParameters(samfile_t *ins, double outlierFraction, int qOff) {

  int i,j;
  long *mapLens[MATE_ORIENTATION_MAX];

  // allocate memory
  for(i=0; i < MATE_ORIENTATION_MAX; i++)
    mapLens[i] = malloc(sizeof(long)*mapLens_MAX);
  libraryParametersT *libParams = malloc(sizeof(libraryParametersT));

  // initialize
  for(j=0; j < MATE_ORIENTATION_MAX; j++) {
    for(i = 0; i < mapLens_MAX; i++)
      mapLens[j][i] = 0;
    libParams->mateParameters[j].insertLength = 0.0;
    libParams->mateParameters[j].insertStd = 0.0;
    libParams->mateParameters[j].zNormalizationInsert = 1.0;
    libParams->mateParameters[j].count = 0;
    libParams->mateParameters[j].mapped = 0;
    libParams->mateParameters[j].placed = 0;
    libParams->mateParameters[j].isValid = 0;
  }
  // set SINGLE_READ insert parameters
  libParams->mateParameters[SINGLE_READ].insertStd = 1.0;
  libParams->mateParameters[SINGLE_READ].zNormalizationInsert = zNormalizationInsertStd(1.0);

  libParams->qOff = qOff;
  libParams->avgReadSize = 0;
  libParams->numReads = 0;
  libParams->isSortedByName = -1; // undefined
  libParams->primaryOrientation = SINGLE_READ;  // can be only SINGLE_READ or the most abundant of the MAPPED_PAIRED_ORIENTATION


  bam1_t *thisRead = bam_init1();
  long readCount = 0;
  long unmappedReads = 0;
  long chimericReads = 0;
  long improperReads = 0;
  long mappedPairedReads = 0;
  long newNames = 0;
  long newPairedNames = 0;
  long pairedReads = 0;
  char *lastName = strdup("");

  while(1){
    enum MATE_ORIENTATION orientation = readNextBAM(ins, thisRead);

    if (orientation == NO_READS){
      break;
    }

    libraryMateParametersT *mateParams = &libParams->mateParameters[orientation];

    int mapLen = -1;
    switch(orientation) {
      case(VALID_FR):
      case(VALID_RF):
      case(VALID_FF):
      case(NOT_PROPER_FR):
      case(NOT_PROPER_RF):
      case(NOT_PROPER_FF):
        mapLen = getFragmentMapLenBAM(thisRead);
        mappedPairedReads++;
        break;

      case(SINGLE_READ):
      case(UNMAPPED_SINGLE):
        break;

      case(SINGLE_UNMAPPED_MATE):
        break;

      case(CHIMER):
        chimericReads++;
        break;

      case(UNMAPPED_PAIR):
      case(UNMAPPED_MATE):
        break;

      default:
        // should never get here
        printf("Improper read with invalid orientation/pairing: %s\n", lastName);
        improperReads++;
    }

    readCount++;
    mateParams->count++;
    libParams->avgReadSize += thisRead->core.l_qseq;
    libParams->numReads++;
    if (libParams->qOff < 0) {
      libParams->qOff = guessQualityOffset(thisRead);
    }

    if ( (thisRead->core.flag & BAM_FUNMAP) == BAM_FUNMAP ) {
        unmappedReads++;
    } else {
    	mateParams->mapped++;
    }

    if ( (thisRead->core.flag & BAM_FPAIRED) == BAM_FPAIRED ) {
    	pairedReads++;
    }
    if (strcmp(lastName, bam1_qname(thisRead)) != 0){
        newNames++;
        if ( (thisRead->core.flag & BAM_FPAIRED) == BAM_FPAIRED ) {
        	newPairedNames++;
        }
    }
    free(lastName);
    lastName = strdup(bam1_qname(thisRead));

    if (mapLen > 0) {
      mateParams->insertLength += mapLen;
      if (mapLen < mapLens_MAX){
        ++mapLens[orientation][mapLen];
      }
    }

    if (readCount%1000000 == 0){
      printf("Read %ld reads...\n", readCount);
    }
  }
  bam_destroy1(thisRead);
  free(lastName);

  if (pairedReads > 0 && newNames <= readCount/2 && newPairedNames == pairedReads / 2) {
    libParams->isSortedByName = 1;
    printf("Setting library to be sorted by name (%ld new sequential names vs %ld reads)\n", newNames, readCount);
  } else {
    libParams->isSortedByName = 0;
  }

  // zero out top and bottom outliers
  long totalValidMateReads = 0;
  long totalValidSingleReads = 0;
  double maximumFraction = 0.0;
  double pairedReadsFraction = (double) mappedPairedReads / (double) libParams->numReads;
  for(j = 0; j < MATE_ORIENTATION_MAX; j++) {
    libraryMateParametersT *mateParams = &libParams->mateParameters[j];
    //printf("Evaluating %s orientation with %ld reads\n", MATE_ORIENTATION_LABELS[j], mateParams->count);
    long observed = 0;
    long purged = 0;
    long lengthTotal = 0;
    for(i = 0; i < mapLens_MAX; i++) {
      if (observed > (1.0-outlierFraction) * mateParams->count) {
        purged += mapLens[j][i];
        mapLens[j][i] = 0;
      } else {
        observed += mapLens[j][i];
        lengthTotal += i*mapLens[j][i];
      }
      if (observed < outlierFraction * mateParams->count) {
        purged += mapLens[j][i];
        mapLens[j][i] = 0;
      }
    }
    mateParams->libraryFraction = (double) mateParams->count / (double) libParams->numReads;
    if (mateParams->libraryFraction > maximumFraction && j <= MAPPED_PAIRED_ORIENTATION) {
        maximumFraction = mateParams->libraryFraction;
        libParams->primaryOrientation = j;
    }

    if (mateParams->count > 0 && (j == VALID_FR || j == VALID_RF || j == VALID_FF || j == NOT_PROPER_FR || j == NOT_PROPER_RF || j == NOT_PROPER_FF )) {

      // TODO better test for significance and normal distribution for a valid orientation
      if (mateParams->libraryFraction > SIGNIFICANT_LIBRARY_FRACTION * pairedReadsFraction ) {
        mateParams->isValid = 1;
      } else {
        improperReads += mateParams->count;
        mateParams->isValid = 0;
      }

    } else if (mateParams->count > 0 && (j == SINGLE_READ || j == SINGLE_UNMAPPED_MATE)) {
        totalValidSingleReads += mateParams->count;
        mateParams->isValid = 1;
    } else {
        improperReads += mateParams->count;
        mateParams->isValid = 0;
    }

    if (mateParams->isValid == 1) {
      /*printf("Read %ld %s oriented sequences and purged %ld %0.1lf%% & %0.1lf%% outliers.  This %s a valid orientation (%01lf%%)\n",
          mateParams->count,
          MATE_ORIENTATION_LABELS[j],
          purged,
          (outlierFraction*100.0),
          ((1.0-outlierFraction)*100.0),
          mateParams->isValid ? "is" : "is NOT",
          mateParams->libraryFraction*100.0);*/
    }

    long modifiedReadCount = mateParams->count - purged;

    if (mateParams->isValid && j <= MAPPED_PAIRED_ORIENTATION) {
      mateParams->insertLength = (double)lengthTotal/(double)modifiedReadCount;
      for(i = 0; i < mapLens_MAX; i++){
        if(mapLens[j][i] > 0){
          double tmp = mapLens[j][i]*((double)i - mateParams->insertLength)*((double)i - mateParams->insertLength);
          mateParams->insertStd += tmp;
          //printf("i : %s mapLens[i] :: %i : %ld\n", MATE_ORIENTATION_LABELS[j], i, mapLens[j][i]);
        }
      }
      mateParams->insertStd = sqrt(mateParams->insertStd/(double)(modifiedReadCount-1));
      mateParams->zNormalizationInsert = zNormalizationInsertStd(mateParams->insertStd);
      printf("Found %s sample avg insert length to be %lf from %ld mapped reads\n", MATE_ORIENTATION_LABELS[j], mateParams->insertLength, modifiedReadCount);
      printf("Found %s sample insert length std to be %lf\n", MATE_ORIENTATION_LABELS[j], mateParams->insertStd);
      if (mateParams->insertStd > MAXIMUM_STD_DEV_FRACTION * mateParams->insertLength) {
    	  mateParams->isValid = 0;
    	  improperReads += mateParams->count;
    	  printf("stddev is too large (> %lf * insertLength).  This orientation is invalid\n", MAXIMUM_STD_DEV_FRACTION );
      } else {
          totalValidMateReads += mateParams->count;
      }
    }
  }
  printf("There were %ld total reads, %ld paired (%ld properly mated), %ld proper singles, %ld improper reads (%ld chimeric). (%ld reads were unmapped)\n", readCount, pairedReads, totalValidMateReads, totalValidSingleReads, improperReads, chimericReads, unmappedReads);

  libParams->avgReadSize = libParams->avgReadSize / libParams->numReads;

  libParams->totalValidMateFraction = (double) (totalValidMateReads) / (double) libParams->numReads;
  libParams->totalValidSingleFraction = (double) totalValidSingleReads / (double) libParams->numReads;
  libParams->totalChimerMateFraction = (double) (improperReads) / (double) libParams->numReads;
  libParams->totalUnmappedFraction = (double) unmappedReads / (double) libParams->numReads;

  /*printf("ValidMates %0.3lf%%, Single(+half-mapped mate pairs) %0.3lf%%, Improper/ChimerMates %0.3lf%% (Unmapped %0.3lf%%)\n",
          libParams->totalValidMateFraction*100,
          libParams->totalValidSingleFraction*100,
          libParams->totalChimerMateFraction*100,
          libParams->totalUnmappedFraction*100);*/

  // release memory
  for(i=0; i < MATE_ORIENTATION_MAX; i++)
    free(mapLens[i]);

  return libParams;
}

