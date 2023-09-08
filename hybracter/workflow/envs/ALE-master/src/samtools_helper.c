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


#include "samtools_helper.h"

// only slight modifications from function found in bam_md.c
void bam_fillmd1_core_ALE(bam1_t *b, char *ref)
{
	int flag = (UPDATE_MD | UPDATE_NM), max_nm = 0;
	uint8_t *seq = bam1_seq(b);
	uint32_t *cigar = bam1_cigar(b);
	bam1_core_t *c = &b->core;
	int i, x, y, u = 0;
	kstring_t *str;
	int32_t old_nm_i = -1, nm = 0;

	str = (kstring_t*)calloc(1, sizeof(kstring_t));
	for (i = y = 0, x = c->pos; i < c->n_cigar; ++i) {
		int j, l = cigar[i]>>4, op = cigar[i]&0xf;
		if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
			for (j = 0; j < l; ++j) {
				int z = y + j;
				int c1 = bam1_seqi(seq, z), c2 = bam_nt16_table[(int)ref[x+j]];
				if (ref[x+j] == 0) break; // out of boundary
				if ((c1 == c2 && c1 != 15 && c2 != 15) || c1 == 0) { // a match
					if (flag&USE_EQUAL) seq[z/2] &= (z&1)? 0xf0 : 0x0f;
					++u;
				} else {
					kputw(u, str); kputc(ref[x+j], str);
					u = 0; ++nm;
				}
			}
			if (j < l) break;
			x += l; y += l;
		} else if (op == BAM_CDEL) {
			kputw(u, str); kputc('^', str);
			for (j = 0; j < l; ++j) {
				if (ref[x+j] == 0) break;
				kputc(ref[x+j], str);
			}
			u = 0;
			if (j < l) break;
			x += l; nm += l;
		} else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
			y += l;
			if (op == BAM_CINS) nm += l;
		} else if (op == BAM_CREF_SKIP) {
			x += l;
		}
	}
	kputw(u, str);
	// apply max_nm
	if (max_nm > 0 && nm >= max_nm) {
		for (i = y = 0, x = c->pos; i < c->n_cigar; ++i) {
			int j, l = cigar[i]>>4, op = cigar[i]&0xf;
			if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
				for (j = 0; j < l; ++j) {
					int z = y + j;
					int c1 = bam1_seqi(seq, z), c2 = bam_nt16_table[(int)ref[x+j]];
					if (ref[x+j] == 0) break; // out of boundary
					if ((c1 == c2 && c1 != 15 && c2 != 15) || c1 == 0) { // a match
						seq[z/2] |= (z&1)? 0x0f : 0xf0;
						bam1_qual(b)[z] = 0;
					}
				}
				if (j < l) break;
				x += l; y += l;
			} else if (op == BAM_CDEL || op == BAM_CREF_SKIP) x += l;
			else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) y += l;
		}
	}
	// update NM
	if (flag & UPDATE_NM) {
		uint8_t *old_nm = bam_aux_get(b, "NM");
		if (c->flag & BAM_FUNMAP) return;
		if (old_nm) old_nm_i = bam_aux2i(old_nm);
		if (!old_nm) bam_aux_append(b, "NM", 'i', 4, (uint8_t*)&nm);
		else if (nm != old_nm_i) {
			//fprintf(stderr, "[bam_fillmd1] different NM for read '%s': %d -> %d\n", bam1_qname(b), old_nm_i, nm);
			bam_aux_del(b, old_nm);
			bam_aux_append(b, "NM", 'i', 4, (uint8_t*)&nm);
		}
	}
	// update MD
	if (flag & UPDATE_MD) {
		uint8_t *old_md = bam_aux_get(b, "MD");
		if (c->flag & BAM_FUNMAP) return;
		if (!old_md) bam_aux_append(b, "MD", 'Z', str->l + 1, (uint8_t*)str->s);
		else {
			int is_diff = 0;
			if (strlen((char*)old_md+1) == str->l) {
				for (i = 0; i < str->l; ++i)
					if (toupper(old_md[i+1]) != toupper(str->s[i]))
						break;
				if (i < str->l) is_diff = 1;
			} else is_diff = 1;
			if (is_diff) {
				//fprintf(stderr, "[bam_fillmd1] different MD for read '%s': '%s' -> '%s'\n", bam1_qname(b), old_md+1, str->s);
				bam_aux_del(b, old_md);
				bam_aux_append(b, "MD", 'Z', str->l + 1, (uint8_t*)str->s);
			}
		}
	}
	// drop all tags but RG
	if (flag&DROP_TAG) {
		uint8_t *q = bam_aux_get(b, "RG");
		bam_aux_drop_other(b, q);
	}
	// reduce the resolution of base quality
	if (flag&BIN_QUAL) {
		uint8_t *qual = bam1_qual(b);
		for (i = 0; i < b->core.l_qseq; ++i)
			if (qual[i] >= 3) qual[i] = qual[i]/10*10 + 7;
	}
	free(str->s); free(str);
}

