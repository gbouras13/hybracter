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

import numpy # for sci comp
#import matplotlib.pylab as plt # for plotting
#import sys
import heapq

import ProgressBar

#import time # for timing
#import commands # print commands.getoutput(script)

import CmdIn

USAGE = """Usage:
./ALE_low_scores.py [options] input.ale

basic options:
-h : print full usage and options"""

FULL_USAGE = """Usage:
./ALE_low_scores.py [options] input.ale

help options:
-h   : print full usage and options

basic options [default]:
-st  : score types; (p)lacement, (i)nsert, (d)epth and/or (k)mer [p]
-nwp : n worst positions (the number of positions reported) [100]
-idb : ignore depth below a certain threshold (hard stops) [2]
-mcs : min contig size for output [100]
-vf  : variants file (for ROC stats) [None]
-op  : output prefix for file output [input.ale]
-rl  : read length (for ROC stats) [77]
-il  : insert length (for ROC stats) [200]

binary options [default]:
-ie  : ignore edges [off]
-v   : verbose [off]
-co  : condense output (for IGV viewing, less human readable) [off]"""

class Contig():
    """A contig from an ALE assembly

    more info

    Attributes:
        name (str): The name of the contig (defaults to unnamed)

        length (int): The length of the contig

        depth (numpy.array): A numpy vector of depths for each position

        depth_prob (numpy.array): A numpy vector of depth probabilities for each position

        placement_prob (numpy.array): A numpy vector of placement probabilities for each position

        insert_prob (numpy.array): A numpy vector of insert probabilities for each position

        kmer_prob (numpy.array): A numpy vector of kmer probabilities for each position

        total_prob (numpy.array): A numpy vector of total probabilities for each position

    """
    def __init__(self, length=0, name="unnamed"):
        """Inits Contig with name and length (and prob vectors)
        
        Kwargs:
            length (int): The length of the contig (>=0)

            name (str): The name of the contig
            
        Raises:
            ValueError: length must be >= 0
        """

        if length < 0:
          raise ValueError("length must be >= 0")

        self.name = name
        self.length = length

        self.depth = numpy.zeros(length)

        self.prob_vecs={'d':LikelihoodVector(name, length=length, color='r'),
                   'p':LikelihoodVector(name, length=length, color='b'),
                   'i':LikelihoodVector(name, length=length, color='m'),
                   'k':LikelihoodVector(name, length=length, color='g')}

    def find_N_worst(self, typer, N, user_param):
        """Finds the N worst scores in each of the likelihood vectors

        uses a set of min heaps to do this in O(length) time
        ignores edges if user_param.get("ignores_edges") is true
        an edge is defined by user_param.get("min_plot_size")/2
        
        returns an array of PositionViolation(s)
        """
        ignore_edges = user_param.get("ignore_edges")
        edge_length = user_param.get("min_contig_size")/2
        ignore_depth_below = user_param.get("ignore_depth_below")
        read_len = user_param.get("read_len")
        insert_len = user_param.get("insert_len")

        worst = []
        # build heapable list
        prob_heap = []
        if ignore_edges:
            prob_to_iterate = self.prob_vecs[typer].prob[edge_length:-edge_length]
        else:
            prob_to_iterate = self.prob_vecs[typer].prob
        extendable = False
        for i, prob in enumerate(prob_to_iterate):
            if self.depth[i] < ignore_depth_below:
                extendable = False
            else:
                if extendable and prob == prob_heap[-1][0]:
                    prob_heap[-1][2] += 1
                else:
                    prob_heap.append([prob, i, i]) # prob because heapq uses a minheap
                    extendable = True

        heapq.heapify(prob_heap) # make into a heap in O(n) time
        smallest_N = heapq.nsmallest(N, prob_heap) # \^altTheory
        for position in smallest_N:
            tv = ThresholdViolation(typer, position[1], position[2], position[0], self.name)
            tv.within_read_of_hard_stop = self.low_depth_within_len(position[1], position[2], read_len, ignore_depth_below)
            tv.within_insert_of_hard_stop = self.low_depth_within_len(position[1], position[2], insert_len, ignore_depth_below)
            worst.append(tv)
        return worst

    def low_depth_within_len(self, start, end, length, low_val):
        depth_list = self.depth[max((0,start-length)):min((self.length,end+length))]
        if min(depth_list) < low_val:
            return True
        else:
            return False

    def find_and_validate_low_scores(self, N, type_of, igv_variants, user_params, silent=True):
        """Prints out the lowest score sections of a contig and matches them
           up with igv variants from a file, if they exist"""

        # get the worst scores of type_of using a heap and group them
        worst_scores = self.find_N_worst(type_of, N, user_params)
        seeds = group_threshold_violations(worst_scores)
        if len(seeds) == 0:
            return seeds, 0.0, 0.0, 0.0 

        c_exact = 0
        c_read = 0
        c_insert = 0

        # calculate false positive statistics
        if igv_variants:
            for i, seed in enumerate(seeds):

                f_exact = False
                f_read = False
                f_insert = False

                if not silent:
                    print "Group %d" % i
                # apply variants to all the violations in the group
                for vi in seed:
                    if igv_variants and self.name in igv_variants:
                        vi.apply_variants(igv_variants[self.name], user_params)
                    
                    if not silent:
                        print str(vi)

                    # count
                    if vi.igv_variants:
                        f_exact = True
                    if vi.within_read_of_error or vi.within_read_of_hard_stop:
                        f_read = True
                    if vi.within_insert_of_error or vi.within_insert_of_hard_stop:
                        f_insert = True

                c_exact += int(f_exact)
                c_read += int(f_read)
                c_insert += int(f_insert)


        return seeds, float(c_exact)/float(len(seeds)), float(c_read)/float(len(seeds)), float(c_insert)/float(len(seeds))

class LikelihoodVector(object):
    """A likelihood vector of some type in a contig"""
    def __init__(self, contig_name, length=None, color=None):
        if not length or length < 0:
            raise ValueError("length must be >= 0")
        if not color:
            raise ValueError("need to specify color")
        self.contig_name = contig_name
        self.length = length
        self.color = color
        self.prob = numpy.zeros(length)
        self.prob_smoothed = numpy.zeros(length)
        self.smoothing_width = 1

        self.thresh_main_mean = None
        self.thresh_main_std = None

class IGVVariant(object):
    """ A variant from an .igv file"""
    def __init__(self, contig_name, position, type_of):
        self.contig_name = contig_name
        self.position = position - 1 # 0 based distance
        self.type_of = type_of
        self.called_exactly = False
        self.called_within_read = False
        self.called_within_insert = False

    def __str__(self):
        return "%s\t%d\t%d\t%s" % (self.contig_name, self.position+1, self.position+1, self.type_of)

class ThresholdViolation(object):
    """A threshold of an ALE score vector"""
    def __init__(self, type_of, start, end, score, contig_name):
        self.contig_name = contig_name
        self.type_of = type_of
        self.start = start
        self.end = end
        self.score = score
        self.igv_variants = []
        self.within_read_of_error = False
        self.within_insert_of_error = False
        self.within_read_of_hard_stop = False
        self.within_insert_of_hard_stop = False

    def __str__(self):
        return_str = "%s:%d-%d\t%s\t%lf" % (self.contig_name, self.start+1, self.end+1, self.type_of, self.score)
        if self.igv_variants:
            for var in self.igv_variants:
                return_str += " *** %s" % str(var)
        else:
            return_str += " *** nearby: %d %d %d %d" % (int(self.within_read_of_error), int(self.within_insert_of_error), int(self.within_read_of_hard_stop), int(self.within_insert_of_hard_stop))
        return return_str

    def near_error(self, distance, igv_error):
        if self.start - distance <= igv_error.position and self.end + distance >= igv_error.position:
            return True
        return False

    def contains_error(self, igv_error):
        if "INS" in igv_error.type_of: # can be off by 1
            if self.start - 1 <= igv_error.position and self.end >= igv_error.position:
                return True
        else:
            if self.start <= igv_error.position and self.end >= igv_error.position:
                return True
        return False

    def apply_variants(self, variants, user_params):
        read_len = user_params.get("read_len")
        insert_len = user_params.get("insert_len")

        for var in variants:
            assert(var.contig_name == self.contig_name)
            if self.contains_error(var):
                self.igv_variants.append(var)
                var.called_exactly = True
                self.within_read_of_error = True
                var.called_within_read = True
                self.within_insert_of_error = True
                var.called_within_insert = True
            if self.near_error(read_len, var):
                var.called_within_read = True
                if not self.within_read_of_error:
                    self.within_read_of_error = True
            if self.near_error(insert_len, var):
                var.called_within_insert = True
                if not self.within_insert_of_error:
                    self.within_insert_of_error = True



def read_igv_file(file_name, silent=True):
    """reads an igv file and returns a dictionary of IGVErrors keyed on contig_name"""
    fin = open(file_name, "r")
    fin.readline() # skip labels
    igv_variants = {}
    for line in fin:
        contig_name = line.split('\t')[0].strip()
        position = int(line.split('\t')[1].strip())
        type_of = line.split(str(position))[-1].strip()
        if contig_name not in igv_variants:
            igv_variants[contig_name] = []
        igv_variants[contig_name].append(IGVVariant(contig_name, position, type_of))

    if not silent:
        for contig_name in igv_variants:
            for var in igv_variants[contig_name]:
                print str(var)

    return igv_variants

def get_variant_stats(contig_name, igv_variants, silent=True):
    if igv_variants and contig_name in igv_variants:
        contig_variants = igv_variants[contig_name]
    else:
        return 0, 0, 0
    num_var = len(contig_variants)
    num_hit = 0
    num_win_read = 0
    num_win_insert = 0
    for var in contig_variants:
        num_hit += int(var.called_exactly)
        num_win_read += int(var.called_within_read)
        num_win_insert += int(var.called_within_insert)
    if not silent: print "%s: total: %d hit: %d win_read: %d win_insert: %d" % (contig_name, num_var, num_hit, num_win_read, num_win_insert)
    return float(num_hit)/float(num_var), float(num_win_read)/float(num_var), float(num_win_insert)/float(num_var), num_var

def group_threshold_violations(violations, silent=True):
    """Group violations that are within 1 base of each other
    return a array of violation arrays, each one of which contains
    at least one violation, more if they are grouped together.

    input: an array of violations
    output: an array of arrays of violations (grouped by position)
    """
    start_indexed = {}
    end_indexed = {}
    start_to_end = {}
    end_to_start = {}
    for violation in violations:
        start_indexed[violation.start] = violation
        start_to_end[violation.start] = violation.end
        end_indexed[violation.end] = violation
        end_to_start[violation.end] = violation.start

    seeds = []
    while len(violations) > 0:
        seeds.append([violations.pop(0)]) # top violation
        # build bottom
        while True:
            if seeds[-1][0].start - 1 in end_indexed:
                for pos, violation in enumerate(violations):
                    if seeds[-1][0].start - 1 == violation.end:
                        seeds[-1].insert(0, violations.pop(pos)) # delete from violations
                        break
                #del end_indexed[seeds[-1][0].start - 1] # delete from end_index
                #del start_indexed[end_to_start[seeds[-1][0].start - 1]] # delete from start_index 
            else:
                break
        # build top
        while True:
            if seeds[-1][-1].end + 1 in start_indexed:
                for pos, violation in enumerate(violations):
                    if seeds[-1][-1].end + 1 == violation.start:
                        seeds[-1].append(violations.pop(pos)) # delete from violations
                        break
                #del start_indexed[seeds[-1][-1].end + 1] # delete from end_index
                #del end_indexed[start_to_end[seeds[-1][-1].end + 1]] # delete from start_index              
            else:
                break

    if not silent:
        for num, seed in enumerate(seeds):
            print "Violation %d" % num
            for vi in seed:
                print str(vi)

    return seeds

def read_in_info(user_params):
    """Reads in an ALE placement file, returns a list of Contigs.

    Args:
        placement_file: An ALE placement file (*.ale)
            must be in the following format::

                # comments/metadata
                # can have multiple lines, all starting with #
                # Reference: gi|170079663|ref|NC_010473.1| 350000 29.3
                # contig position depth ln(depthLike) ln(placeLike) ln(insertLike) ln(kmerLike) 
                0 0 1.000000 -60.000000 0.194888 -5.760798 -65.565910
                0 1 3.000000 -60.000000 0.466271 -5.608334 -65.142063
                0 2 5.000000 -60.000000 0.010585 -5.541655 -65.531071
                0 3 12.000000 -60.000000 -0.057731 -5.380759 -65.438491

            Specific lines (using the above as an example):
                0. Any number of comment lines starting with #
                1. The length of the contig is::

                       length = int(line.split(' ')[3]) == 350000

                   The name of the contig is::

                       name = line.split(' ')[2] == gi|170079663|ref|NC_010473.1|

                   name **cannot** be 'position'

                2. The following line is ignored and lists what is in the columns of following lines

                3. The data corresponding to the column headers for each position in the contig

                4. See 2.

    Returns:
        A list of Contigs (see class :py:mod:`plotter3.Contig`)

    Raises:
        IOError: An error occured accessing the placement file.

        FormattingError: The placement file was not formatted correctly.
    """

    MINIMUM_VALUE = -120.0 # minimum value we allow a position to have

    placement_file = user_params.get_input("ale_file")
    specific_contig = user_params.get("specific_contig")

    ale_placement_file = open(placement_file, 'r')

    contigs = []
    previous_line_one = ""
    previous_line_two = ""

    for line in ale_placement_file:
        if line[0] == '#':
            if previous_line_one == "":
                previous_line_one = line
            else:
                previous_line_two = previous_line_one
                previous_line_one = line           
        else:
            if previous_line_two != "":
                tName = previous_line_two.split(' ')[2]              
                tLen = int(previous_line_two.split(' ')[-2])
                if specific_contig and len(contigs) > 0 and contigs[-1].name == specific_contig:
                    print "\n"
                    return contigs # found what we are looking for, go to next step
                contigs.append(Contig(tLen, name = tName))
                place = 0
                print "\nReading in contig: " + tName + " len " + str(tLen)
                print ""
                bar = ProgressBar.progressBar(0, tLen, 42)
                previous_line_two = ""
            data = line.split(' ')
            contigs[-1].depth[place] = numpy.double(data[2])
            for i in range(1,7):
                if "-nan"==data[i] or "nan"==data[i] or "inf"==data[i] or "-inf"==data[i] or numpy.double(data[i]) != numpy.double(data[i]):
                    data[i] = MINIMUM_VALUE # Predefined threshold
            contigs[-1].prob_vecs['d'].prob[place] = numpy.double(data[3])
            contigs[-1].prob_vecs['p'].prob[place] = numpy.double(data[4])
            contigs[-1].prob_vecs['i'].prob[place] = numpy.double(data[5])
            contigs[-1].prob_vecs['k'].prob[place] = numpy.double(data[6])           
            place += 1
            if tLen > 40:
                if (place)%(int(tLen)/40) == 0:
                    bar(place)
    print "\n"
    return contigs


def main():
    user_params = CmdIn.CommandLineParameters()
    user_params.__usage__ = USAGE
    user_params.__full_usage__ = FULL_USAGE
    user_params.add_parameter("score_type", "-st", str, 'p') # some combo of pikd
    user_params.add_parameter("specific_contig", "-sc", str, None)
    user_params.add_parameter("N_worst_positions", "-nwp", int, 100)
    user_params.add_parameter("ignore_edges", "-ie", None, False)
    user_params.add_parameter("ignore_depth_below", "-idb", int, 2)
    user_params.add_parameter("min_contig_size", "-mcs", int, 100)
    user_params.add_parameter("variants_file", "-vf", str, None)
    user_params.add_parameter("output_prefix", "-op", str, None)
    user_params.add_parameter("condense_output", "-co", None, False)
    user_params.add_parameter("verbose", "-v", None, False)
    user_params.add_parameter("read_len", "-rl", int, 77)
    user_params.add_parameter("insert_len", "-il", int, 200)
    user_params.add_input("ale_file")

    # read in command line arguments
    user_params.read_sys_args()

    verbose = user_params.get("verbose")

    if user_params.get("min_contig_size") < user_params.get("N_worst_positions"):
        user_params.set("min_contig_size", user_params.get("N_worst_positions"))

    # read in variants file, if it exists
    variants_file = user_params.get("variants_file")
    if variants_file:
        igv_variants = read_igv_file(variants_file)
    else:
        igv_variants = None

    # read in contigs
    contigs = read_in_info(user_params)

    score_type = user_params.get("score_type")
    placements = {}
    group_rates = {}
    variant_rates = {}
    for typer in score_type:
        placements[typer] = {}
        group_rates[typer] = {}
    
    for contig in contigs:
        if contig.length >= user_params.get("min_contig_size"):
            if not user_params.get("specific_contig") or user_params.get("specific_contig") == contig.name:

                
                n = user_params.get("N_worst_positions")
                for typer in score_type:
                    print "Finding %d worst positions (of %d) of type %s in contig %s" % (n, contig.length, typer, contig.name)
                    seeds, g_exact, g_read, g_insert = contig.find_and_validate_low_scores(n, typer, igv_variants, user_params)
                    placements[typer][contig.name] = seeds

                    group_rates[typer][contig.name] = {}
                    group_rates[typer][contig.name]["Exact Match"] = g_exact
                    group_rates[typer][contig.name]["Within Read Length"] = g_read
                    group_rates[typer][contig.name]["Within Insert Length"] = g_insert
                    group_rates[typer][contig.name]["num_groups"] = len(seeds)

                # possibly print out some statistics
                if variants_file:
                    v_exact, v_read, v_insert, v_num = get_variant_stats(contig.name, igv_variants)
                    variant_rates[contig.name] = {}
                    variant_rates[contig.name]["Exact Match"] = v_exact
                    variant_rates[contig.name]["Within Read Length"] = v_read
                    variant_rates[contig.name]["Within Insert Length"] = v_insert
                    variant_rates[contig.name]["total_num"] = v_num

    # make output files...

    output_prefix = user_params.get("output_prefix")
    if not output_prefix:
        output_prefix = user_params.get_input("ale_file")

    fout = open("%s.thresholds" % (output_prefix), 'w')
    for typer in placements:
        tmp_str = "----------------------- TYPE %s -------------------------" % typer
        if verbose: print tmp_str
        fout.write("%s\n" % tmp_str)
        for contig in placements[typer]:
            tmp_str = "->Contig: %s" % contig
            if verbose: print tmp_str
            fout.write("%s\n" % tmp_str) 
            seeds = placements[typer][contig]
            for i, seed in enumerate(seeds):
                tmp_str = "--->Group %d" % i
                if verbose: print tmp_str
                fout.write("%s\n" % tmp_str)
                for vi in seed:
                    tmp_str = "----->%s" % str(vi)
                    if verbose: print tmp_str
                    fout.write("%s\n" % tmp_str)
    # groups

    condense_output = user_params.get("condense_output")

    # .igv files
    for typer in placements:
        fout = open("%s.%s.igv" % (output_prefix, typer), 'w')
        for contig in placements[typer]:
            seeds = placements[typer][contig]
            for seed in seeds:
                if condense_output:
                    tmp_str = "%s:%d-%d\t%s" % (contig, seed[0].start, seed[-1].end, typer)
                    if verbose: print tmp_str
                    fout.write("%s\n" % tmp_str)
                else:
                    for vi in seed:
                        tmp_str = str(vi)
                        if verbose: print tmp_str
                        fout.write("%s\n" % tmp_str)
        fout.close()

    read_len = user_params.get("read_len")
    insert_len = user_params.get("insert_len")
    for typer in score_type:
        for contig in group_rates[typer]:       
            print "%s Type: %s true positive rate (near variant or hard stop) (of %d groups):" % (contig, typer, group_rates[typer][contig]["num_groups"])
            print "Exact Hit: %f" % group_rates[typer][contig]["Exact Match"]
            print "Within Read Length (%d): %f" % (read_len, group_rates[typer][contig]["Within Read Length"])
            print "Within Insert Length (%d): %f" % (insert_len, group_rates[typer][contig]["Within Insert Length"])
    for contig in variant_rates:
        print "%s precision of %d variants:" % (contig, variant_rates[contig]["total_num"])
        print "Exact Hit: %f" % variant_rates[contig]["Exact Match"]
        print "Within Read Length (%d): %f" % (read_len, variant_rates[contig]["Within Read Length"])
        print "Within Insert Length (%d): %f" % (insert_len, variant_rates[contig]["Within Insert Length"])

if __name__ == '__main__':
    main()

