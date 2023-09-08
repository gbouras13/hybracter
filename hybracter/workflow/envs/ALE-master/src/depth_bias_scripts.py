#!/usr/bin/python

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
import matplotlib.pylab as plt # for plotting
import scipy.stats
from scipy.optimize import fmin
from matplotlib.backends.backend_pdf import PdfPages
import sys
#import time # for timing
#import commands # print commands.getoutput(script)

import plotter3

def get_first_contig_seq(fasta_file):
    """Returns the sequence (in a list) of the first contig in a fasta file"""
    file_in = open(fasta_file, 'r')
    started = False
    seq = []
    for line in file_in:
        if line[0] == '>': # new contig
            if not started:
                started = True
            else:
                break
        else: # contig sequence
            seq.extend(list(line[:-1])) # add it to the sequence
    return seq

def get_GC_content(seq, window):
    GC_content_single_vec = numpy.zeros(len(seq) - window)
    GC_total = 0
    for i, base in enumerate(seq):
        if not i%100000: print i
        if i < window:
            if base in ['G','C','g','c']:
                GC_total += 1
        else:
            GC_content_single_vec[i - window] = float(GC_total)/float(window)
            if base in ['G','C','g','c']:
                GC_total += 1
            if seq[i - window] in ['G','C','g','c']:
                GC_total -= 1
    #print GC_content_single_vec

    GC_content_window_vec = numpy.zeros(len(seq))

    for pos, val in enumerate(GC_content_single_vec):
        if not pos%100000: print pos
        if pos == 0:
            GC_content_window_vec[pos] = GC_content_single_vec[0]
        elif pos < window:
            GC_content_window_vec[pos] = numpy.sum(GC_content_single_vec[:pos+1])/float(pos+1)
        else:
            GC_content_window_vec[pos] = numpy.sum(GC_content_single_vec[pos-window:pos+1])/float(window)

    # get the end
    for pos in range(len(seq) - window, len(seq)):
        GC_content_window_vec[pos] = numpy.sum(GC_content_single_vec[pos-window:])/float(len(seq) - pos)
    return GC_content_window_vec

def gamma_poisson_opt2(in_vec, depth_dist, error_type="MLE"):
    p = in_vec[0]
    r = in_vec[1]
    error = 0.0
    neg_binom_dist = numpy.zeros(len(depth_dist))
    for k in range(len(depth_dist)):
        neg_binom_dist[k] = scipy.stats.distributions.nbinom.pmf(k, r, p)
        if error_type == "MLE":
            error -= depth_dist[k]*numpy.log(neg_binom_dist[k])
    return error

def gamma_poisson_opt(p, r, depth_dist, error_type="MLE"):
    error = 0.0
    neg_binom_dist = numpy.zeros(len(depth_dist))
    for k in range(len(depth_dist)):
        neg_binom_dist[k] = scipy.stats.distributions.nbinom.pmf(k, r, p)
        if error_type == "MLE":
            error -= depth_dist[k]*numpy.log(neg_binom_dist[k])
    return error

def basic_poisson_opt(poisson_mean, depth_dist, error_type="MLE"):
    error = 0.0
    poisson_dist = numpy.zeros(len(depth_dist))
    for k in range(len(depth_dist)):
        poisson_dist[k] = scipy.stats.distributions.poisson.pmf(k, poisson_mean)
        if error_type == "MLE":
            error -= depth_dist[k]*numpy.log(poisson_dist[k])
    return error

def get_histogram(input_x, input_y):
    max_val = numpy.max(input_x)
    min_val = numpy.min(input_x)
    bin_size = (max_val - min_val)/100.0
    x_vec = numpy.arange(min_val, max_val, bin_size)
    if bin_size < 1e-6:
        bin_size = 0.0001
        x_vec = numpy.arange(0.0, 0.01, bin_size)
    histogram = numpy.zeros(len(x_vec))
    for i, value in enumerate(input_x):
        try:
            histogram[min((max((0,int(numpy.floor((value - min_val)/bin_size)))),99))] += input_y[i]
        except:
            pass
    
    return x_vec, histogram

def plot_depth_dists(GC_content_window_vec, depth_vec, GC_content=0.5, width=0.01, error_type="MLE"):
    depth_dist = numpy.zeros(numpy.max(depth_vec)+1)
    depth_set = []
    depth_norm = 0.0

    fig = plt.figure(figsize=(12,8))
    ax1 = fig.add_subplot(221)
    ax1.set_title("Depth distributions")
    ax1.set_xlabel("Depth")

    ax4 = fig.add_subplot(222)
    ax4.set_title("Depth distributions")
    ax4.set_xlabel("Depth")


    # depth distribution
    for pos, GC_cont in enumerate(GC_content_window_vec):
        if numpy.floor(GC_content/width) == numpy.floor(GC_cont/width):
            depth_set.append(depth_vec[pos])
            depth_dist[depth_vec[pos]] += 1
            depth_norm += 1.0
    depth_dist = depth_dist/depth_norm
    ax1.plot(depth_dist)
    ax4.plot(depth_dist)

    # for basic poisson
    poisson_mean = numpy.mean(depth_set)
    poisson_var = pow(numpy.std(depth_set), 2.0)
    poisson_dist = numpy.zeros(numpy.max(depth_vec)+1)
    for k in range(len(depth_dist)):
        poisson_dist[k] = scipy.stats.distributions.poisson.pmf(k, poisson_mean)
    ax1.plot(poisson_dist)
    ax4.plot(poisson_dist)
    print "depth mean: %f" % poisson_mean
    print "depth var: %f" % poisson_var
    print "Basic poisson error: %s" % basic_poisson_opt(poisson_mean, depth_dist, error_type=error_type)

    # neg binom with moment matching
    nbinom_r = poisson_mean/(poisson_var/poisson_mean - 1.0)
    nbinom_p = poisson_mean/(nbinom_r + poisson_mean)
    print "mm gamma param: %f %f" % (nbinom_r, nbinom_p)
    print "Basic mm gamma error: %s" % gamma_poisson_opt(nbinom_p, nbinom_r, depth_dist, error_type=error_type)
    nbinom_dist_mm = numpy.zeros(len(depth_dist))
    for k in range(len(depth_dist)):
        nbinom_dist_mm[k] = scipy.stats.distributions.nbinom.pmf(k, nbinom_r, nbinom_p)
    ax1.plot(nbinom_dist_mm)
    ax4.plot(nbinom_dist_mm)

    # for fitted neg binom of p (fmin)
    nbinom_p_start = 0.01
    nbinom_r = poisson_mean
    nbinom_p_opt = poisson_mean/(nbinom_r + poisson_mean) # by definition
    print "gamma r-fixed param start: %f" % (nbinom_p_start)
    print "gamma r-fixed param: %f" % (nbinom_p_opt)
    print "Basic r-fixed gamma error: %s" % gamma_poisson_opt(nbinom_p_opt, nbinom_r, depth_dist, error_type=error_type)
    nbinom_dist_opt = numpy.zeros(len(depth_dist))
    for k in range(len(depth_dist)):
        nbinom_dist_opt[k] = scipy.stats.distributions.nbinom.pmf(k, nbinom_r, nbinom_p_opt)
    ax1.plot(nbinom_dist_opt)
    ax4.plot(nbinom_dist_opt)

    # for fitted neg binom or both r and p (fmin)
    nbinom_p_start = 0.01
    nbinom_r_start = poisson_mean
    [nbinom_p_opt, nbinom_r_opt] = fmin(gamma_poisson_opt2, (nbinom_p_start, nbinom_r_start), args = [depth_dist, error_type])
    print "gamma free param start: %f %f" % (nbinom_r_start, nbinom_p_start)
    print "gamma free param: %f %f" % (nbinom_r_opt, nbinom_p_opt)
    print "double free gamma error: %s" % gamma_poisson_opt2([nbinom_p_opt, nbinom_r_opt], depth_dist, error_type=error_type)
    nbinom_dist_opt2 = numpy.zeros(len(depth_dist))
    for k in range(len(depth_dist)):
        nbinom_dist_opt2[k] = scipy.stats.distributions.nbinom.pmf(k, nbinom_r_opt, nbinom_p_opt)
    ax1.plot(nbinom_dist_opt2)
    ax4.plot(nbinom_dist_opt2)

    

    # differences
    #ax2 = fig.add_subplot(222)
    #ax2.set_title("Difference between analytic and emperic")
    #ax2.set_xlabel("Depth")
    #ax2.plot(depth_dist - depth_dist)
    #ax2.plot(poisson_dist - numpy.log(depth_dist))
    #ax2.plot(nbinom_dist_opt - numpy.log(depth_dist))
    #ax2.plot(nbinom_dist_opt2 - numpy.log(depth_dist))

    # distribution of scores
    ax3 = fig.add_subplot(223)
    ax3.set_title("Distribution of scores")
    ax3.set_xlabel("Score")
    x, y = get_histogram(poisson_dist, depth_dist)
    ax3.plot(x, numpy.zeros(len(x))) # blank one to keep the colors right
    ax3.plot(x, y)
    x, y = get_histogram(nbinom_dist_mm, depth_dist)
    ax3.plot(x, y)
    x, y = get_histogram(nbinom_dist_opt, depth_dist)
    ax3.plot(x, y)
    x, y = get_histogram(nbinom_dist_opt2, depth_dist)
    ax3.plot(x, y)
    leg = ax1.legend(('Actual Depth', 'Poisson', 'Neg Binom mm', 'Neg Binom (fixed r)', 'Neg Binom'), 'upper left')

    return fig

def plot_depth_vs_pos(depth_vec):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title("total depth distribution")
    ax.plot(depth_vec)

    return fig

def plot_total_depth_dist(depth_vec):
    depth_dist = numpy.zeros(numpy.ceil(numpy.max(depth_vec) + 1))
    depth_norm = float(len(depth_vec))
    for depth in depth_vec:
        depth_dist[numpy.floor(depth)] += 1.0
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title("depth vs position")
    ax.plot(depth_dist/depth_norm)

    return fig

def plot_GC_dist(GC_content_window_vec, width=0.01):
    GC_dist = numpy.zeros(1.0/width)
    GC_norm = float(len(GC_content_window_vec))
    for GC_cont in GC_content_window_vec:
        GC_dist[numpy.floor(GC_cont/width)] += 1.0
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_title("GC content distribution")
    ax.plot(GC_dist/GC_norm)

    return fig

def get_depth_vector(placement_file):
    """Returns the depth vector of the first contig of an ale file"""
    contigs = plotter3.read_in_info(placement_file)
    # only return the depth vector of the first contig
    return contigs[0].depth

def usage(exit_value):
    print >> sys.stderr, "depth_bias_scripts.py window_len ale_file seq_file"
    sys.exit(exit_value)

def main():
    import os
    if len(sys.argv) == 1: usage(1)
    if sys.argv[1] == '--help' or sys.argv[1] == '-h': usage(0)
    if len(sys.argv) > 4: usage(1)
    try:
        window_len = int(sys.argv[1])
    except ValueError: 
        print >> sys.stderr, 'invalid literal <%s> for window len' %(sys.argv[1])
        sys.exit(1)
    ale_file = sys.argv[2]
    seq_file = sys.argv[3]
    if not os.path.isfile(ale_file):
        print >> sys.stderr, 'no such file or directory', ale_file
        sys.exit(1)
    if not os.path.isfile(seq_file):
        print >> sys.stderr, 'no such file or directory', seq_file
        sys.exit(1)
    depth_vec = get_depth_vector(ale_file)
    seq = get_first_contig_seq(seq_file)
    GC_content_window_vec = get_GC_content(seq, window_len)

    figure_name = ale_file + 'GC_cont.pdf'
    pdf_stream = PdfPages(figure_name)

    plot_depth_vs_pos(depth_vec)
    pdf_stream.savefig()

    plot_total_depth_dist(depth_vec)
    pdf_stream.savefig()

    plot_GC_dist(GC_content_window_vec, width=0.01)
    pdf_stream.savefig()

    for GC_content in numpy.arange(0.3,0.6,0.1):
        print "Plotting for content %f (+.1)" % GC_content
        plot_depth_dists(GC_content_window_vec, depth_vec, GC_content=GC_content, width=0.1)
        pdf_stream.savefig()

    pdf_stream.close()
    print "File saved in %s" % figure_name

if __name__ == '__main__':
    main()

