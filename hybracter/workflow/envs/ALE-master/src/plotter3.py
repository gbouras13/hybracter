#!/usr/bin/env python

# (C) 2011 Scott Clark

"""plotter3 - a plotting package for ALE scoring output

    Depends upon:
        matplotlib
            http://matplotlib.sourceforge.net/
        numpy
            http://numpy.scipy.org/
        mpmath
            http://code.google.com/p/mpmath/
        python
            2.6+
            http://www.python.org/
        pymix
            http://www.pymix.org/pymix/
"""

__version__ = "0.2"
__usage__ = """Usage: ./ALE_plotter.py [-options] <inputfile.ale>

where basic options are:
  -h      : show brief help on version and full usage
  -nosave : do not save the figure as a pdf (instead plot to screen)
"""
__full_usage__="""Usage: ./ALE_plotter.py [-options] <inputfile.ale>

where basic options are:
  -h      : show brief help on version and full usage
  -nosave : do not save the figure as a pdf (instead plot to screen)

parameter options accepting <f>loats and <i>ntegers and <s>trings (default):
  -s <i>   : the starting position to plot (0)
  -e <i>   : the ending position of the plot (contig length)
  -spo     : sub plots on, recursive search for errors (off)
  -pt <s>  : plot type 'i'nsert 'k'mer 'p'lacement 'd'epth (-pt dpkt)
  -dsw <i> : depth smoothing window, averaging over position (-dsw 10000)
  -psw <i> : placement smoothing window (-psw 1000)
  -ksw <i> : kmer smoothing window (-ksw 1000)
  -isw <i> : insert smoothing window (-ksw 1000)
  -tp <f>  : threshold percentage, see paper (-tp 0.01)
  -tw <f>  : threshold width, see paper (-tw 1000)
  -td <f>  : threshold depth, see paper (-td -5.0)
  -std <f> : subplot threshold depth (-std -30.0)
  -sl <i>  : subplot length (-sl 5000)
  -plt <f> : plot threshold, only plot if more than % of errors (-plt 0.0)
  -fn <s>  : figure name (default: contig name)
  -mps <i> : minimum plot size in bp (-mps 100)
  -sc <s>  : plot only a specific contig (ie -sc contigName213)
  -pmo     : plot meta information only (off)
  -dpm     : don't plot meta information at all (off)
  -mgm <i> : maximum gaussian mixtures (5)
"""
__author__ = "Scott Clark <sc932 at cornell dot edu>"
__copyright__ = """
                                   ALE

                             COPYRIGHT NOTICE
                               (insert GPL)

Scott Clark Copyright 2011
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy
import sys
import heapq
import mixture # http://www.pymix.org/pymix/
import Smooth
import ProgressBar
import CmdIn

READ_LEN = 36
INSERT_LEN = 200

#import logging

class ALEFigure(object):
    """an ALE figure"""

    def __init__(self, start, end):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.ax2 = self.ax.twinx()

        self.start = start
        self.end = end
        
        self.threshold_window_set = []

    def add_threshold_windows(self, windows):
        self.threshold_window_set.extend(windows)

    def threshold_windows_in_last_half(self):
        for thresh_window in self.threshold_window_set:
            if thresh_window.start < (self.end-self.start)/2:
                return False
        return True

    def threshold_windows_in_first_half(self):
        for thresh_window in self.threshold_window_set:
            if thresh_window.start > (self.end-self.start)/2:
                return False
        return True

    def plot_threshold_windows(self):        
        # plot the thresholds
        for thresh_window in self.threshold_window_set:
            print "thresh window [%d,%d] %s" % (thresh_window.start, thresh_window.end, thresh_window.type_of)
            self.ax.axvspan(thresh_window.start, thresh_window.end, facecolor='r', alpha=0.1)          

    def get_percent_thresholded(self):
        total_below_threshold = numpy.zeros(self.end - self.start)
        for thresh_window in self.threshold_window_set:
            total_below_threshold[thresh_window.start:thresh_window.end] += 1
        summer = 0
        for pos in total_below_threshold:
            if pos:
                summer += 1
        return float(summer)/float(self.end - self.start)

    def set_labels(self, ax, plot_type, prob_vecs, std_witdh=2, twin=False, std_width=2):
        """Sets labels on y-axis for each subplot"""
        number_marks=3
        ticks = []
        labels = []
        i = -1
        for typer in plot_type:
            i += 1
            for j in range(-number_marks, number_marks + 1):
                ticks.append(4 + 7*i + j)
                if not twin:
                    if j < 0:
                        labels.append(str(j*std_width) + '$\sigma$')
                        #labels.append(str(j) + '$\sigma$ = ' + str(data_dict[typer] - j*sigma)[0:5])
                        #labels.append(str(data_dict[typer] - j*sigma)[0:5])
                    else:
                        labels.append(' ') # leave positive labels blank
                        #labels.append('+' + str(j) + '$\sigma$')
                        #labels.append('+' + str(j) + '$\sigma$ = ' + str(data_dict[typer] + j*sigma)[0:5])
                        #labels.append(str(data_dict[typer] + j*sigma)[0:5])
                else:
                    # TODO clean up
                    if j == 0:
                        labels.append(str(prob_vecs[typer].thresh_main_mean)[:7])
                    elif j < 0:
                        labels.append('(' + str(-j*std_width*prob_vecs[typer].thresh_main_std)[:6] + ')')
                    else:
                        labels.append(' ')

        i = 0
        if not twin:
            special_labels = {'i':'Insert', 'd':'Depth', 'p':'Place', 'k':'K-mer'}
            for typer in plot_type:
                labels[3 + 7*i] = special_labels[typer] + ' ' + '$\mu$' #str(data_dict[typer])[0:5]
                i += 1

        ax.set_yticks(ticks)
        ax.set_yticklabels(labels)

        labels = ax.get_yticklabels()
        for i, typer in enumerate(list(plot_type)):
            for j in range(2*number_marks + 1):
                labels[i*number_marks*2 + i + j].set_color(prob_vecs[typer].color)

    def plot_std_marks(self, current_subplot, length, color, number=3):
        """Plots std marks every sigma from the mean for each subplot"""
        for i in range(-number, number + 1):
            self.ax.plot([0, length], [4 + 7*current_subplot + i, 4 + 7*current_subplot + i], color + '--', alpha = 0.1)

    def format_data_for_plot(self, current_subplot, mean, sigma, data, start, end, std_width=2.0):
        """Formats data to be plotted on the predefined grid"""
        # scale
        data = 1.0/(sigma*std_width)*data[start:end]

        #mu = numpy.mean(data)

        # translate
        # want to put a line every 7 (3 sig buffer), with a 1 sig buffer on top/bottom
        data = (4 + 7*current_subplot) - 1.0/(sigma*std_width)*mean + data

        for point in data:
            if point < 0.0:
                point = 0.0

        return data

    def plot_cleanup(self, start, end, num_subplots):
        self.ax.set_ylabel('Avg (log) Likelihood')
        self.ax.set_xlim((0, end-start))
        self.ax.set_ylim((0,num_subplots*7))

class PositionViolation(object):
    """A threshold of an ALE score vector"""
    def __init__(self, type_of, pos, score, contig_name):
        self.contig_name = contig_name
        self.type_of = type_of
        self.pos = pos
        self.score = score
        self.igv_variants = []

    def __str__(self):
        if self.igv_variants:
            return_str = "%s:%d-%d\t%s\t%lf" % (self.contig_name, self.pos+1, self.pos+1, self.type_of, self.score)
            for var in self.igv_variants:
                return_str += " *** %s" % str(var)
        else:
            return "%s:%d-%d\t%s\t%lf" % (self.contig_name, self.pos+1, self.pos+1, self.type_of, self.score)

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

    def apply_variants(self, variants):
        for var in variants:
            assert(var.contig_name == self.contig_name)
            if self.contains_error(var):
                self.igv_variants.append(var)
                var.called_exactly = True
                self.within_read_of_error = True
                var.called_within_read = True
                self.within_insert_of_error = True
                var.called_within_insert = True
            if self.near_error(READ_LEN, var):
                var.called_within_read = True
                if not self.within_read_of_error:
                    self.within_read_of_error = True
            if self.near_error(INSERT_LEN, var):
                var.called_within_insert = True
                if not self.within_insert_of_error:
                    self.within_insert_of_error = True

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

    def smooth_prob(self, smoothing_width=None):
        """Smooth a specific prob vector"""
        if smoothing_width and smoothing_width < 0:
            raise ValueError("smoothing_width must be >= 0")

        if not smoothing_width:
            self.smoothing_width = int(numpy.max((1, self.length/250.0)))
        else:
            self.smoothing_width = smoothing_width

        self.prob_smoothed = Smooth.smooth(self.prob, self.smoothing_width)[self.smoothing_width:-self.smoothing_width]

    def get_threshold_windows(self, data, user_params, typer='?', thresh_mult=-5.0):
        """Returns the start and end points of the windows that cross the threshold with some constraints"""
        # TODO make more pythonic
        # get the starts
        total_sigma = self.thresh_main_std
        data_mean = self.thresh_main_mean
        # TODO base off of user_params
        threshold = thresh_mult*total_sigma
        cross_thresh = user_params.get("threshold_percent")
        len_thresh = user_params.get("threshold_width")

        line_threshold = threshold/total_sigma + data_mean/total_sigma
        end_point = end_all = len(data)
        def get_thresh(data, limits=[]):
            end_started = False
            ends = []
            end_started = False

            if limits:
                limit_on_num = len(limits) - 1
                limit_on = limits[-1]
            else:
                limit_on_num = 0
                limit_on = 0
            
            for position, point in enumerate(data):
                if end_all - position <= limit_on and end_started:
                    end_started = False
                    ends.append(end_point)
                    if limits:
                        limit_on_num -= 1
                        if limit_on_num < 0:
                            break
                        limit_on = limits[limit_on_num]
                else:
                    if point/total_sigma < line_threshold:
                        if not end_started:
                            end_started = True
                            end_point = position
                            end_window_len = 1
                            end_cross_total = 1
                        else:
                            end_window_len += 1
                            end_cross_total += 1
                    else:
                        if not end_started:
                            pass
                        else:
                            end_window_len += 1
                            if float(end_cross_total)/float(end_window_len) < cross_thresh:
                                end_started = False
                                if end_window_len - 1 > len_thresh:
                                    ends.append(end_point)
                                    if limits:
                                        limit_on_num -= 1
                                        if limit_on_num < 0:
                                            break
                                        limit_on = limits[limit_on_num]

            if end_started and limit_on_num >= 0:
                ends.append(end_point)
            return ends

        starts = get_thresh(data, limits=[])
        print starts
        reverse_data = data[::-1] # reverse
        ends = get_thresh(reverse_data, limits=starts)
        ends.reverse()
        ends = end_all - numpy.array(ends)

        print "starts and ends of windowing threshold"
        print starts, ends

        threshold_windows = []

        for i in range(len(starts)):
            #score = numpy.sum(data[starts[i]:ends[i]])/float(ends[i]-starts[i])
            score = -1
            threshold_windows.append(ThresholdViolation(typer, starts[i], ends[i], score, self.contig_name))

        return threshold_windows

    def find_threshold(self, user_params):
        """Finds the thresholds for errors given the data using Gaussian Mixture Model

        Args:
            data: The data to fit

        Kwargs:
            method: Whether to us [min,median,mean] of data in each bin
            thresh: Threshold for find_alpha
            bins: Number of pieces of the data we look at
            plot: Whether to plot the cdf and the two alpha cutoffs

        Returns:
            A soft threshold (alpha0) and A strong threshold (alpha1)

        Raises:
            
        """

        max_gauss_mixtures = user_params.get("max_gauss_mixtures")
        data = self.prob_smoothed

        #print data

        # http://www.pymix.org/pymix/index.php?n=PyMix.Tutorial

        # make two gaussains
        gaussian_one = mixture.NormalDistribution(numpy.mean(data),numpy.std(data))
        gaussian_two = mixture.NormalDistribution(10.0*numpy.mean(data),numpy.std(data))

        mixture_model = mixture.MixtureModel(2, [0.99,0.01], [gaussian_one, gaussian_two])

        # print mixture_model

        EM_tuned = False
        while not EM_tuned:
            # make mix_data from a random 10% of the original data
            index_array = numpy.arange(data.size)
            numpy.random.shuffle(index_array)
            mix_data = mixture.DataSet()
            data_size = numpy.min((int(numpy.floor(data.size/10.0)),50000))
            mix_data.fromArray(data[index_array[:data_size]])

            try:
                mixture_model.randMaxEM(mix_data, max_gauss_mixtures, 40, 0.001, silent=True)
                EM_tuned = True
            except AssertionError:
                # pymix likes to throw assertion errors when it has small machine precision errors...
                print "Caught an assertion error in pymix, randomizing input and trying again"
            except:
                print "pymix failed to find mixture model, using single gaussian"
                gaussian_two = mixture.NormalDistribution(numpy.mean(data),numpy.std(data))
                EM_tuned = True

        #print mixture_model

        # hacky, no good api access to the model components
        gauss_one_mean = float(str(mixture_model.components[0][0]).split('[')[1].split(',')[0])
        gauss_one_std = float(str(mixture_model.components[0][0]).split(', ')[1].split(']')[0])

        gauss_two_mean = float(str(mixture_model.components[1][0]).split('[')[1].split(',')[0])
        gauss_two_std = float(str(mixture_model.components[1][0]).split(', ')[1].split(']')[0])

        print "Gauss1: mu: %f, std: %f" % (gauss_one_mean, gauss_one_std)
        print "Gauss2: mu: %f, std: %f" % (gauss_two_mean, gauss_two_std)

        #print "Using threshold %f" % threshold

        # inv normal cdf
        if gauss_one_mean > gauss_two_mean or mixture_model.pi[1] < 0.60:
            self.thresh_main_mean = gauss_one_mean
            self.thresh_main_std = gauss_one_std
        else:
            self.thresh_main_mean = gauss_two_mean
            self.thresh_main_std = gauss_two_std

        # may be an error in pymix
        #if self.thresh_main_std == 0.1:
        #    self.thresh_main_std = numpy.std(data)

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

        self.main_figure = ALEFigure(0, length)

        self.pre_plot_run = False

    def find_N_worst(self, typer, N, user_param):
        """Finds the N worst scores in each of the likelihood vectors

        uses a set of min heaps to do this in O(length) time
        ignores edges if user_param.get("ignores_edges") is true
        an edge is defined by user_param.get("min_plot_size")/2
        
        returns an array of PositionViolation(s)
        """
        ignore_edges = user_param.get("ignore_edges")
        edge_length = user_param.get("min_plot_size")/2
        ignore_depth_below = user_param.get("ignore_depth_below")
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
            tv.within_read_of_hard_stop = self.low_depth_within_len(position[1], position[2], READ_LEN, ignore_depth_below)
            tv.within_insert_of_hard_stop = self.low_depth_within_len(position[1], position[2], INSERT_LEN, ignore_depth_below)
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
                    if igv_variants:
                        vi.apply_variants(igv_variants[self.name])
                    
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


        return float(c_exact)/float(len(seeds)), float(c_read)/float(len(seeds)), float(c_insert)/float(len(seeds))

    def pre_plot(self, user_params):
        """smooth the prob for the main figure and find the thresholds"""
        if self.pre_plot_run:
            print "Only need to run pre_plot() once! Ignoring..."
            return 0
        else:
            self.pre_plot_run = True
        # load parameters
        start = user_params.get("start")
        self.start = start
        end = user_params.get("end")
        self.end = end
        if end == 0:
            end = self.length

        # SMOOTHING
        # if there are no smoothing widths set them to a fraction of the window
        smoothing_dict = {'d':user_params.get("depth_smoothing_width"), 'p':user_params.get("placement_smoothing_width"), 'i':user_params.get("insert_smoothing_width"), 'k':user_params.get("kmer_smoothing_width")}

        largest_smooth = 0
        for typer in self.prob_vecs:
            if typer in user_params.get("plot_type"):
                print "smoothing %s" % typer
                self.prob_vecs[typer].smooth_prob(smoothing_width=smoothing_dict[typer])
                if self.prob_vecs[typer].smoothing_width > largest_smooth:
                    largest_smooth = self.prob_vecs[typer].smoothing_width

        if start < largest_smooth:
            self.start = largest_smooth
        if end > self.length - largest_smooth:
            self.end = self.length - largest_smooth

        self.largest_smooth = largest_smooth

        # THRESHOLDING
        for typer in user_params.get("plot_type"):
            print "Thresholding for %s" % typer
            self.prob_vecs[typer].find_threshold(user_params)
            print "mean, std", self.prob_vecs[typer].thresh_main_mean, self.prob_vecs[typer].thresh_main_std

    def sub_plots(self, user_params, start, end, pdf_stream = None):
        # make a sub plot
        increment = user_params.get("subplot_length")
        exists_plot = False
        plotted_previous = False
        if end == 0:
            end = self.length
        for s in numpy.arange(start, end, increment):

            sub_plot_start = s
            sub_plot_end = numpy.min((s+increment, end))
            if exists_plot and not self.main_figure.threshold_windows_in_last_half():
                if self.main_figure.threshold_windows_in_first_half() and plotted_previous:
                    plotted_previous = False
                else:
                    self.save_plot(user_params, pdf_stream=pdf_stream)
                    plotted_previous = True
            else:
                plotted_previous = False
            print "making figure %d-%d of %d" % (sub_plot_start, sub_plot_end, end)
            self.plot(user_params, sub_plot_start, sub_plot_end, is_subplot=True)
            exists_plot = True

            sub_plot_start = numpy.min((s+increment/2, end))
            sub_plot_end = numpy.min((s+3*increment/2, end))
            if exists_plot and not self.main_figure.threshold_windows_in_last_half():
                if self.main_figure.threshold_windows_in_first_half() and plotted_previous:
                    plotted_previous = False
                else:
                    self.save_plot(user_params, pdf_stream=pdf_stream)
                    plotted_previous = True
            else:
                plotted_previous = False
            print "making figure %d-%d of %d" % (sub_plot_start, sub_plot_end, end)
            self.plot(user_params, numpy.min((s+increment/2, end)), numpy.min((s+3*increment/2, end)), is_subplot=True)
            exists_plot = True

    def save_plot(self, user_params, pdf_stream=None):
        save_figure = user_params.get("save_figure")
        plot_threshold = user_params.get("plot_threshold")
        percent_thresholded = self.main_figure.get_percent_thresholded()

        if percent_thresholded > plot_threshold:
            if save_figure:
                pdf_stream.savefig()
            else:
                plt.show()
    
    def plot(self, user_params, start, end, is_subplot=False):
        if not self.pre_plot_run:
            print "Running pre_plot..."
            self.pre_plot(user_params)

        # load parameters
        if end == 0:
            end = self.length

        # make a new figure
        self.main_figure = ALEFigure(start, end)

        plot_type = user_params.get("plot_type")
        min_plot_size = user_params.get("min_plot_size")

        if end-start < min_plot_size:
            print "lower than min plot size"
            return 0, []

        if is_subplot:
            thresh_mult = user_params.get("sub_threshold_depth")
            std_width = int(-2.0*thresh_mult/5.0)
        else:
            thresh_mult = user_params.get("threshold_depth")
            std_width = int(-2.0*thresh_mult/5.0)
        
        # Main plotting code

        # THRESHOLD WINDOWS
        plot_type = list(plot_type)
        plot_type.sort() # make sure t comes last, if at all
        for typer in plot_type:
            if not is_subplot:
                threshold_windows = self.prob_vecs[typer].get_threshold_windows(self.prob_vecs[typer].prob_smoothed[start:end], user_params, typer=typer, thresh_mult=thresh_mult)
            else:
                threshold_windows = self.prob_vecs[typer].get_threshold_windows(self.prob_vecs[typer].prob[start:end], user_params, typer=typer, thresh_mult=thresh_mult)
            self.main_figure.add_threshold_windows(threshold_windows)
        self.main_figure.plot_threshold_windows()
        percent_thresholded = self.main_figure.get_percent_thresholded()

        # DATA LINES
        current_subplot = 0
        for typer in plot_type:
            if not is_subplot:
                self.main_figure.ax.plot(self.main_figure.format_data_for_plot(current_subplot, self.prob_vecs[typer].thresh_main_mean, self.prob_vecs[typer].thresh_main_std, self.prob_vecs[typer].prob_smoothed, start, end), self.prob_vecs[typer].color)
            else:
                self.main_figure.ax.plot(self.main_figure.format_data_for_plot(current_subplot, self.prob_vecs[typer].thresh_main_mean, self.prob_vecs[typer].thresh_main_std, self.prob_vecs[typer].prob, start, end, std_width=std_width), self.prob_vecs[typer].color)
            self.main_figure.ax.plot([0, end - start], [4 + 7*current_subplot - 2.5,4 + 7*current_subplot - 2.5], 'black')
            self.main_figure.plot_std_marks(current_subplot, end - start, self.prob_vecs[typer].color)
            current_subplot += 1

        # LABELS, TITLES, ETC
        self.main_figure.set_labels(self.main_figure.ax2, plot_type, self.prob_vecs, twin=True, std_width=std_width)
        self.main_figure.set_labels(self.main_figure.ax, plot_type, self.prob_vecs, std_width=std_width)

        self.main_figure.ax.set_title('Average Likelihoods (' + str(percent_thresholded)[:6] + ' below threshold)')
        self.main_figure.ax.set_xlabel('Position (bp) (+' + str(numpy.max([start, self.largest_smooth])) + ')')
        self.main_figure.plot_cleanup(start, end, current_subplot)

        return percent_thresholded, self.main_figure.threshold_window_set

def plot_histogram(input_data, save_figure=False, pdf_stream=None):
    max_val = numpy.max(input_data)
    min_val = numpy.min(input_data)
    bin_size = (max_val - min_val)/100.0
    if bin_size < 1e-6: bin_size = 0.01
    histogram = numpy.zeros(100)
    for value in input_data:
        histogram[min((max((0,int(numpy.floor((value - min_val)/bin_size)))),99))] += 1

    fig = plt.figure()
    ax = fig.add_subplot(111)

    ax.hist(histogram, 100, normed=1, facecolor='green', alpha=0.75)

    ax.set_xlabel('threshold')
    ax.set_ylabel('distribution')

    if save_figure:
        pdf_stream.savefig()
    else:
        plt.show()

def read_in_info(placement_file):
    """Reads in an ALE placement file, returns a list of Contigs.

    Args:
        placement_file: An ALE placement file (*.ale)
            must be in the following format::

                # comments/metadata
                # can have multiple lines, all starting with #
                # Reference: gi|170079663|ref|NC_010473.1| 350000 24.3
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

    MINIMUM_VALUE = -60.0 # minimum value we allow a position to have

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
                contigs.append(Contig(tLen, name = tName))
                place = 0
                print "Reading in contig: " + tName + " len " + str(tLen)
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
    print "\nYou now have a list of contigs, try contig[0].plot()"
    return contigs

def get_variant_stats(contig_name, igv_variants):
    contig_variants = igv_variants[contig_name]
    num_var = len(contig_variants)
    num_hit = 0
    num_win_read = 0
    num_win_insert = 0
    for var in contig_variants:
        num_hit += int(var.called_exactly)
        num_win_read += int(var.called_within_read)
        num_win_insert += int(var.called_within_insert)
    print "%s: total: %d hit: %d win_read: %d win_insert: %d" % (contig_name, num_var, num_hit, num_win_read, num_win_insert)
    return float(num_hit)/float(num_var), float(num_win_read)/float(num_var), float(num_win_insert)/float(num_var)

def plot_snp_rate(name,ge,gr,gi,ve,vr,vi):
    fig = plt.figure()
    ax = fig.add_subplot(111)

    print ge, gr, gi
    print ve, vr, vi

    ax.plot(ge)
    ax.plot(gr)
    ax.plot(gi)
    ax.plot(ve)
    ax.plot(vr)
    ax.plot(vi)

    plt.show()

    #pdf_stream.savefig()


def main():
    # default parameter values
    user_params = CmdIn.CommandLineParameters()
    user_params.add_parameter("start", "-s", int, 0)
    user_params.add_parameter("end", "-e", int, 0)
    user_params.add_parameter("plot_type", "-pt", str, "idpk")
    user_params.add_parameter("save_figure", "-nosave", None, True)
    user_params.add_parameter("subplots_on", "-spo", None, False)
    user_params.add_parameter("depth_smoothing_width", "-dsw", int, 0)
    user_params.add_parameter("placement_smoothing_width", "-psw", int, 5)
    user_params.add_parameter("insert_smoothing_width", "-isw", int, 0)
    user_params.add_parameter("kmer_smoothing_width", "-ksw", int, 0)
    user_params.add_parameter("threshold_depth", "-td", float, -5.0)
    user_params.add_parameter("threshold_percent", "-tp", float, 0.01)
    user_params.add_parameter("threshold_width", "-tw", int, 1000)
    user_params.add_parameter("sub_threshold_depth", "-std", float, -30.0)
    user_params.add_parameter("subplot_length", "-sl", int, 5000)
    user_params.add_parameter("plot_threshold", "-plt", float, 0.0)
    user_params.add_parameter("figure_name", "-fn", str, "")
    user_params.add_parameter("plot_meta", "-dpm", None, True)
    user_params.add_parameter("plot_meta_only", "-pmo", None, False)
    user_params.add_parameter("specific_contig", "-sc", str, None)
    user_params.add_parameter("min_plot_size", "-mps", int, 100)
    user_params.add_parameter("max_gauss_mixtures", "-mgm", int, 5)
    user_params.add_parameter("N_worst_positions", "-nwp", int, 50)
    user_params.add_parameter("ignore_edges", "-ie", None, False)
    user_params.add_parameter("ignore_depth_below", "-idb", int, 2)
    user_params.add_input("ale_file")

    user_params.__full_usage__ = __full_usage__
    user_params.__usage__ = __usage__

    # read in command line arguments
    user_params.read_sys_args()

    # read in contigs
    contigs = read_in_info(user_params.get_input("ale_file"))

    save_figure = user_params.get("save_figure")
    # open up a pdf_stream and file for output
    if save_figure:
        figure_name = user_params.get("figure_name")
        if figure_name == "":
            figure_name = sys.argv[-1] + '.pdf'
        pdf_stream = PdfPages(figure_name)
    else:
        pdf_stream = None

    print "Generating figures..."

    meta_percent = []

    fout = open(sys.argv[-1] + ".thresholds", 'w')

    # generate the data, given the user params
    if not user_params.get("plot_meta_only"):
        for contig in contigs:
            if contig.length >= user_params.get("min_plot_size"):
                if not user_params.get("specific_contig") or user_params.get("specific_contig") == contig.name:
                    contig.pre_plot(user_params)
                    percent_thresholded, threshold_windows = contig.plot(user_params, user_params.get("start"), user_params.get("end"))
                    contig.save_plot(user_params, pdf_stream=pdf_stream)
                    if user_params.get("subplots_on"):
                        contig.sub_plots(user_params, user_params.get("start"), user_params.get("end"), pdf_stream=pdf_stream)
                    print "%s had (%f) thresholded." % (contig.name, percent_thresholded)
                    meta_percent.append(percent_thresholded)
                    for window in threshold_windows:
                        fout.write("%s:%i-%i\t%s\n" % (contig.name, window.start, window.end, window.type_of))
                    
    fout.close()

    
    # plot the meta data
    if len(meta_percent) > 1:
        if user_params.get("plot_meta") or user_params.get("plot_meta_only"):
            plot_histogram(meta_percent, save_figure=save_figure, pdf_stream=pdf_stream)

    # save and close the output
    if save_figure:
        print "saved file %s" % figure_name
        pdf_stream.close()

    print "Executed Sucessfully"

if __name__ == '__main__':
    main()
