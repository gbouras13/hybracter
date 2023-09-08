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

import sys
import os
import subprocess
import commands
import time
import random

ALE = '../src/ALE'
synthreadgen = '../src/synthReadGen'
error_script = '../src/artificial_errors.py'
plotter = '../src/plotter3.py'
original_file = 'Ecoli_first350k'

def run_it_through(file_name, error_opts):
    """Perform the following (if needed):

            1. make a fasta file using artificial_errors
            2. make a bowtie db
            3. run bowtie on the db and fasta file
            4. run ALE on the map and fasta
            5. run plotter3 on the ALE file
       """
    # make fasta file
    if not os.path.exists("%s.fna" % file_name):
        fasta_script = "%s %s %s.fna" % (error_script, error_opts, original_file)
        print fasta_script
        print commands.getoutput(fasta_script)
    # run bowtie on it
    if not os.path.exists("%s.1.ebwt" % file_name):
        bowtie_build = "bowtie-build %s.fna %s" % (file_name, file_name)
        print bowtie_build
        print commands.getoutput(bowtie_build)
    if not os.path.exists("%s.map.sam" % file_name):
        bowtie_script = "bowtie -t -I 0 -X 300 --fr -S --threads 2 %s -1 part1_Ecoli_first350k.fastq  -2 part2_Ecoli_first350k.fastq %s.map.sam" % (file_name, file_name)
        print bowtie_script
        print commands.getoutput(bowtie_script)
    # run ALE on it
    if not os.path.exists("%s.ale" % file_name):
        ALE_script = "%s %s.map.sam %s.fna %s.ale" % (ALE, file_name, file_name, file_name)
        print ALE_script
        print commands.getoutput(ALE_script)
    # run the plotter
    if not os.path.exists("%s.ale.pdf" % file_name):
        plot_script = "%s %s.ale" % (plotter, file_name)
        print plot_script
        print commands.getoutput(plot_script)

def main():
    """Perform the following (if needed):
            
            1. Download an E.Coli genome
            2. Truncate it to 350k bases
            3. Synthesize 2M reads
            4. Run ALE on some basic transformations
            5. Generate output plots
    """
    os.chdir("../example")
    
    if not os.path.exists("CP000948.fna"):
        # wget from ncbi
        wget_command = "wget --user anonymous ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Bacteria/Escherichia_coli_K_12_substr__DH10B_uid20079/CP000948.fna"
        print wget_command
        print commands.getoutput(wget_command)
        
    if not os.path.exists("%s.fna" % original_file):
        # wget from ncbi and cat first 350k reads
        head_command = "head -5001 CP000948.fna > %s.fna.temp" % original_file
        print head_command
        print commands.getoutput(head_command)

        # clean it of ambiguity codes
        fin = open("%s.fna.temp" % original_file, "r")
        fout = open("%s.fna" % original_file, "w")
        line = fin.readline()
        fout.write(line)
        for line in fin:
            line_string = ""
            for base in line:
                if base not in ['A','T','C','G','\n']:
                    line_string += random.choice(['A','T','C','G'])
                else:
                    line_string += base
            fout.write(line_string)
        fin.close()
        fout.close()


    if not os.path.exists("part1_%s.fastq" % original_file) or not os.path.exists("part2_%s.fastq" % original_file):
        # make reads
        synth_command = "%s -ip 1.0 -nr 200000 -ps 10 -b %s.fna %s.fastq" % (synthreadgen, original_file, original_file)
        print synth_command
        print commands.getoutput(synth_command)

    ### 1 sub/indel errors
    file_name = '%s_sub_errors_1' % original_file
    error_opts = "-ase 75000 1 -ade 175000 1 -aie 275000 1 -o %s.fna" % file_name
    run_it_through(file_name, error_opts)
    
    ### 2 sub/indel errors
    file_name = '%s_sub_errors_2' % original_file
    error_opts = "-ase 75000 2 -ade 175000 2 -aie 275000 2 -o %s.fna" % file_name
    run_it_through(file_name, error_opts)

    ### 3 sub/indel errors
    file_name = '%s_sub_errors_3' % original_file
    error_opts = "-ase 75000 3 -ade 175000 3 -aie 275000 3 -o %s.fna" % file_name
    run_it_through(file_name, error_opts)
    
    ### transpose errors
    file_name = '%s_transpose_error' % original_file
    error_opts = "-trp 175000 175500 174500 -o %s.fna" % file_name
    run_it_through(file_name, error_opts)

    ### inversion errors
    file_name = '%s_inversion_error_10' % original_file
    error_opts = "-inv 175000 10 -o %s.fna" % file_name
    run_it_through(file_name, error_opts)

    ### inversion errors
    file_name = '%s_inversion_error_50' % original_file
    error_opts = "-inv 175000 50 -o %s.fna" % file_name
    run_it_through(file_name, error_opts)

    ### inversion errors
    file_name = '%s_inversion_error_100' % original_file
    error_opts = "-inv 175000 100 -o %s.fna" % file_name
    run_it_through(file_name, error_opts)

    ### inversion errors
    file_name = '%s_inversion_error_200' % original_file
    error_opts = "-inv 175000 200 -o %s.fna" % file_name
    run_it_through(file_name, error_opts)

    ### sub/indel errors
    #file_name = '%s_broken_apart' % original_file
    #error_opts = "-ab 50000 -ab 150000 -o %s.fna" % file_name
    #run_it_through(file_name, error_opts)

    ### copy error
    file_name = '%s_copy_error_77' % original_file
    error_opts = "-cip 150000 77 -o %s.fna" % file_name
    run_it_through(file_name, error_opts)

    ### copy error
    file_name = '%s_copy_error_200' % original_file
    error_opts = "-cip 150000 200 -o %s.fna" % file_name
    run_it_through(file_name, error_opts)

    ### copy error
    file_name = '%s_copy_error_400' % original_file
    error_opts = "-cip 150000 400 -o %s.fna" % file_name
    run_it_through(file_name, error_opts)

    ### copy error
    file_name = '%s_copy_error' % original_file
    error_opts = "-cip 150000 50000 -o %s.fna" % file_name
    run_it_through(file_name, error_opts)

    ### figure_2
    file_name = '%s_figure_two' % original_file
    error_opts = "-ase 69000 1 -ade 70000 1 -aie 71000 1 -inv 140000 10 -trp 210000 210200 209800 -cip 280000 77 -o %s.fna" % file_name
    run_it_through(file_name, error_opts)

if __name__ == '__main__':
    main()

