#!/usr/bin/python

"""
/*
 * Copyright (C) 2012 Scott Clark. All rights reserved.
 *
 * Developed by:
 * Scott Clark
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
 *   3. The names of its contributors may NOT be used to endorse or promote 
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

# (C) 2012 Scott Clark

#import numpy # for sci comp
#import matplotlib.pylab as plt # for plotting
#import time # for timing
#import commands # print commands.getoutput(script)
import sys

class CommandLineParameters(object):
    """User params

    EXAMPLE:
    
    # default parameter values
    user_params = CommandLineParameters()
    user_params.add_parameter("start", "-s", int, 0)
    user_params.add_input("ale_file")

    # read in command line arguments
    user_params.read_sys_args()

    # get using

    start = user_params.get("start")
    ale_file = user_params.get("ale_file")
    """
    def __init__(self):

        self.__full_usage__ = """full usage"""
        self.__usage__ = """usage"""

        # parameters
        self.parameter_list = []
        self.type_list = []
        self.lookup_by_name = {}
        self.lookup_by_command = {}

        # inputs
        self.input_list = []
        self.input_names = {}

    def read_sys_args(self):
        if len(sys.argv) < 2:
            print self.__usage__
            sys.exit(0)
            
        if sys.argv[1] == '--help' or sys.argv[1] == '-h' or sys.argv[1] == '-help' or sys.argv[1] == '--h':
            print self.__full_usage__
            sys.exit(0)

        if len(sys.argv) == 2:
            arg_on = 1
        else:
            # read in command line arguments
            arg_on = 1
            while(arg_on + 1 < len(sys.argv)):            
                if sys.argv[arg_on] in self.get_val_set():
                    # user defined value
                    self.set_by_command(sys.argv[arg_on], sys.argv[arg_on + 1])
                    arg_on += 2
                elif sys.argv[arg_on] in self.get_bool_set():
                    # flip the switch
                    self.set_by_command(sys.argv[arg_on], None)
                    arg_on += 1
                else:
                    print "Did not recognize command line argument %s." % sys.argv[arg_on]
                    print "Try -h for help."
                    exit(0)

        for i, ins in enumerate(sys.argv[arg_on:]):
            self.input_list[i] = ins        

    def add_input(self, name):
        number_on = len(self.input_list)
        self.input_list.append("")
        self.input_names[name] = number_on

    def get_input(self, name):
        if name not in self.input_names:
            raise ValueError("%s not in input_names" % name)
        return self.input_list[self.input_names[name]]

    def add_parameter(self, name=None, command=None, val_type=str, init_val=None):
        if not name:
            raise ValueError("add_parameter needs a name supplied")
        if not command:
            raise ValueError("add_parameter needs a command supplied")
        number_on = len(self.parameter_list)
        self.parameter_list.append(init_val)
        self.type_list.append(val_type)
        self.lookup_by_name[name] = number_on
        self.lookup_by_command[command] = number_on

    def get_val_set(self):
        val_set = []
        for command in self.lookup_by_command:
            if self.type_list[self.lookup_by_command[command]] in (str, int, float):
                val_set.append(command)
        return val_set

    def get_bool_set(self):
        bool_set = []
        for command in self.lookup_by_command:
            if self.type_list[self.lookup_by_command[command]] not in (str, int, float):
                bool_set.append(command)
        return bool_set

    def set_by_lookup_val(self, lookup_val, value):
        if value != None:
            try:
                value = self.type_list[lookup_val](value) # try to cast to type
            except:
                raise ValueError("Error setting parameter of %s with %s", (str(self.type_list[lookup_val]), str(value)))
            self.parameter_list[lookup_val] = value
        else:
            self.parameter_list[lookup_val] = not self.parameter_list[lookup_val]

    def set_by_command(self, command=None, value=None):
        if command not in self.lookup_by_command:
            raise ValueError("%s command not in param table")
        else:
            lookup_val = self.lookup_by_command[command]

        self.set_by_lookup_val(lookup_val, value)

    def get_by_command(self, command=None):
        if command not in self.lookup_by_command:
            raise ValueError("%s command not in param table")
        return self.parameter_list[self.lookup_by_command[command]]

    def set(self, name=None, value=None):
        if name not in self.lookup_by_name:
            raise ValueError("%s name not in param table")
        else:
            lookup_val = self.lookup_by_name[name]

        self.set_by_lookup_val(lookup_val, value)

    def get(self, name=None):
        if name not in self.lookup_by_name:
            raise ValueError("%s name not in param table")
        return self.parameter_list[self.lookup_by_name[name]]

def main():
    pass

if __name__ == '__main__':
    main()

