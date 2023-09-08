################################################################################
# 
#       This file is part of the Python Mixture Package
#
#       file:    mixture.py  
#       author: Benjamin Georgi 
#  
#       Copyright (C) 2004-2009 Benjamin Georgi
#       Copyright (C) 2004-2009 Max-Planck-Institut fuer Molekulare Genetik,
#                               Berlin
#
#       Contact: georgi@molgen.mpg.de
#
#       This library is free software; you can redistribute it and/or
#       modify it under the terms of the GNU Library General Public
#       License as published by the Free Software Foundation; either
#       version 2 of the License, or (at your option) any later version.
#
#       This library is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#       Library General Public License for more details.
#
#       You should have received a copy of the GNU Library General Public
#       License along with this library; if not, write to the Free
#       Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#
#
#
################################################################################

"""   
 
PyMix - Python Mixture Package

The PyMix library implements algorithms and data structures for data mining 
with finite mixture models. The framework is object oriented and 
organized in a hierarchical fashion.
 


"""
import random
import math
import copy
from alphabet import Alphabet,IntegerRange
from string import split,replace,join
import numpy
from numpy import linalg as la
import setPartitions
import sys
import re
import cStringIO, tokenize

import _C_mixextend  # import C extension


# ToDo: use logging package for verbosity control, makes all 'silent' parameters superfluous,
# needs to trickle down into the whole library 
import logging
log = logging.getLogger("PyMix")

# creating StreamHandler to stderr
hdlr = logging.StreamHandler(sys.stderr)

# setting message format
#fmt = logging.Formatter("%(name)s %(asctime)s %(filename)s:%(lineno)d %(levelname)s %(thread)-5s - %(message)s")
fmt = logging.Formatter("%(name)s %(filename)s:%(lineno)d - %(message)s")
hdlr.setFormatter(fmt)

# adding handler to logger object
log.addHandler(hdlr)

# Set the minimal severity of a message to be shown. The levels in
# increasing severity are: DEBUG, INFO, WARNING, ERROR, CRITICAL
log.setLevel(logging.ERROR)



# By default numpy produces a warning whenever we call  numpy.log with an array
# containing zero values. Usually this will happen a lot, so we change  numpys error handling
# to ignore this error. Since  numpy.log returns -inf for zero arguments the computations run
# through just fine.
numpy.seterr(divide="ignore",invalid="raise")

class DataSet:
    """ 
        Class DataSet is the central data object.
    """
    def __init__(self):
        """
        Creates and returns an empty DataSet object
        """
        self.N = None   # number of samples
        self.p = None   # number of dimensions 
        self.seq_p = None  # number of GHMM sequence features        
        self.complex = None
        
        self.sampleIDs = []   # unique sample ids, row label
        self.headers = []     # label for each column
        self.dataMatrix = []
        
        # attributes for internal data representation by sufficient statistics
        # these attributes are context specific (e.g. they depend on the MixtureModel)
        # and are initialised by the internalInit method.
        self.internalData = None
        
        self._internalData_views = None  # list of feature-wise views on internalData
        
        self.suff_dataRange = None
        self.suff_p = None
        self.suff_p_list = []


        # for each feature we can define a symbol (or value in case of continuous features)
        # to represent missing values in the data. This is used for instance in modelInitialization
        # to prevent the placeholders for missing data to influence the initial parameters.
        self.missingSymbols = {}

        
    def __len__(self):
        """
        Returns the number of samples in the DataSet.
        
        @return: Number of samples in the DataSet.
        """
        return self.N


    def __copy__(self):
        """
        Interface to copy.copy function.
        
        @return: deep copy of 'self'
        """
        
        cop = DataSet()
        
        cop.N = self.N    
        cop.p = self.p
        cop.sampleIDs = copy.copy(self.sampleIDs)
        cop.headers = copy.copy(self.headers)
        cop.dataMatrix = copy.deepcopy(self.dataMatrix)
        cop.internalData = copy.deepcopy(self.internalData)
        
        cop._internalData_views = copy.deepcopy(self._internalData_views) # XXX
        
        cop.suff_dataRange = copy.copy(self.suff_dataRange)
        cop.suff_p = self.suff_p
        cop.suff_p_list = copy.deepcopy(self.suff_p_list)
        return cop
    
    
    def fromArray(self,array, IDs = None, col_header = None):
        """
        Initializes the data set from a 'numpy' object.
        
        @param array: 'numpy' object containing the data
        @param IDs: sample IDs (optional)
        @param col_header: feature headers (optional)
        """
        
        self.complex = 0  # DataSet is not complex
        self.seq_p = 0
        
        self.N = len(array)
        try:
            self.p = len(array[0])
        except TypeError:  # if len() raises an exception array[0] is not a list -> p = 1
            self.p = 1    
        
        if not IDs:
            self.sampleIDs = range(self.N)
        else:
            self.sampleIDs = IDs

        if not col_header:
            self.headers = range(self.p)
        else:
            self.headers = col_header          

        self.dataMatrix = array.tolist()  
        

    def fromList(self,List, IDs = None, col_header = None):
        """
        Initializes the data set from a Python list.
        
        @param List: Python list containing the data
        @param IDs: sample IDs (optional)
        @param col_header: feature headers (optional)
        """
        self.complex = 0  # DataSet is not complex
        self.seq_p = 0
        
        self.N = len(List)
        try:
            self.p = len(List[0])
        except TypeError:  # if len() raises an exception array[0] is not a list -> p = 1
            self.p = 1    
        
        if IDs:
            self.sampleIDs = IDs
        else:
            self.sampleIDs = range(self.N)
        if col_header:
            self.headers = col_header   
        else:
            self.headers = range(self.p)            

        self.dataMatrix = List

        
    def fromFiles(self, fileNames, sep = "\t", missing="*",fileID = None, IDheader=False, IDindex=None):
        """
        Initializes the data set from a list of data flat files.
        
        @param fileNames: list of data flat files
        @param sep: separator string between values in flat files, tab is default
        @param missing: symbol for missing data '*' is default
        @param fileID: optional prefix for all features in the file
        @param IDheader: flag whether the sample ID column has a header in the first line of the flat files
        @param IDindex: index where the sample ids can be found, 0 by default
        """
        
        if IDindex == None:
            IDindex = [0] * len(fileNames)
        
        self.complex = 0  # DataSet is not complex
        self.seq_p = 0
        
        splitter = re.compile(sep)
        
        # store data in dictionary with sampleIDs as keys
        data_dict = {}
        data_nrs = [0] * len(fileNames)
        for q,fname in enumerate(fileNames):
            f = open(fname,"r")

            # parse header line
            l1 = f.next().rstrip()
            
            l1 = replace(l1,'"','')  # remove " characters, if present
            
            # split at separation characters
            #list1= split(l1,sep)   
            list1 = splitter.split(l1)   
            
            # prepending file identifier to column labels
            if fileID:
                print "File ident: ", fileID
                for i in range(len(list1)):
                    list1[i] = str(fileID)+"-"+str(list1[i])
            
            print fname,":",len(list1),"features"

            #print list1

            if IDheader == True:  # remove header for sample ID column, if present
                tt = list1.pop(IDindex[q])
                print tt
            #print list1            

            data_nrs[q] = len(list1)
            self.headers = self.headers + list1

            

            for h,line in enumerate(f):
                line = chomp(line)
                line = replace(line,'"','')  # remove " characters, if present 

                # cast data values into numerical types if possible
                sl =  splitter.split(line) 
                
                l = numerize(sl)
                sid = l.pop(IDindex[q])

                #print '->',sid
               
                if len(l) != data_nrs[q]:
                    print l
                    print list1
                    raise RuntimeError, "Different numbers of headers and data columns in files " + str(fname)+", sample "+str(sid) +" ,"+str(len(l))+" != "+str(data_nrs[q])
                
                if not data_dict.has_key(sid):
                    data_dict[sid] = {}
                    data_dict[sid][fname] = l

                else:
                    data_dict[sid][fname] = l

        # assembling data set from data dictionary
        for k in data_dict:
            self.sampleIDs.append(k)
            citem = []
            for q,fname in enumerate(fileNames):
                if data_dict[k].has_key(fname):
                    citem += data_dict[k][fname]
                else:
                    incomplete = 1
                    print "Warning: No data for sample "+str(k)+" in file "+str(fname)+"."
                    citem += [missing] * data_nrs[q]

            self.dataMatrix.append(citem)
           
        self.N = len(self.dataMatrix)
        # checking data-label consistency
        for i in range(self.N):
            assert len(self.dataMatrix[i]) == len(self.headers), "Different numbers of headers and data columns in files " + str(fileNames)+", sample "+str(self.sampleIDs[i]) +" ,"+str(len(self.dataMatrix[i])) +" != "+ str( len(self.headers))
            
        self.p = len(self.dataMatrix[0])


    def __str__(self):
        """
        String representation of the DataSet
        
        @return: string representation
        """
        strout = "Data set overview:\n"
        strout += "N = "+ str(self.N) + "\n"
        strout += "p = "+str(self.p) + "\n\n"       
        strout +=  "sampleIDs = " + str(self.sampleIDs)+ "\n\n"
        strout += "headers = "+ str(self.headers)+ "\n\n"
        #strout += "dataMatrix = "+str(self.dataMatrix) + "\n"

        return strout

    def printClustering(self,c,col_width=None):
        """
        Pretty print of a clustering .
        
        @param c: numpy array of integer cluster labels for each sample
        @param col_width: column width in spaces (optional)

        """
        if self.complex:
            raise NotImplementedError, "Needs to be done..."
        
        # get number of clusters in 'c'
        d = {}
        for lab in c:
            if lab == -1: # unassigned samples are handled seperately below
                continue
            d[lab]= ""
        G = len(d.keys())

        max_h = 0
        for h in self.headers:
            if len(str(h)) > max_h:
                max_h = len(str(h))
        max_sid = 0
        for s in self.sampleIDs:
            if len(str(s)) > max_sid:
                max_sid = len(str(s))
        if not col_width:
            space = max_h+2 
        else:
            space = col_width 
        
        for i in d:
            t = numpy.where(c == i)
            index = t[0]
            print "\n----------------------------------- cluster ",i,"------------------------------------"
            print ' ' * (max_sid+3),
            for k in range(len(self.headers)):
                hlen = len(str(self.headers[k]))
                print str(self.headers[k])+ " " * (space-hlen),
            print
            for j in range(len(index)):
                print '%-*s' % ( max_sid+3, self.sampleIDs[index[j]]),
                for k in range(len(self.dataMatrix[index[j]])):
                    dlen = len(str(self.dataMatrix[index[j]][k]))
                    print str(self.dataMatrix[index[j]][k] ) + " " * (space-dlen),
                print
        
        t = numpy.where(c == -1)
        index = t[0]
        if len(index) > 0:
            print "\n----------- Unassigned ----------------"
            space = max_h+2  
            print ' ' * (max_sid+3),
            for k in range(len(self.headers)):
                hlen = len(str(self.headers[k]))
                print self.headers[k]+ " " * (space-hlen),
            print
            for j in range(len(index)):
                print '%-*s' % ( max_sid+3, self.sampleIDs[index[j]]),
                for k in range(len(self.dataMatrix[index[j]])):
                    dlen = len(str(self.dataMatrix[index[j]][k]))
                    print str(self.dataMatrix[index[j]][k] ) + " " * (space-dlen),
                print

    def internalInit(self,m):
        """
        Initializes the internal representation of the data
        used by the EM algorithm .
        
        @param m: MixtureModel object
        """
        assert m.p == self.p,"Invalid dimensions in data and model."+str(m.p)+' '+str(self.p)
       
        templist = []
        for i in range(len(self.dataMatrix)):
            try:
                [t,dat] = m.components[0].formatData(self.dataMatrix[i])
            except InvalidDistributionInput, ex:
                ex.message += ' ( Sample '+str(self.sampleIDs[i])+', index = '+str(i)+' )'
                raise
            
            templist.append(dat)

        self.internalData = numpy.array(templist,dtype='Float64')
       
        if m.dist_nr > 1:
            self.suff_dataRange = copy.copy(m.components[0].suff_dataRange) 
        else:    
            self.suff_dataRange = [m.suff_p]
            
        self.suff_p = m.components[0].suff_p
        
        self._internalData_views = []
        for i in range(m.components[0].dist_nr):
            self.suff_p_list.append(m.components[0][i].suff_p)  # XXX suff_p_list should go away, only need in singleFeatureSubset
            if i == 0:
                prev_index = 0
            else:    
                prev_index =  self.suff_dataRange[i-1]  

            this_index =  self.suff_dataRange[i]
            if self.p == 1:   # only a single feature
                self._internalData_views.append(self.internalData)
            else:
                self._internalData_views.append( self.internalData[:,prev_index:this_index ] )

    def getInternalFeature(self, i):
        """
        Returns the columns of self.internalData containing the data of the feature with index 'i'
        
        @param i: feature index
        @return: numpy containing the data of feature 'i'
        """
        #assert self.suff_dataRange is not None,'DataSet needs to be initialized with .internalInit()'  
        if i < 0 or i >= len(self.suff_dataRange):
            raise IndexError, "Invalid index " + str(i)

        return self._internalData_views[i]

    def removeFeatures(self, ids,silent = 0):
        """
        Remove a list of features from the data set.
        
        @param ids: list of feature identifiers 
        @param silent: verbosity control
        """
        # removing columns from data matrix
        for i in ids:
            try:
                ind = self.headers.index(i)
            except ValueError:       
                sys.stderr.write("\nERROR:  Feature ID "+str(i)+ " not found.\n")
                raise

            for k in range(self.N):
                self.dataMatrix[k].pop(ind)
            
            r = self.headers.pop(ind)    
            self.p -= 1
            if not silent:
                print "Feature "+str(r)+" has been removed."

    def removeSamples(self, ids,silent = 0):
        """
        Remove a list of samples from the data set.
        
        @param ids: list of sample identifiers 
        @param silent: verbosity control
        """
        if self.internalData:
            print "Warning: internalInit has to be rerun after removeSamples."
            self.internalData = None
            self.suff_dataRange = None
            self.suff_p = None
            self.suff_p_list = None

        for si in ids:
            sind = self.sampleIDs.index(si)
            self.dataMatrix.pop(sind)
            self.sampleIDs.pop(sind)

        if not silent:
            print 'Samples '+str(ids)+' removed'

        self.N = self.N - len(ids)


    def filterSamples(self,fid,min_value,max_value):
        """
        Removes all samples with values < 'min_value' or > 'max_value' in feature 'fid'.
        
        @param fid: feature ID in self.headers
        @param min_value: minimal required value
        @param max_value: maximal required value
        
        """
        if self.internalData:
            print "Warning: internalInit has to be rerun after removeSamples."
            self.internalData = None
            self.suff_dataRange = None
            self.suff_p = None
            self.suff_p_list = None

        ind = self.headers.index(fid)
        
        print "Removing samples with "+fid +" < "+str(min_value)+" or > "+str(max_value)+" ...",
        i = 0  # current index in dataMatrix
        c = 0  # number of samples already considered
        r = 0  # number of removed samples
        while c < self.N:
            if self.dataMatrix[i][ind] < min_value or self.dataMatrix[i][ind] > max_value:
                # remove sample
                self.dataMatrix.pop(i)
                self.sampleIDs.pop(i)
                c += 1
                r += 1
            else:
                i += 1
                c += 1
        
        print str(r)+" samples removed"
        self.N = self.N - r
        

    def maskDataSet(self, valueToMask, maskValue, silent=False):
        """
        Allows the masking of a value with another in the entire data matrix.
        
        @param valueToMask: value to be masked
        @param maskValue: value which is to be substituted
        @param silent: verbosity control (False is default)
        """
        count = 0
        for i in range(self.N):
            for j in range(self.p):
                if self.dataMatrix[i][j] == valueToMask:
                    self.dataMatrix[i][j] = maskValue
                    count += 1

        if not silent:
            print str(count), "values '"+str(valueToMask)+"' masked with '"+str(maskValue)+"' in all features."


    def maskFeatures(self, headerList, valueToMask, maskValue):
        """
        Equivalent to maskDataSet but constrained to a subset of features
        
        @param headerList: list of features IDs
        @param valueToMask: value to be masked
        @param maskValue: value which is to be substituted
        """
        count = 0
        for h in headerList:
            try:
                ind = self.headers.index(h)
            except ValueError:       
                sys.stderr.write("\nERROR:  Feature ID "+str(h)+ " not found.\n")
                raise
            
            for j in range(self.N):
                if str(self.dataMatrix[j][ind]) == str(valueToMask):
                    self.dataMatrix[j][ind] = maskValue
                    count += 1

        print str(count), "values '"+str(valueToMask)+"' masked with '"+str(maskValue)+"' in "+str(len(headerList))+" features."



    def getExternalFeature(self, fid):
        """
        Returns the external data representation of a given feature

        @param fid: feature ID in self.headers

        @return: list of data samples for feature fid
        """
        index = self.headers.index(fid)
        res = []
        for i in range(self.N):
            res.append(self.dataMatrix[i][index])

        return res

    def extractSubset(self,ids):
        """
        Remove all samples in 'ids' from 'self' and return a new DataSet initialised with these samples
        
        @param ids: list of sample indices
        
        @return: DataSet object containing the samples in ids
        """
        res = DataSet()

        res.N = len(ids)
        res.p = self.p
        res.suff_dataRange = copy.copy(self.suff_dataRange)
        res.suff_p = self.suff_p
        res.suff_p_list = self.suff_p_list
        res.sampleIDs = ids
        res.headers = copy.copy(self.headers)        
        
        if self.internalData is not None:
            res.internalData = numpy.zeros((res.N,res.suff_p),self.internalData.type())

        else:
            res.internalData = None
        
        #remove subset entries from self.internalData 
        new_intData = None
        if self.internalData is not None:
            new_intData = numpy.zeros(( (self.N-res.N),res.suff_p),self.internalData.type())
            new_sampleIDs = []
            ni = 0
            for i in range(self.N):
                if self.sampleIDs[i] not in ids:
                    new_intData[ni] = self.internalData[i]
                    ni+=1

        for i,d in enumerate(ids):
            ind = self.sampleIDs.index(d)

            dat = self.dataMatrix.pop(ind)
            res.dataMatrix.append(dat)
    
            # fill internalData matrix
            if self.internalData is not None:
                res.internalData[i]  = copy.deepcopy(self.internalData[ind])

            self.sampleIDs.pop(ind)
            

        self.internalData = new_intData
        self.N -= res.N
        
        return res

    def singleFeatureSubset(self, index):
        """
        Returns a DataSet for the feature with internal index 'index' in 'self'.
        For internal use.
        
        @param index: feature index
        
        @return: DataSet object
        """
        res = DataSet()

        res.N = self.N    # number of samples
        res.p = 1    # number of dimensions 
        res.seq_p = self.seq_p  # number of GHMM sequence features        
        res.suff_p = self.suff_p_list[index]
        res.suff_p_list = [ self.suff_p_list[index] ]
        res.internalData = self.getInternalFeature(index)

        res._internalData_views = [self._internalData_views[index]] # XXX
        
        if self.headers:
            res.headers.append(self.headers[index])
        if self.sampleIDs:
            res.sampleIDs = self.sampleIDs 

        res.missingSymbols = {}
        if self.missingSymbols.has_key(index):
            res.missingSymbols[0] = self.missingSymbols[index]
        
        res.suff_dataRange = [self.suff_dataRange[index]-self.suff_dataRange[index-1]]
        return res

    def setMissingSymbols(self, findices, missing):
        """
        Assigns missing value placeholders to features.        
        
        @param findices: list of internal feature indices
        @param missing: list of missing symbols/values
        """
        assert len(findices) == len( missing)
        for i,h in enumerate(findices):
            self.missingSymbols[h] = missing[i]

    def getMissingIndices(self, ind):
        """
        Get indices of missing values in one feature
        
        @param ind: feature index
        
        @return: list of indices of missing values
        """
        assert self.suff_dataRange is not None
        
        if not self.missingSymbols.has_key(ind):
            return []
        else:
            m_ind = []
            dat_ind = self.getInternalFeature(ind)
            for i,v in enumerate(dat_ind):
                # check sample 'i' for missing symbol
                if  numpy.all(v == self.missingSymbols[ind]):
                    m_ind.append(i)
            return m_ind                


    def writeClusteringFasta(self,fn_pref,m):
        """
        Writes a clustering based on model 'm' into files in FASTA format.
        Note that this implies sequence data.

        @param fn_pref: Filename prefix. The full name of each output file consists 
        of the prefix, the cluster number and the extension .fa
        @param m: MixtureModel object 
        """
        c = m.classify(self,silent=1)

        # get number of clusters in 'c'
        d = {}
        for a in c:
            d[a]= ""
        G = len(d)
 
        # write data to file
        for k in d:
            t = numpy.where(c == k)
            index = t[0]
            f = open(fn_pref+'_'+str(k)+'.fa','w')
            for i in index:
                f.write(">"+str(self.sampleIDs[i])+"\n")
                s = join(self.dataMatrix[i],'') # convert to sequence
                f.write(s+"\n")
            f.close()
                    


def numerize(data):
    """
    Cast all elements in a list to numeric values.
    
    @param data: list of data
    
    @return: list of processed data
    """
    
    for i in range(len(data)):
        try:
            data[i] = int(data[i])

        except ValueError:
           try:
                data[i] = float(data[i]) 
           except ValueError:
                pass
    return data

def remove_col(matrix, index):
    """
    Removes a column in a Python matrix (list of lists)

    @param matrix: Python list of lists
    @param index: index of column to be removed
    
    @return: matrix with column deleted
    """
    for i in range(len(matrix)):
        del(matrix[i][index])
    return matrix



#-------------------------------------------------------------------------------------------


class MixtureError(Exception):
    """Base class for mixture exceptions."""
    def __init__(self, message):
        self.message = message
    def __str__(self):
        return str(self.message)

class InvalidPosteriorDistribution(MixtureError):
    """
    Raised if an invalid posterior distribution occurs.
    """
    def __init__(self,message):
       self.message = message
    def __str__(self):
        return str(self.message)

class InvalidDistributionInput(MixtureError):
    """
    Raised if a DataSet is found to be incompatible with a given MixtureModel.
    """
    def __init__(self,message):
       self.message = message
    def __str__(self):
        return str(self.message)

class ConvergenceFailureEM(MixtureError):
    """
    Raised if a DataSet is found to be incompatible with a given MixtureModel.
    """
    def __init__(self,message):
       self.message = message
    def __str__(self):
        return str(self.message)




#--------------------------------------------------------------------------------------------

class ProbDistribution:
    """
    Base class for all probability distributions.
    """
    def __init__(self):
        """
        Constructor
        """
        pass

    def __eq__(self,other):
        """
        Interface for the '==' operation
        
        @param other: object to be compared
        """
        raise NotImplementedError

    def __str__(self):
        """
        String representation of the DataSet
        
        @return: string representation
        """
        raise NotImplementedError

    def __copy__(self):
        "Interface for the copy.copy function"
        raise NotImplementedError
    
    def pdf(self,data):
        """
        Density function. 
        MUST accept either numpy or DataSet object of appropriate values. We use numpys as input
        for the atomar distributions for efficiency reasons (The cleaner solution would be to construct
        DataSet subset objects for the different features and we might switch over to doing that eventually).
        
        @param data: DataSet object or numpy array 
        
        @return: log-value of the density function for each sample in 'data'
        """
        raise NotImplementedError
        

    def MStep(self,posterior,data,mix_pi=None):       
        """
        Maximization step of the EM procedure. Reestimates the distribution parameters
        using the posterior distribution and the data.

        MUST accept either numpy or DataSet object of appropriate values. numpys are used as input
        for the atomar distributions for efficiency reasons 
        
        @param posterior: posterior distribution of component membership
        @param data: DataSet object or 'numpy' of samples
        @param mix_pi: mixture weights, necessary for MixtureModels as components.
        
        """
        raise NotImplementedError

        
    def sample(self):
        """
        Samples a single value from the distribution.
        
        @return: sampled value
        """
        raise NotImplementedError, "Needs implementation"
        
    def sampleSet(self,nr):
        """
        Samples several values from the distribution.
        
        @param nr: number of values to be sampled.

        @return: sampled values
        """
        raise NotImplementedError, "Needs implementation"

    def sufficientStatistics(self, posterior, data):
        """
        Returns sufficient statistics for a given data set and posterior. 
        
        @param posterior: numpy vector of component membership posteriors
        @param data: numpy vector holding the data 
        
        @return: list with dot(posterior, data) and dot(posterior, data**2)
        """
        raise NotImplementedError, "Needs implementation"


    def isValid(self,x):
        """
        Checks whether 'x' is a valid argument for the distribution and raises InvalidDistributionInput
        exception if that is not the case.
        
        @param x: single sample in external representation, i.e.. an entry of DataSet.dataMatrix
        
        @return: True/False flag
        """
        raise NotImplementedError
        
    def formatData(self,x):
        """
        Formats samples 'x' for inclusion into DataSet object. Used by DataSet.internalInit()
        
        @param x: list of samples
        
        @return: two element list: first element = dimension of self, second element = sufficient statistics for samples 'x'
        """
        return [self.p,x]
    
    
    def flatStr(self,offset):
        """
        Returns the model parameters as a string compatible
        with the WriteMixture/ReadMixture flat file 
        format.
        
        @param offset: number of '\t' characters to be used in the flatfile.
        """
        raise NotImplementedError
                
    def posteriorTraceback(self,x):
        """
        Returns the decoupled posterior distribution for each
        sample in 'x'. Used for analysis of clustering results.
        
        @param x: list of samples
        
        @return: decoupled posterior
        """
        raise NotImplementedError

    def update_suff_p(self):
        """
        Updates the .suff_p field.
        """
        return self.suff_p

    def merge(self,dlist,weights):
        """
        Merges 'self' with the distributions in'dlist' by an 
        convex combination of the parameters as determined by 'weights'
        
        @param dlist: list of distribution objects of the same type as 'self'
        @param weights: list of weights, need not to sum up to one
        """
        raise NotImplementedError


class PriorDistribution(ProbDistribution):
    """
    Prior distribution base class for the Bayesian framework
    """
    def pdf(self, m):    
        """
        Returns the log-density of the ProbDistribution object(s) 'm' under the
        prior.
        
        @param m: single appropriate ProbDistribution object or list of ProbDistribution objects
        """
        raise NotImplementedError, "Needs implementation"

    def posterior(self,m,x):    
        raise NotImplementedError, "Needs implementation"

    def marginal(self,x):
        raise NotImplementedError, "Needs implementation"

    def mapMStep(self, dist, posterior, data, mix_pi=None, dist_ind = None):       
        """
        Maximization step of the maximum aposteriori EM procedure. Reestimates the distribution parameters
        of argument 'dist' using the posterior distribution, the data and a conjugate prior.

        MUST accept either numpy or DataSet object of appropriate values. numpys are used as input
        for the atomar distributions for efficiency reasons.
        
        @param dist: distribution whose parameters are to be maximized
        @param posterior: posterior distribution of component membership
        @param data: DataSet object or 'numpy' of samples
        @param mix_pi: mixture weights, necessary for MixtureModels as components.
        @param dist_ind: optional index of 'dist', necessary for ConditionalGaussDistribution.mapMStep (XXX)
        """
        raise NotImplementedError


    def mapMStepMerge(self, group_list):
        """
        Computes the MAP parameter estimates for a candidate merge in the structure
        learning based on the information of two CandidateGroup objects.
        
        @param group_list: list of CandidateGroup objects
        @return: CandidateGroup object with MAP parameters
        """
        raise NotImplementedError

    def mapMStepSplit(self, toSplitFrom, toBeSplit):
        """
        Computes the MAP parameter estimates for a candidate merge in the structure
        learning based on the information of two CandidateGroup objects.
        
        @return: CandidateGroup object with MAP parameters

        """
        raise NotImplementedError

    def updateHyperparameters(self, dists, posterior, data): 
        """
        Update the hyperparameters in an empirical Bayes fashion.
        
        @param dists: list of ProbabilityDistribution objects
        @param posterior: numpy matrix of component membership posteriors
        @param data: DataSet object
        """
        raise NotImplementedError

#------------------------------------------------------------------------------------------------------

class NormalDistribution(ProbDistribution):
    """
    Univariate Normal Distribution    
    
    """
    
    def __init__(self, mu ,sigma):
        """
        Constructor

        @param mu: mean parameter
        @param sigma: standard deviation parameter
        """
        self.p = 1
        self.suff_p = 1
        self.mu = mu
        self.sigma = sigma
        
        self.freeParams = 2
        
        self.min_sigma = 0.1  # minimal standard deviation
    
    def __eq__(self,other):
        res = False
        if isinstance(other,NormalDistribution):
            if numpy.allclose(other.mu, self.mu) and numpy.allclose(other.sigma, self.sigma):
                res = True
        return res        

    def __copy__(self):
        return NormalDistribution(copy.deepcopy(self.mu), copy.deepcopy(self.sigma))
    
    def __str__(self):
        return "Normal:  ["+str(self.mu)+", "+str(self.sigma) + "]"

   
    def pdf(self, data):

        # Valid input arrays will have the form [[sample1],[sample2],...] or
        # [sample1,sample2, ...], the latter being the input format to the extension function,
        # so we might have to reformat the data
        if isinstance(data, DataSet ):  
            assert data.internalData is not None, "Internal data not initialized."
            nr = len(data.internalData)
            assert data.internalData.shape == (nr,1), 'shape = '+str(data.internalData.shape)
            
            x = numpy.transpose(data.internalData)[0]
            
        elif isinstance(data, numpy.ndarray): 
            nr = len(data)
            
            if data.shape == (nr,1):  # data format needs to be changed
                x = numpy.transpose(data)[0]
            elif data.shape == (nr,): 
                x = data
            else:
                raise TypeError,'Invalid data shape: '+str(data.shape)
        else:
            raise TypeError,"Unknown/Invalid input type:"+str(type(data))
        
        # computing log likelihood 
        res = _C_mixextend.wrap_gsl_ran_gaussian_pdf(self.mu, self.sigma, x)

        return numpy.log(res)

    def sample(self):
        return random.normalvariate(self.mu, self.sigma)

    
    def sampleSet(self,nr):
        res = numpy.zeros(nr,dtype='Float64')
        
        for i in range(nr):
            res[i] = self.sample()
            
        return res    
    
    def sufficientStatistics(self, posterior, data): 
        """
        Returns sufficient statistics for a given data set and posterior. In case of the Normal distribution
        this is the dot product of a vector of component membership posteriors with the data and the square
        of the data.
        
        @param posterior: numpy vector of component membership posteriors
        @param data: numpy vector holding the data 
        
        @return: list with dot(posterior, data) and dot(posterior, data**2)
        """
        return numpy.array( [numpy.dot(posterior, data)[0], numpy.dot(posterior, data**2)[0]], dtype='Float64')


    def MStep(self,posterior,data,mix_pi=None):           
        # data has to be reshaped for parameter estimation
        if isinstance(data,DataSet):
            x = data.internalData[:,0]
        elif isinstance(data,numpy.ndarray):
            x = data[:,0]
            
        else:
            raise TypeError, "Unknown/Invalid input to MStep."    
        nr = len(x)
        
        sh = x.shape
        assert sh == (nr,)  # XXX debug 
        
        post_sum = numpy.sum(posterior)

        # checking for valid posterior: if post_sum is zero, this component is invalid
        # for this data set
        if post_sum != 0.0:
            # computing ML estimates for mu and sigma
            new_mu =  numpy.dot(posterior, x) / post_sum
            new_sigma = math.sqrt(numpy.dot(posterior, (x - new_mu)**2 ) / post_sum)
        else:
            raise InvalidPosteriorDistribution, "Sum of posterior is zero: "+str(self)+" has zero likelihood for data set."
                        
        if new_sigma < self.min_sigma:
           # enforcing non zero variance estimate
            new_sigma = self.min_sigma 

        # assigning updated parameter values
        self.mu = new_mu
        self.sigma = new_sigma

    def isValid(self,x):
        try:
            float(x)
        except (ValueError):
            #print "Invalid data: ",x,"in NormalDistribution."
            raise InvalidDistributionInput, "\n\tInvalid data: "+str(x)+" in NormalDistribution."

    def formatData(self,x):
        if isinstance(x,list) and len(x) == 1: 
            x = x[0]
        self.isValid(x)  # make sure x is valid argument
        return [self.p,[x]]

        
    def flatStr(self,offset):
        offset +=1
        return "\t"*+offset + ";Norm;"+str(self.mu)+";"+str(self.sigma)+"\n"

    def posteriorTraceback(self,x):
        return self.pdf(x)

    def merge(self,dlist, weights):
        raise DeprecationWarning, 'Part of the outdated structure learning implementation.'
        assert len(dlist)+1 == len(weights)

        norm = sum(weights)
        m_mu = self.mu * weights[0]
        #print self.mu," * ", coeff ,"= ",m_mu
        
        for i in range(len(dlist)):
            assert isinstance(dlist[i],NormalDistribution)
            
            #print m_mu, " += ",
            
            m_mu += weights[i+1] * dlist[i].mu
            #print dlist[i].mu," * ", coeff ,"= ",m_mu
            #m_sigma += coeff * dlist[i].sigma
        
        m_mu = m_mu / norm
        m_sigma =  weights[0] *  ( self.sigma + ( self.mu -m_mu)**2 )
        for i in range(len(dlist)):    
            m_sigma +=  weights[i+1] * ( dlist[i].sigma + ( dlist[i].mu - m_mu)**2 )
 
        self.mu = m_mu
        self.sigma = m_sigma / norm



class MultiNormalDistribution(ProbDistribution):
    """
    Multivariate Normal Distribution    

    """
    def __init__(self,p, mu, sigma):
        """
        Constructor

        @param p: dimensionality of the distribution
        @param mu: mean parameter vector
        @param sigma: covariance matrix
        """

        assert len(mu) == len(sigma) == len(sigma[0]) == p, str(len(mu))+ ' == ' + str(len(sigma)) + ' == ' + str(len(sigma[0])) + ' == '+ str(p)
        self.p = p
        self.suff_p = p
        self.mu = numpy.array( mu, dtype='Float64')
        self.sigma = numpy.array(sigma,dtype='Float64')
        self.freeParams = p + p**2

 
    def __copy__(self):
        return MultiNormalDistribution(self.p,self.mu, self.sigma)

    
    def __str__(self):
        return "Normal:  ["+str(self.mu)+", "+str(self.sigma.tolist()) + "]"

    def __eq__(self,other):
        if not isinstance(other,MultiNormalDistribution):
            return False
        if self.p != other.p:
            return False
        if not numpy.allclose( self.mu, other.mu ) or not numpy.allclose( self.sigma, other.sigma):        
            return False
        return True

    def pdf(self, data):
        if isinstance(data, DataSet ):
            x = data.internalData
        elif isinstance(data, numpy.ndarray):
            x = data
        else:
            raise TypeError,"Unknown/Invalid input type." 
            
    	# initial part of the formula
    	# this code depends only on the model parameters ... optmize?
    	dd = la.det(self.sigma);
        inverse = la.inv( self.sigma);
        ff = math.pow(2*math.pi,-self.p/2.0)*math.pow(dd,-0.5);

    	# centered input values
        centered = numpy.subtract(x,numpy.repeat([self.mu],len(x), axis=0))

    	res = ff * numpy.exp(-0.5*numpy.sum(numpy.multiply(centered,numpy.dot(centered,inverse)),1))

        return numpy.log(res)

    def MStep(self,posterior,data,mix_pi=None):       

        if isinstance(data,DataSet):
            x = data.internalData
        elif isinstance(data,numpy.ndarray):
            x = data
        else:
            raise TypeError, "Unknown/Invalid input to MStep."    
        
        post = posterior.sum() # sum of posteriors
        self.mu = numpy.dot(posterior, x)/post

        # centered input values (with new mus)
        centered = numpy.subtract(x,numpy.repeat([self.mu],len(x), axis=0));
        self.sigma = numpy.dot(numpy.transpose(numpy.dot(numpy.identity(len(posterior))*posterior,centered)),centered)/post


    def sample(self, A = None):
        """
        Samples from the mulitvariate Normal distribution.
        
        @param A: optional Cholesky decomposition of the covariance matrix self.sigma, can speed up
        the sampling
        """
        if A == None:
            A = la.cholesky(self.sigma)
        
        z = numpy.zeros(self.p,dtype='Float64')
        for i in range(self.p):
            z[i] = random.normalvariate(0.0,1.0)  # sample p iid N(0,1) RVs

        X = numpy.dot(A,z) + self.mu
        return X.tolist()  # return value of sample must be Python list

    def sampleSet(self,nr):
        A = la.cholesky(self.sigma)
        res = numpy.zeros((nr,self.p),dtype='Float64')
        for i in range(nr):
            res[i,:] = self.sample(A=A)
        return res

    def isValid(self,x):
        if not len(x) == self.p:
            raise InvalidDistributionInput, "\n\tInvalid data: wrong dimension(s) "+str(len(x))+" in MultiNormalDistribution(p="+str(self.p)+")."
        for v in x:
            try:
                float(v)
            except (ValueError):
                raise InvalidDistributionInput, "\n\tInvalid data: "+str(x)+" in MultiNormalDistribution."
        
    def flatStr(self, offset):
        offset +=1
        return "\t"*offset+";MultiNormal;"+str(self.p)+";"+str(self.mu.tolist())+";"+str(self.sigma.tolist())+"\n"        



    
class ConditionalGaussDistribution(ProbDistribution):
    """
    Constructor for conditional Gauss distributions. Conditional Gaussians
    use a sparse formulation of the covariance matrix, which allows computationally
    efficient modeling of covariance for high-dimensional data.

    Also, the conditional Gaussians implement a tree dependency structure.

    """

    def __init__(self, p, mu, w, sigma, parents):
        """
        Constructor

        @param p: dimensionality of the distribution
        @param mu: mean parameter vector
        @param w: covariance weights (representing off-diagonal entries in the full covariance matrix)
        @param sigma: standard deviations (diagonal entries of the covariance matrix)
        @param parents: parents in the tree structure implied by w
        """
        assert p == len(mu) == len(w) == len(sigma) == len(parents)
        
        self.p = p
        self.suff_p = p
        self.freeParams = p * 3 
        self.mu = mu  # mean values
        self.w = w    # conditional weights
        self.sigma = sigma  # standard deviations
        
        self.parents = parents  # tree structure encoded by parent index relationship
        
    def __str__(self):
        return 'ConditionalGaussian: \nmu='+str(self.mu)+', \nsigma='+str(self.sigma)+', \nw='+str(self.w)+', \nparents='+str(self.parents)


    def sample(self):
        s = [None] * self.p
        s[0] = random.normalvariate(self.mu[0], self.sigma[0]) 
        
        for i in range(1,self.p):
            pid = self.parents[i]
            assert s[pid] != None   # XXX assumes that the features are in topological order
            shift_mu = self.mu[i] - (self.w[i]* self.mu[pid])
            s[i] = random.normalvariate( shift_mu + (self.w[i]* s[pid]), self.sigma[i]) 

        return s
        
    def sampleSet(self, nr):
        s = numpy.zeros((n,self.p))
        for i in range(nr):
            s[i,:] = self.sample()

        return s
        

    def pdf(self, data):
        
        # XXX assume root as first index
        assert self.parents[0] == -1
        assert self.w[0] == 0.0
        
        res = numpy.zeros(len(data))
        
        for i in range(len(data)):
            res[i] = math.log( (1.0 / (math.sqrt(2.0*math.pi) * self.sigma[0])) * math.exp(  ( data[i,0] - self.mu[0]  )**2 / (-2.0*self.sigma[0]**2) ))
            for j in range(1,self.p):
                pind = self.parents[j]
                res[i] += math.log( (1.0 / (math.sqrt(2.0*math.pi) * self.sigma[j])) * math.exp(  ( data[i,j] - self.mu[j] - self.w[j]* ( data[i,pind]-self.mu[pind] )  )**2 / (-2.0*self.sigma[j]**2) ))            

        return res


    def MStep(self,posterior,data,mix_pi=None):       
        var = {}
        post_sum = numpy.sum(posterior)

        # checking for valid posterior: if post_sum is zero, this component is invalid
        # for this data set
        if post_sum != 0.0:
            # reestimate mu
            for j in range(self.p):
                self.mu[j] =  numpy.dot(posterior, data[:,j]) / post_sum
                var[j] = numpy.dot(posterior, (data[:,j] - self.mu[j])**2 ) / post_sum
           
            for j in range(self.p):
                # computing ML estimates for w and sigma
                pid = self.parents[j]
                cov_j = numpy.dot(posterior, (data[:,j] - self.mu[j]) * (data[:,pid] - self.mu[pid])) / post_sum 
        
                if pid <> -1:  # has parents
                    self.w[j] = cov_j / var[pid]
                    print  var[j], self.w[j]**2, var[pid], var[j] - (self.w[j]**2 * var[pid])
                    self.sigma[j] = math.sqrt( var[j] - (self.w[j]**2 * var[pid]) )
                else:
                    self.sigma[j] = math.sqrt( var[j])
                
        else:
            raise ValueError, 'Invalid posterior.'

        
    def isValid(self,x):
        if not len(x) == self.p:
            raise InvalidDistributionInput, "\n\tInvalid data: wrong dimension(s) "+str(len(x))+" in MultiNormalDistribution(p="+str(self.p)+")."
        for v in x:
            try:
                float(v)
            except (ValueError):
                raise InvalidDistributionInput, "\n\tInvalid data: "+str(x)+" in MultiNormalDistribution."

class DependenceTreeDistribution(ConditionalGaussDistribution):
    """
    This class implemements a tree of conditional Gaussians, including the
    tree topology learning. 
    """
    def __init__(self, p, mu, w, sigma):
        """
        Constructor

        @param p: dimensionality of the distribution
        @param mu: mean parameter vector
        @param w: covariance weights (representing off-diagonal entries in the full covariance matrix)
        @param sigma: standard deviations (diagonal entries of the covariance matrix)
        """

        parents = self.randomStructure(p)
        # linear initialization of tree structure
        struct = {}
        for i in range(p-1):
          struct[i+1] = i
        struct[0]=-1
    	ConditionalGaussDistribution.__init__(self,p, mu, w, sigma, struct)


    def MStep(self,posterior,data,mix_pi=None):
        if isinstance(data,DataSet):
            x = data.internalData
        elif isinstance(data,numpy.ndarray):
            x = data
        else:
            raise TypeError, "Unknown/Invalid input to MStep."    

        post = posterior.sum() # sum of posteriors
        self.mu = numpy.dot(posterior, x)/post

        # centered input values (with new mus)
        centered = numpy.subtract(x,numpy.repeat([self.mu],len(x), axis=0)); 


        # estimating correlation factor
        sigma = numpy.dot(numpy.transpose(numpy.dot(numpy.identity(len(posterior))*posterior,centered)),centered)/post # sigma/covariance matrix

        diagsigma = numpy.diagflat(1.0/numpy.diagonal(sigma)) # vector with diagonal entries of sigma matrix
        correlation = numpy.dot(numpy.dot(diagsigma,numpy.multiply(sigma,sigma)),diagsigma) # correlation matrix with entries sigma_xy^2/(sigma^2_x * sigma^2_y)
        
        correlation = correlation - numpy.diagflat(numpy.diagonal(correlation)) # making diagonal entries = 0

        # XXX - check this
        parents = self.maximunSpanningTree(correlation) # return maximun spanning tree from the correlation matrix
        self.parents = self.directTree(parents,0) # by default direct tree from 0


        # XXX note that computational time could be saved as these functions share same suficient statistics
        ConditionalGaussDistribution.MStep(self,posterior,data,mix_pi)

         
    def maximunSpanningTree(self,weights):
        """
        Estimates the MST given a fully connected graph defined by the symetric matrix weights.
        using Prim`s algorithm.
        
        @param weights: edge weights
        """

        # start with an empty tree and random vertex  
        edgestree = {}   
        for i in range(self.p): 
          edgestree[i] = []
        verticestree = [0]
       
        while len(verticestree) < self.p:
            # possible edges = only ones form vertices at the current tree 
            candidates = weights[verticestree,:]    
          
            # look for maximal candidate edges       
            indices = numpy.argmax(candidates,1) # max neighboors in verticestrrees
            values = numpy.max(candidates,1) 
            uaux = numpy.argmax(values)
            u = verticestree[uaux] 
            v = indices[uaux]

            # add (u,v) att tree
            edgestree[u].append(v)
            edgestree[v].append(u)
            #edgestree[v] = u

            #zeroing all vertices between v and verticestree   
            for i  in  verticestree:
                weights[v,i] = 0
                weights[i,v] = 0

            # add (v) at tree
            verticestree.append(v)
	
	    return edgestree

    def directTree(self,tree,root):
        parent = {}
        queue = []
        # directing the tree from the root
        parent[root] = -1
        visited = numpy.zeros((self.p,1))      
        for u in tree[root]: 
            queue.append((root,u))
            visited[root] = 1
        while len(queue) > 0:
            (u,v) = queue.pop()
            parent[v] = u
            for newv in tree[v]:
               if not visited[newv]:
                   queue.append((v,newv))
            visited[v] = 1
	    return parent                                 
        
  
    def randomStructure(self,p):
        # linear initialization of tree structure
        struct = {}
        for i in range(p-1):
            struct[i+1] = i
        struct[0]=-1
        return struct

        
    def __str__(self):
        return 'Dependence Tree: \nmu='+str(self.mu)+' \nsigma='+str(self.sigma)+'\nw='+str(self.w)+'\nparents='+str(self.parents)

    
            
class ExponentialDistribution(ProbDistribution):
    """
    Exponential distribution
    """
    def __init__(self, lambd):
        """
        Constructor

        @param lambd: shape parameter lambda
        """
        
        self.p = self.suff_p = 1
        self.lambd = lambd  # lambd is a rate: 0.0 < lambd <= 1.0
        self.freeParams = 1
    
    def __copy__(self):
        return ExponentialDistribution(self.lambd)

    
    def __str__(self):
        return "Exponential:  ["+str(self.lambd)+"]"
    
    def __eq__(self, other):
        if not isinstance(other, ExponentialDistribution):
            return False
        if not self.lambd == other.lambd:
            return False
        return True 
                

    def pdf(self, data):
        if isinstance(data, DataSet ):  
            assert data.internalData is not None, "Internal data not initialized."
            nr = len(data.internalData)
            assert data.internalData.shape == (nr,1), 'shape = '+str(data.internalData.shape)
            
            x = numpy.transpose(data.internalData)[0]
            
        elif isinstance(data, numpy.ndarray): 
            nr = len(data)
            
            if data.shape == (nr,1):  # data format needs to be changed
                x = numpy.transpose(data)[0]
            elif data.shape == (nr,): 
                x = data
            else:
                raise TypeError,'Invalid data shape: '+str(data.shape)
        else:
            raise TypeError,"Unknown/Invalid input type:"+str(type(data))
        
        return math.log(self.lambd) + (-self.lambd * x)  # XXX pure Python implementation for now


    def sample(self):
        return random.expovariate(self.lambd)

    def MStep(self,posterior,data,mix_pi=None):       
        # data has to be reshaped for parameter estimation
        if isinstance(data,DataSet):
            x = data.internalData[:,0]
        elif isinstance(data,numpy.ndarray):
            x = data[:,0]
        else:
            raise TypeError, "Unknown/Invalid input to MStep."    
        
        self.lambd = posterior.sum() / numpy.dot(posterior, x)
        


    def isValid(self,x):
        "Checks whether 'x' is a valid argument for the distribution."
        try:
            float(x)
        except ValueError:
            #print "Invalid data: ",x,"in ExponentialDistribution."
            raise InvalidDistributionInput,"\n\tInvalid data: "+str(x)+" in ExponentialDistribution."

        if x < 0:
            raise InvalidDistributionInput,"\n\tInvalid data: negative float "+str(x)+" in ExponentialDistribution."
            
    def formatData(self,x):
        """
        """
        if type(x) == list and len(x) == 1:
            x = x[0]
        self.isValid(x)
        return [self.p,[x]]

    def flatStr(self,offset):
        offset +=1
        return "\t"*offset+";Exp;"+str(self.lambd)+"\n"

    def posteriorTraceback(self,x):
        return self.pdf(x)

class UniformDistribution(ProbDistribution):
    """
    Uniform distribution over a given intervall.
    """
    def __init__(self, start, end):
        """
        Constructor
        
        @param start: begin of interval
        @param end: end of interval
        """
        assert start < end

        self.p = self.suff_p = 1
        self.freeParams = 0
        
        self.start = start
        self.end = end
        self.density = numpy.log(1.0 / (end - start))   # compute log density value only once

    def __eq__(self,other):
        raise NotImplementedError

    def __str__(self):
        return "Uniform:  ["+str(self.start)+","+str(self.end)+"]"

    def __copy__(self):
        raise NotImplementedError
    
    def pdf(self,data):
        if isinstance(data, DataSet ):
            x = data.internalData
        elif isinstance(data, numpy.ndarray):
            x = data
        else:
            raise TypeError,"Unknown/Invalid input type." 
        res = numpy.zeros(len(x),dtype='Float64')
        for i in range(len(x)):
            # density is self.density inside the interval and -inf (i.e. 0) outside
            if self.start <= x[i][0] <= self.end:
                res[i] = self.density
            else:
                res[i] = float('-inf')

        return res
        
    def MStep(self,posterior,data,mix_pi=None):       
        # nothing to be done...
        pass
        
    def sample(self):
        return random.uniform( self.start, self.end)
        
        
    def sampleSet(self,nr):
        set = []
        for i in range(nr):
            set.append(self.sample())
        return set
            
    def isValid(self,x):
        try:
            float(x)
        except (ValueError):
            raise InvalidDistributionInput, "\n\tInvalid data in "+str(x)+" in UniformDistribution."
        
    def formatData(self,x):
        if isinstance(x,list) and len(x) == 1: 
            x = x[0]
        self.isValid(x)  # make sure x is valid argument
        return [self.p,[x]]

    
    def flatStr(self,offset):
        raise NotImplementedError, "Boom !"
                
    def posteriorTraceback(self,x):
        raise NotImplementedError, "Kawoom !"

    def update_suff_p(self):
        return self.suff_p

    def merge(self,dlist, weights):
        raise NotImplementedError, "Kawoom !"


# ------------------------------------------------------------------------
           
class MultinomialDistribution(ProbDistribution):
    """
    Multinomial Distribution
    """
    
    def __init__(self,p, M, phi, alphabet = None,parFix = None):
        """
        Constructor

        @param M: number of possible outcomes (0 to M-1)
        @param p: number of values in each sample
        @param phi:= discrete distribution of N objects
        @param alphabet: Alphabet object (optional)
        @param parFix: list of flags to determine if any elements of phi should be fixed
        """
        assert len(phi) == M, "Invalid number of parameters for MultinomialDistribution."
        assert abs((1.0 - sum(phi))) < 1e-12, str(phi)+": "+ str(1.0 - sum(phi)) # check parameter validity

        self.p = p # lenght of input vectors, corresponds to p in MixtureModel
        self.M = M
        self.suff_p = M  # length of the sufficient statistics, equal to size of alphabet
        
        # in case there is no alphabet specified IntegerRange is used
        if alphabet:
            assert len(alphabet) == self.M, "Size of alphabet and M does not match: "+str(len(alphabet))+" != "+str(self.M)
            self.alphabet = alphabet
        else:
            self.alphabet = IntegerRange(0, self.M)    

        if parFix == None:
            self.parFix = numpy.array([0] * self.M)
        else:
            assert len(parFix) == self.M, "Invalid length of parFix vector."
            self.parFix = numpy.array(parFix)
            
        # number of free parameters is M-1 minus the number of fixed entries in phi
        self.freeParams = M-1 - sum(self.parFix)

        self.phi = numpy.array(phi,dtype='Float64')
       
        # minimal value for any component of self.phi, enforced in MStep
        self.min_phi =  ( 1.0 / self.M ) * 0.001

    def __eq__(self,other):
        res = False
        if isinstance(other,MultinomialDistribution):
            if  other.p == self.p and other.M == self.M and numpy.allclose( other.phi,self.phi):
                res = True
        return res        
    
    def __copy__(self):
        "Interface for the copy.copy function"
        return MultinomialDistribution(self.p,self.M,copy.deepcopy(self.phi),self.alphabet,parFix=self.parFix)
    
    def __str__(self):
        outstr = "Multinom(M = "+ str(self.M)+", N = "+ str(self.p) +" ) : " + str(self.phi) #+"\n"
        #outstr += str(self.alphabet) + "\n"
        return outstr
    
    def pdf(self, data):
        # Note: The multinomial coefficient is omitted in the implementation.
        # Result is proportional to the true log densitiy which is sufficient for 
        # the EM.
        # gsl computes the true density, including the multinomial coefficient normalizing constant
        # therefore it is less efficient than the implementation below
        if isinstance(data, DataSet ):
            x = data.internalData
        elif isinstance(data, numpy.ndarray):
            x = data
        else:
            raise TypeError,"Unknown/Invalid input type." 

        # switch to log scale for density computation
        log_phi = numpy.log(self.phi) 

        # computing un-normalized density
        res = numpy.zeros(len(x),dtype='Float64')
        for j in range(len(x)):
            for i in range(self.M):
                res[j] +=  (log_phi[i] * x[j,i])

        return res

    def sample(self):
        sample = []
        for i in range(self.p): 
            sum = 0.0
            p = random.random()
            for k in range(self.M):
                sum += self.phi[k]
                if sum >= p:
                    break
            sample.append(k)

        return map(self.alphabet.external,sample)   
    
    def sampleSet(self,nr):
        res = []
        
        for i in range(nr):
            res.append(self.sample())
            
        return res    
    
    def MStep(self,posterior,data,mix_pi=None):       
        if isinstance(data,DataSet):
            x = data.internalData
        elif isinstance(data,numpy.ndarray):
            x = data
        else:
            raise TypeError, "Unknown/Invalid input to MStep."    

        ind = numpy.where(self.parFix == 0)[0]
        fix_flag = 0
        fix_phi = 1.0     
        dsum = 0.0
        
        # reestimating parameters
        for i in range(self.M):
            if self.parFix[i] == 1:
               fix_phi -= self.phi[i] 
               fix_flag = 1
               continue
            else:    
                est = numpy.dot(x[:,i], posterior)
                self.phi[i] = est
                dsum += est

        if dsum == 0.0: 
            raise InvalidPosteriorDistribution, "Invalid posterior in MStep."
        
        # normalzing parameter estimates
        self.phi[ind] = (self.phi[ind] * fix_phi) / dsum

        adjust = 0  # adjusting flag
        for i in range(self.M):
            if self.parFix[i] == 0 and self.phi[i] < self.min_phi:  
                adjust = 1
                self.phi[i] =  self.min_phi
        
        # renormalizing the adjusted parameters if necessary
        if adjust:
            dsum = sum(self.phi[ind])
            self.phi[ind] = (self.phi[ind] * fix_phi) / dsum

    def isValid(self,x):
        if sum(map(self.alphabet.isAdmissable,x)) != self.p:
            raise InvalidDistributionInput, "\n\tInvalid data: "+str(x)+" in MultinomialDistribution("+str(self.alphabet.listOfCharacters)+")."

    def formatData(self,x):
        count = [0] * self.M #  numpy.zeros(self.M)        

        # special case of p = 1
        if len(x) == 1:
            self.isValid( str(x[0]) )
            count[self.alphabet.internal(str(x[0]))] = 1
                        
            return [self.M, count]

        for i in range(self.M):
            self.isValid(x)
            f = lambda x: x == self.alphabet.listOfCharacters[i]
            count[i] = sum(map(f,x))
       
        return [self.M, count] 

    def flatStr(self,offset):
        offset +=1
        return "\t"*offset+";Mult;"+str(self.p)+";"+str(self.M)+";"+str(self.phi.tolist())+";"+str(self.alphabet.listOfCharacters)+";"+str(self.parFix.tolist())+"\n"

    def posteriorTraceback(self,x):
        return self.pdf(x)

    def merge(self,dlist, weights):
        raise DeprecationWarning, 'Part of the outdated structure learning implementation.'
        
        norm = sum(weights)
        m_phi = self.phi * weights[0]
        for i in range(len(dlist)):
            assert isinstance(dlist[i],MultinomialDistribution)
            assert dlist[i].p == self.p
            assert dlist[i].M == self.M
            for j in range(self.M):
                m_phi[j] += weights[i+1] * dlist[i].phi[j]

        f = lambda x: x/norm
        self.phi = numpy.array(map(f, m_phi) , dtype='Float64'   )
 


class DiscreteDistribution(MultinomialDistribution):
    """
    This is the special case of a MultinomialDistribution with p = 1, that is a simple univariate discrete
    distribution. Certain key functions are overloaded for increased efficiency.
    """
    def __init__(self, M, phi, alphabet = None,parFix = None):
        """
        @param M: size of alphabet
        @param phi: distribution parameters
        @param alphabet: Alphabet object (optional)
        @param parFix: list of flags to determine if any elements of phi should be fixed
        """
        
        MultinomialDistribution.__init__(self, 1, M, phi, alphabet = alphabet, parFix = parFix)
        self.suff_p = 1
        
    def __str__(self):
        outstr = "DiscreteDist(M = "+ str(self.M)+"): " + str(self.phi) #+"\n"
        #outstr += str(self.alphabet) + "\n"
        return outstr

    def __copy__(self):
        return DiscreteDistribution(self.M,copy.deepcopy(self.phi),self.alphabet,parFix=self.parFix)
    
    def pdf(self, data):
        if isinstance(data, DataSet ):
            assert data.p == 1
            x = data.getInternalFeature(0)
        elif isinstance(data, numpy.ndarray):
            x = data
        else:
            raise TypeError,"Unknown/Invalid input type." 

        # switch to log scale for density computation
        log_phi = numpy.log(self.phi) 

        # computing un-normalized density
        res = log_phi[x[:,0].astype('Int32')]
        return res

    def MStep(self,posterior,data,mix_pi=None):       
        if isinstance(data,DataSet):
            x = data.internalData
        elif isinstance(data,numpy.ndarray):
            x = data
        else:
            raise TypeError, "Unknown/Invalid input to MStep."    

        ind = numpy.where(self.parFix == 0)[0]
        fix_flag = 0
        fix_phi = 1.0     
        dsum = 0.0
        # reestimating parameters
        for i in range(self.M):
            if self.parFix[i] == 1:
               fix_phi -= self.phi[i] 
               fix_flag = 1
               continue
            else:    
                i_ind = numpy.where(x == i)[0]
                est = numpy.sum(posterior[i_ind])
                self.phi[i] = est
                dsum += est

        if dsum == 0.0: 
            print self
            print posterior
            
            raise InvalidPosteriorDistribution, "Invalid posterior in MStep."
        
        # normalzing parameter estimates
        self.phi[ind] = (self.phi[ind] * fix_phi) / dsum
        
        adjust = 0  # adjusting flag
        for i in range(self.M):
            if self.parFix[i] == 0 and self.phi[i] < self.min_phi:  
                #print "---- enforcing minimal phi -----"
                adjust = 1
                self.phi[i] =  self.min_phi
        
        # renormalizing the adjusted parameters if necessary
        if adjust:
            dsum = sum(self.phi[ind])
            self.phi[ind] = (self.phi[ind] * fix_phi) / dsum

    def sample(self):
        for i in range(self.p): 
            sum = 0.0
            p = random.random()
            for k in range(self.M):
                sum += self.phi[k]
                if sum >= p:
                    break
        return self.alphabet.external(k)
    
    def sampleSet(self,nr):
        res = []
        for i in range(nr):
            res.append(self.sample())
        return res    

    def sufficientStatistics(self, posterior, data):
        stat = numpy.zeros(self.M,dtype='Float64')
        for i in range(self.M):
            i_ind = numpy.where(data == i)[0]
            stat[i] = numpy.sum(posterior[i_ind])
        return stat


    def formatData(self,x):
        self.isValid(x)
        if type(x) == list:
            assert len(x) == 1
            internal = self.alphabet.internal( str(x[0]))
        else:
              internal = self.alphabet.internal( str(x) )  
        return [1, [internal] ]
    
    def isValid(self,x):
        if type(x) == str or type(x) == int or type(x) == float: 
            if not self.alphabet.isAdmissable(str(x)):
                raise InvalidDistributionInput, "\n\tInvalid data: "+str(x)+" in DiscreteDistribution("+str(self.alphabet.listOfCharacters)+")."
        else:
            if type(x) == list and len(x) == 1:
                self.isValid(x[0])
            else:            
                raise InvalidDistributionInput, "\n\tInvalid data: "+str(x)+" in DiscreteDistribution("+str(self.alphabet.listOfCharacters)+")."

    def flatStr(self,offset):
        offset +=1
        return "\t"*offset+";Discrete;"+str(self.M)+";"+str(self.phi.tolist())+";"+str(self.alphabet.listOfCharacters)+";"+str(self.parFix.tolist())+"\n"



class DirichletPrior(PriorDistribution):  # DirichletDistribution, 
    """
    Dirichlet distribution as Bayesian prior for MultinomialDistribution and derived .
    """
    def __init__(self, M ,alpha):
        """
        @param M: number of dimensions
        @param alpha: distribution parameters
        """
        assert M == len(alpha)
        #for a in alpha:
        #    assert a > 0.0, "Invalid parameter."
        
        self.M = M
        self.alpha = numpy.array(alpha, dtype='Float64' )
        self.alpha_sum = numpy.sum(alpha) # assumes alphas to remain constant !
        self.p = M
        self.suff_p = M
        self.freeParams = M
        
        self.constant_hyperparams = 1  # hyperparameters are constant

    def __copy__(self):
        cp_alpha = copy.deepcopy(self.alpha)
        return DirichletPrior(self.M,cp_alpha)

    def __str__(self):
        return "DirichletPrior: "+str(self.alpha)
 
    def __eq__(self,other):
        if isinstance(other,DirichletPrior):
            if self.M == other.M and numpy.alltrue( self.alpha == other.alpha ):
                return True
            else:
                return False    
        else:
            return False

    def sample(self):
        """
        Samples from Dirichlet distribution
        """
        phi = _C_mixextend.wrap_gsl_dirichlet_sample(self.alpha,self.M)

        d = DiscreteDistribution(self.M, phi)
        return d
        

    def pdf(self,m): 

        # XXX should be unified ...
        if isinstance(m, MultinomialDistribution):
            # use GSL implementation
            #res = pygsl.rng.dirichlet_lnpdf(self.alpha,[phi])[0] XXX
            try:
                res = _C_mixextend.wrap_gsl_dirichlet_lnpdf(self.alpha, [m.phi])
            except ValueError:
                print m
                print self
                raise
            
            return res[0]

        elif isinstance(m, list):
            in_l = [ d.phi for d in m ]
            # use GSL implementation
            res = _C_mixextend.wrap_gsl_dirichlet_lnpdf(self.alpha, in_l)
            return res
        else:
            raise TypeError            
        
    def posterior(self,m,x):  
        """
        Returns the posterior for multinomial distribution 'm' for multinomial count data 'x'
        The posterior is again Dirichlet.
        """
        assert isinstance(m, MultinomialDistribution)
        res = numpy.ones(len(x),dtype='Float64')
        for i,d in enumerate(x):
            post_alpha = self.alpha + d
            res[i] = numpy.log( _C_mixextend.wrap_gsl_dirichlet_pdf(post_alpha,[m.phi]))

        return res

    def marginal(self,x):
        """ 
        Returns the log marginal likelihood of multinomial counts 'x' (sufficient statistics)
        with Dirichlet prior 'self' integrated over all parameterizations of the multinomial.    
        """
        # XXX should be eventually replaced by more efficient implementation
        # in Dirchlet mixture prior paper (K. Sjoelander,Karplus,..., D.Haussler) 
         
        x_sum = sum(x)

        term1 =  _C_mixextend.wrap_gsl_sf_lngamma(self.alpha_sum) - _C_mixextend.wrap_gsl_sf_lngamma(self.alpha_sum + x_sum)
        term2 = 0.0
        for i in range(self.p):
            term2 +=  _C_mixextend.wrap_gsl_sf_lngamma(self.alpha[i] + x[i] ) -  _C_mixextend.wrap_gsl_sf_lngamma(self.alpha[i])

        res = term1 + term2
        return res

    def mapMStep(self, dist, posterior, data, mix_pi=None, dist_ind = None):
        # Since DiscreteDistribution is a special case of MultinomialDistribution
        # the DirichletPrior applies to both. Therefore we have to distinguish the
        # two cases here. The cleaner alternative would be to derive specialized prior
        # distributions but that would be unnecessarily complicated at this point.
        if isinstance(dist, DiscreteDistribution):
            ind = numpy.where(dist.parFix == 0)[0]
            fix_phi = 1.0     
            dsum = 0.0
            for i in range(dist.M):
                if dist.parFix[i] == 1:
                   fix_phi -= dist.phi[i] 
                   continue
                else:
                    i_ind = numpy.where(data == i)[0]
                    est = numpy.sum(posterior[i_ind]) + self.alpha[i] -1 
                    dist.phi[i] = est
                    dsum += est

            # normalizing parameter estimates
            dist.phi[ind] = (dist.phi[ind] * fix_phi) / dsum
        elif isinstance(dist, MultinomialDistribution):

            fix_phi = 1.0     
            dsum = 0.0
            # reestimating parameters
            for i in range(dist.M):
                if dist.parFix[i] == 1:
                   #print "111"
                   fix_phi -= dist.phi[i] 
                   continue
                else:       
                    est = numpy.dot(data[:,i], posterior) + self.alpha[i] -1 
                    dist.phi[i] = est
                    dsum += est
    
            if dsum == 0.0: 
                raise InvalidPosteriorDistribution, "Invalid posterior in MStep."
    
            ind = numpy.where(dist.parFix == 0)[0]
            # normalzing parameter estimates
            dist.phi[ind] = (dist.phi[ind] * fix_phi) / dsum
    
        else:
            raise TypeError,'Invalid input '+str(dist.__class__)


    def mapMStepMerge(self, group_list):
        #XXX only for DiscreteDistribution for now, MultinomialDistribution to be done
        assert isinstance(group_list[0].dist, DiscreteDistribution), 'only for DiscreteDistribution for now'

        pool_req_stat = copy.copy(group_list[0].req_stat)
        pool_post_sum = group_list[0].post_sum
        pool_pi_sum = group_list[0].pi_sum
        
        for i in range(1,len(group_list)):
            pool_req_stat += group_list[i].req_stat
            pool_post_sum += group_list[i].post_sum
            pool_pi_sum += group_list[i].pi_sum
        

        new_dist = copy.copy(group_list[0].dist)  # XXX copy necessary ?

        ind = numpy.where(group_list[0].dist.parFix == 0)[0]
        fix_phi = 1.0     
        dsum = 0.0
        for i in range(group_list[0].dist.M):
            if group_list[0].dist.parFix[i] == 1:
               assert group_list[1].dist.parFix[i] == 1  # incomplete consistency check of parFix (XXX)
               
               fix_phi -= new_dist.phi[i] 
               continue
            else:
                est = pool_req_stat[i] + self.alpha[i] -1 
                new_dist.phi[i] = est
                
                dsum += est

        # normalizing parameter estimates
        new_dist.phi[ind] = (new_dist.phi[ind] * fix_phi) / dsum
            
        return CandidateGroup(new_dist, pool_post_sum, pool_pi_sum, pool_req_stat )


    def mapMStepSplit(self, toSplitFrom, toBeSplit):
        #XXX only for DiscreteDistribution for now, MultinomialDistribution to be done
        assert isinstance(toSplitFrom.dist, DiscreteDistribution), 'only for DiscreteDistribution for now'

        split_req_stat = copy.copy(toSplitFrom.req_stat)
        split_req_stat -= toBeSplit.req_stat
        
        split_post_sum = toSplitFrom.post_sum - toBeSplit.post_sum
        split_pi_sum = toSplitFrom.pi_sum - toBeSplit.pi_sum
        
        

        new_dist = copy.copy(toSplitFrom.dist)  # XXX copy necessary ?

        ind = numpy.where(toSplitFrom.dist.parFix == 0)[0]
        fix_phi = 1.0     
        dsum = 0.0
        for i in range(toSplitFrom.dist.M):
            if toSplitFrom.dist.parFix[i] == 1:
               
               fix_phi -= new_dist.phi[i] 
               continue
            else:
                est = split_req_stat[i] + self.alpha[i] -1 
                new_dist.phi[i] = est
                dsum += est

        # normalizing parameter estimates
        new_dist.phi[ind] = (new_dist.phi[ind] * fix_phi) / dsum
            
        return CandidateGroup(new_dist, split_post_sum, split_pi_sum, split_req_stat )

        
    def isValid(self,x):
        if not isinstance(x,MultinomialDistribution):
            raise InvalidDistributionInput, "in DirichletPrior: "+str(x)
        else:
            if self.M != x.M:
                raise InvalidDistributionInput, "in DirichletPrior: "+str(x)

    def flatStr(self,offset):
        offset +=1
        return "\t"*offset+";DirichletPr;"+str(self.M)+";"+str(self.alpha.tolist())+"\n"



class NormalGammaPrior(PriorDistribution): 
    """
    Inverse-Gamma Normal distribution prior for univariate Normal distribution.
    """

    def __init__(self, mu, kappa, dof, scale  ):
        """
        Constructor
        
        @param mu: hyper-parameter mu
        @param kappa: hyper-parameter kappa 
        @param dof: hyper-parameter dof
        @param scale: hyper-parameter scale
        """
        # parameters of the Normal prior on mu | sigma
        self.mu_p = float(mu)
        self.kappa = float(kappa)
        
        # parameters on the inverse-gamma prior on sigma
        self.dof = float(dof)
        self.scale = float(scale)

        self.constant_hyperparams = 1  # hyperparameters are constant

    def __str__(self):
        outstr = "NormalGamma: mu_p="+str(self.mu_p)+", kappa="+str(self.kappa)+", dof="+str(self.dof)+", scale="+str(self.scale)
        return outstr

    def __eq__(self,other):
        if not isinstance(other, NormalGammaPrior):
            return False
        if self.mu_p != other.mu_p:
            return False    
        if self.kappa != other.kappa:
            return False    
        if self.dof != other.dof:
            return False    
        if self.scale != other.scale:
            return False    
        return True

    def __copy__(self):
        return NormalGammaPrior(self.mu_p, self.kappa, self.dof, self.scale )
        
    
    def pdf(self, n):    

        if isinstance(n,NormalDistribution):
            res = _C_mixextend.get_log_normal_inverse_gamma_prior_density( self.mu_p, self.kappa, self.dof, self.scale,[n.mu], [n.sigma] )[0]
            return res        

        elif isinstance(n,list):
            # extract parameters, XXX better way to do this ?
            d_sigma = numpy.zeros(len(n))                
            d_mu = numpy.zeros(len(n))
            for i,d in enumerate(n):
                d_sigma[i] = d.sigma
                d_mu[i] = d.mu
            
            # call to extension function
            return _C_mixextend.get_log_normal_inverse_gamma_prior_density( self.mu_p, self.kappa, self.dof, self.scale, d_mu, d_sigma )
        else:
            raise TypeError

    def posterior(self,m,x):  
        raise NotImplementedError, "Needs implementation"

    def marginal(self,x):
        raise NotImplementedError, "Needs implementation"

    def mapMStep(self, dist, posterior, data, mix_pi=None, dist_ind = None):

        assert isinstance(dist,NormalDistribution)
        # data has to be reshaped for parameter estimation
        if isinstance(data,DataSet):
            x = data.internalData[:,0]
        elif isinstance(data,numpy.ndarray):
            x = data[:,0]
        else:
            raise TypeError, "Unknown/Invalid input to MStep."    
        nr = len(x)
        sh = x.shape
        assert sh == (nr,)
        
        post_sum = numpy.sum(posterior)  # n_k
        if post_sum == 0.0:
            print dist
            raise InvalidPosteriorDistribution, "Sum of posterior is zero."

        # computing ML estimates for mu and sigma
        ml_mu = numpy.dot(posterior, x) / post_sum  # ML estimator for mu
        new_mu =  ( (post_sum*ml_mu) + (self.kappa * self.mu_p) )/ ( post_sum + self.kappa)

        n_sig_num_1 = self.scale + ( (self.kappa * post_sum) / ( self.scale+post_sum ) ) * (ml_mu - self.mu_p)**2
        n_sig_num_2 = numpy.dot(posterior, (x - ml_mu)**2 )
        n_sig_num = n_sig_num_1 + n_sig_num_2
        n_sig_denom = self.dof + post_sum + 3.0

        new_sigma = math.sqrt( n_sig_num / n_sig_denom)

        # assigning updated parameter values
        dist.mu = new_mu
        dist.sigma = new_sigma
     

    def mapMStepMerge(self, group_list):

        pool_req_stat = copy.copy(group_list[0].req_stat)
        pool_post_sum = group_list[0].post_sum
        pool_pi_sum = group_list[0].pi_sum
        
        for i in range(1,len(group_list)):
            pool_req_stat += group_list[i].req_stat
            pool_post_sum += group_list[i].post_sum
            pool_pi_sum += group_list[i].pi_sum
        
        
        new_mu = (pool_req_stat[0] + (self.kappa * self.mu_p) ) / ( pool_post_sum + self.kappa)

        y = (pool_req_stat[0] ) / ( pool_post_sum)
        n_sig_num_1 = self.scale + ( (self.kappa * pool_post_sum) / ( self.scale+pool_post_sum ) ) * (y - self.mu_p)**2
        n_sig_num_2 = (pool_req_stat[1]) - 2 * y * (pool_req_stat[0]) + y**2 * pool_post_sum
        n_sig_num = n_sig_num_1 + n_sig_num_2 
        n_sig_denom = self.dof + pool_post_sum + 3.0
        
        new_sigma = math.sqrt( n_sig_num / n_sig_denom)
        new_dist = NormalDistribution(new_mu, new_sigma) 
        
        return CandidateGroup(new_dist, pool_post_sum, pool_pi_sum, pool_req_stat )

    def mapMStepSplit(self, toSplitFrom, toBeSplit):

        split_req_stat = copy.copy(toSplitFrom.req_stat)
        split_req_stat -= toBeSplit.req_stat
        
        split_post_sum = toSplitFrom.post_sum - toBeSplit.post_sum
        split_pi_sum = toSplitFrom.pi_sum - toBeSplit.pi_sum
        
        new_mu = (split_req_stat[0] + (self.kappa * self.mu_p) ) / ( split_post_sum + self.kappa)

        y = (split_req_stat[0] ) / ( split_post_sum)
        n_sig_num_1 = self.scale + ( (self.kappa * split_post_sum) / ( self.scale+split_post_sum ) ) * (y - self.mu_p)**2
        n_sig_num_2 = (split_req_stat[1]) - 2 * y * (split_req_stat[0]) + y**2 * split_post_sum
        n_sig_num = n_sig_num_1 + n_sig_num_2 
        n_sig_denom = self.dof + split_post_sum + 3.0
        new_sigma = math.sqrt( n_sig_num / n_sig_denom)
        
        new_dist = NormalDistribution(new_mu, new_sigma) 
        
        return CandidateGroup(new_dist, split_post_sum, split_pi_sum, split_req_stat )

    def flatStr(self,offset):
        offset +=1
        return "\t"*offset+";NormalGamma;"+str(self.mu_p)+";"+str(self.kappa)+";"+str(self.dof)+";"+str(self.scale)+"\n"
        
    def isValid(self,x):
        if not isinstance(x,NormalDistribution):
            raise InvalidDistributionInput, "NormalGammaPrior: "+ str(x)

    def setParams(self, x, K):
        """
        Get guesses for hyper-parameters according to the heuristics used in "Bayesian Regularization for Normal
        Mixture Estimation and Model-Based Clustering" (C.Fraley and A.E. Raftery)
        
        @param x: numpy data vector
        @param K: number of components
        """
        nr = len(x)
        assert x.shape == (nr,1)
        x = x[:,0]

        self.mu_p = x.mean()
        self.kappa = 0.01       
        self.dof = 3.0
        
        self.scale = x.var() / (K**2)


class DirichletMixturePrior(PriorDistribution):
    """
    Mixture of Dirichlet distributions prior for multinomial data.
    """     
    def __init__(self, G, M, pi, dComp):
        """
        @param G: number of components
        @param M: dimensions of component Dirichlets
        @param pi: mixture weights
        @param dComp: list of DirichletPrior distributions
        """
        assert len(dComp) == len(pi) == G
        for d in dComp:
            assert d.M == M

        self.G = G
        self.M = M
        self.pi = numpy.array(pi,dtype='Float64')
        self.log_pi = numpy.log(self.pi)  # assumes pi is not changed from the outside (XXX accesor functions ?)
        self.dComp = dComp
        self.constant_hyperparams = 1  # hyperparameters are constant

    def __str__(self): 
        s = ['DirichletMixturePrior( G='+str(self.G)+' )']
        s.append( '\tpi='+str(self.pi)+'\n' ) 
        for i in range(self.G):
            s.append('\t\t'+str(self.dComp[i])+'\n')
        return ''.join(s)

    def __eq__(self,other):
        if not isinstance(other,DirichletMixturePrior):
            return False
        if self.G != other.G or self.M != other.M:
            return False
        if not numpy.alltrue(other.pi == self.pi):
            return False
        for i,d1 in enumerate(self.dComp):
            if not d1 == other.dComp[i]:
                return False
        return True                

    def __copy__(self):
        cp_pi = copy.deepcopy(self.pi)
        cp_dC = [ copy.deepcopy( self.dComp[i] ) for i in range(self.G) ]
        return DirichletMixturePrior(self.G,self.M,cp_pi,cp_dC)

    def pdf(self, m):    
        if isinstance(m, MultinomialDistribution):  # XXX debug
            logp_list = numpy.zeros(self.G,dtype='Float64')
            for i in range(self.G):
                logp_list[i] = self.log_pi[i]  + self.dComp[i].pdf(m) 
            res =  sumlogs(logp_list)
            return res
        
        elif type(m) == list:
            logp_mat = numpy.zeros((self.G, len(m)))
            for i in range(self.G):
                logp_mat[i,:] = self.dComp[i].pdf(m)
           
            for i in range(len(m)):  # XXX slow
                logp_mat[:,i] += self.log_pi
            
            res = _C_mixextend.matrix_sum_logs(logp_mat)
            return res
        else:
            raise TypeError


    def marginal(self,dist,posterior,data):
        suff_stat = numpy.zeros(self.M,dtype='Float64')
        if isinstance(dist,DiscreteDistribution):
            for i in range(self.M):
                i_ind = numpy.where(data == i)[0]
                suff_stat[i] = numpy.sum(posterior[i_ind])
        elif isinstance(dist, MultinomialDistribution):
            for i in range(self.M):
                suff_stat[i] = numpy.dot(data[:,i], posterior)
        else:
            raise TypeError,'Invalid input '+str(dist.__class__)

        res = 0.0
        for i in range(self.G):
            res += self.dComp[i].marginal(suff_stat) + numpy.log(self.pi[i])
        return res



    def posterior(self,dist):  
        """
        Component membership posterior distribution of MultinomialDistribution 'dist'.

        @param dist: MultinomialDistribution object

        @return: numpy of length self.G containing the posterior of component membership
        """
        prior_post = numpy.array([ dirich.pdf(dist) + self.log_pi[i] for i,dirich in enumerate(self.dComp) ], dtype='Float64')        
        log_sum = sumlogs(prior_post)
        prior_post -= log_sum

        prior_post = numpy.exp(prior_post)
        return prior_post

    def mapMStep(self, dist, posterior, data, mix_pi=None, dist_ind = None):
        suff_stat = numpy.zeros(self.M,dtype='Float64')
        if isinstance(dist,DiscreteDistribution):
            for i in range(self.M):
                i_ind = numpy.where(data == i)[0]
                suff_stat[i] = numpy.sum(posterior[i_ind])

        elif isinstance(dist, MultinomialDistribution):
            for i in range(self.M):
                suff_stat[i] = numpy.dot(data[:,i], posterior)
        else:
            raise TypeError,'Invalid input '+str(dist.__class__)

        # posterior of the given multinomial distribution 'dist'
        # with respect to the components of the Dirichlet mixture prior
        prior_post = self.posterior(dist)

        fix_flag = 0    
        fix_phi = 1.0     
        dsum = 0.0
        for i in range(self.M):
            if dist.parFix[i] == 1:  # checking for fixed entries in phi
                fix_flag = 1  
                fix_phi -= dist.phi[i] 
            else: # updating phi[i]
                e = numpy.zeros(self.G,dtype='Float64')
                for k in range(self.G):
                    e[k] = (suff_stat[i] + self.dComp[k].alpha[i] ) / ( sum(suff_stat) + self.dComp[k].alpha_sum  )

                est = numpy.dot(prior_post, e)
                dist.phi[i] = est
                dsum += est

        # re-normalizing parameter estimates if necessary
        if fix_flag:
            ind = numpy.where(dist.parFix == 0)[0]
            dist.phi[ind] = (dist.phi[ind] * fix_phi) / dsum

    def mapMStepMerge(self, group_list):

        new_dist = copy.copy(group_list[0].dist)

        prior_post = self.posterior(group_list[0].dist) 
        pool_req_stat = copy.copy(group_list[0].req_stat)
        pool_post_sum = group_list[0].post_sum
        pool_pi_sum = group_list[0].pi_sum
        
        for i in range(1,len(group_list)):
            pool_req_stat += group_list[i].req_stat
            pool_post_sum += group_list[i].post_sum
            pool_pi_sum += group_list[i].pi_sum
            
            prior_post += self.posterior(group_list[i].dist) * group_list[i].pi_sum
            
        prior_post = prior_post / pool_pi_sum  
        
        assert 1.0 - sum(prior_post) < 1e-10, str(prior_post)+' , '+str(sum(prior_post))+', '+ str(1.0 - sum(prior_post))  # XXX debug

        fix_flag = 0    
        fix_phi = 1.0     
        dsum = 0.0
        for i in range(self.M):
            if new_dist.parFix[i] == 1:  # assumes parFix is consistent 
                fix_flag = 1  
                fix_phi -=  new_dist.phi[i] 
            else: # updating phi[i]
                e = numpy.zeros(self.G,dtype='Float64')
                for k in range(self.G):
                    e[k] = (pool_req_stat[i] + self.dComp[k].alpha[i]) / ( pool_post_sum + self.dComp[k].alpha_sum  )

                est = numpy.dot(prior_post, e)
                new_dist.phi[i] = est
                dsum += est
                
                #print i,est

        # re-normalizing parameter estimates if necessary
        if fix_flag:
            ind = numpy.where(new_dist.parFix == 0)[0]
            new_dist.phi[ind] = (new_dist.phi[ind] * fix_phi) / dsum

        return CandidateGroup(new_dist, pool_post_sum, pool_pi_sum, pool_req_stat )


    def isValid(self,x):
        if not isinstance(x,MultinomialDistribution):
            raise InvalidDistributionInput, "DirichletMixturePrior: " + str(x)
            
        if x.M != self.M:
            raise InvalidDistributionInput, "DirichletMixturePrior: unequal dimensions " + str(x.M) + " != "+str(self.M)
            

    def flatStr(self,offset):
        offset +=1
        s = "\t"*offset + ";DirichMixPrior;"+str(self.G)+";"+str(self.M)+";"+str(self.pi.tolist())+"\n"
        for d in self.dComp:
            s+= d.flatStr(offset)
        return s


class ConditionalGaussPrior(PriorDistribution): 
    """
    Prior over ConditionalGaussDistribution. Assumes Normal prior over the covariance parameters w.
    
    """
    
    def __init__(self, nr_comps, p):
        """
            Constructor
            
            @param nr_comps: number of components in the mixture the prior is applied to
            @param p:  number of features in the ConditionalGaussDistribution the prior is applied to
        """
        
        self.constant_hyperparams = 0  # hyperparameters are updated as part of the mapEM
        self.nr_comps = nr_comps    # number of components in the mixture the prior is applied to
        self.p = p   # number of features in the ConditionalGaussDistribution the prior is applied to
        
        # no initial value needed, is updated as part of EM in updateHyperparameters
        self.beta =  numpy.zeros((self.nr_comps, self.p)) 
        self.nu =  numpy.zeros((self.nr_comps, self.p)) 

        # XXX initialization of sufficient statistics, necessary for hyperparameter updates
        self.post_sums = numpy.zeros(self.nr_comps)
        self.var = numpy.zeros( (self.nr_comps,self.p) )
        self.cov = numpy.zeros((self.nr_comps,self.p))
        self.mu = numpy.zeros((self.nr_comps,self.p))


    def __str__(self):
        return 'ConditionalGaussPrior(beta='+str(self.beta)+')'


    def pdf(self, d):    
        if type(d) == list:
            N = numpy.sum(self.post_sums)

            res = numpy.zeros( len(d))
            for i in range(len(d)):
                for j in range(1,d[i].p):
                    pid = d[i].parents[j]
                    res[i] += (1.0/self.cov[i,j]**2) / (self.nu[i,j] * (self.post_sums[i]/N))
                    res[i] +=  numpy.log(_C_mixextend.wrap_gsl_ran_gaussian_pdf(0.0, 
                                         math.sqrt((self.beta[i,j] * self.cov[i,j]**2) / (self.var[i,pid] * (self.post_sums[i]/N) )), 
                                         [d[i].w[j]] ))
        else:
            raise TypeError, 'Invalid input '+str(type(d))

        return res        


    def posterior(self,m,x):  
        raise NotImplementedError, "Needs implementation"

    def marginal(self,x):
        raise NotImplementedError, "Needs implementation"


    def mapMStep(self, dist, posterior, data, mix_pi=None, dist_ind=None):    
        assert not dist_ind == None # XXX debug

        post_sum = numpy.sum(posterior)
        self.post_sums[dist_ind] = post_sum

        # checking for valid posterior: if post_sum is zero, this component is invalid
        # for this data set
        if post_sum != 0.0:

            # reestimate mu
            for j in range(dist.p):
                # computing ML estimates for w and sigma
                self.mu[dist_ind,j] = numpy.dot(posterior, data[:,j]) / post_sum
                #self.var[dist_ind,j] = numpy.dot(posterior, (data[:,j] - dist.mu[j])**2 ) / post_sum
                self.var[dist_ind,j] = numpy.dot(posterior, (data[:,j] - self.mu[dist_ind,j])**2 ) / post_sum

                if j > 0:  # w[0] = 0.0 is fixed
                    pid = dist.parents[j]
                    self.cov[dist_ind,j] = numpy.dot(posterior, (data[:,j] - self.mu[dist_ind,j]) * (data[:,pid] - self.mu[dist_ind,pid])) / post_sum 
                    
                    # update hyperparameters beta
                    self.beta[dist_ind,j] = post_sum / ( (( self.var[dist_ind,j] * self.var[dist_ind,pid]) / self.cov[dist_ind,j]**2)  -1 )

                    # update hyperparameters nu
                    self.nu[dist_ind,j] = - post_sum / (2 * dist.sigma[j]**2)

                    # update regression weights
                    dist.w[j] = self.cov[dist_ind,j] / (dist.sigma[pid]**2 * (1+ self.beta[dist_ind,j]**-1 ) )
                    
                    # update standard deviation
                    dist.sigma[j] = math.sqrt( self.var[dist_ind,j] - (dist.w[j]**2 * dist.sigma[pid]**2 * (1+ (1.0/self.beta[dist_ind,j])) ) - self.nu[dist_ind,j]**-1  ) 
                    # update means
                    dist.mu[j] =  self.mu[dist_ind,j] #- (dist.w[j] * self.mu[dist_ind,pid])

                else:
                    dist.sigma[j] = math.sqrt( self.var[dist_ind,j]  )  # root variance
                    dist.mu[j] =  self.mu[dist_ind,j]


    def updateHyperparameters(self, dists, posterior, data): 
        """
        Updates the hyperparamters in an empirical Bayes fashion as part of the EM parameter estimation.
        
        """
        assert len(dists) == posterior.shape[0]  # XXX debug
        
        # update component-specific hyperparameters
        for i in range(self.nr_comps):
            self.post_sums[i] = numpy.sum(posterior[i,:])
            for j in range(0,self.p):  
                #  var_j = numpy.dot(posterior, (data[:,j] - dist.mu[j])**2 ) / post_sum
                self.var[i,j] = numpy.dot(posterior[i,:], (data[:,j] - dists[i].mu[j])**2 ) / self.post_sums[i]

                if j > 0: # feature 0 is root by convention
                    pid_i_j = dists[i].parents[j]
                    self.cov[i,j] = numpy.dot(posterior[i,:], (data[:,j] - dists[i].mu[j]) * (data[:,pid_i_j] - dists[i].mu[pid_i_j])) / self.post_sums[i] 
                    self.beta[i,j] = self.post_sums[i] / ( (( self.var[i,j] * self.var[i,pid_i_j]) / self.cov[i,j]**2)  -1 )
                    self.nu[i,j] = - self.post_sums[i]  / (2 * dists[i].sigma[j]**2)
           

    def isValid(self,x):
        if not isinstance(x,ConditionalGaussDistribution):
            raise InvalidDistributionInput, "ConditionalGaussPrior: " + str(x)



class ProductDistributionPrior(PriorDistribution):
    """
    Prior for ProductDistribution objects. Basically only holds a list of priors for
    atomar distributions. Necessary for model hierarchy.
    """
    def __init__(self,priorList):
        """
        Constructor
        
        @param priorList: list of PriorDistribution objects
        """
        self.priorList = priorList
        self.dist_nr = len(priorList)

    def __getitem__(self,ind):
        if ind < 0 or ind > self.dist_nr-1:
            raise IndexError, 'Index '+str(ind)
        else:
            return self.priorList[ind]

    def __setitem__(self,ind, item):
        assert isinstance(item,PriorDistribution) 
        
        if ind < 0 or ind > self.dist_nr-1:
            raise IndexError
        else:
            self.priorList[ind] = item
        
    def __eq__(self, other):
        if not isinstance(other,ProductDistributionPrior):
            return False
        if self.dist_nr != other.dist_nr:
            return False
        for i in range(self.dist_nr):
            if not self.priorList[i] == other.priorList[i]:
                return False
        return True        
            
    def pdf(self,dist):
        assert isinstance(dist,ProductDistribution)
        res = 0
        for i in range(self.dist_nr):
            res += self.priorList[i].pdf(dist.distList[i])
        
        return res
        
    def marginal(self,dist,posterior,data):
        assert isinstance(dist,ProductDistribution)
        res = 0
        for i in range(self.dist_nr):
            res += self.priorList[i].marginal(dist.distList[i],posterior,data.getInternalFeature(i))
        
        return res

    def mapMStep(self, dist, posterior, data, mix_pi=None, dist_ind = None):    
        assert dist.suff_dataRange and dist.suff_p, "Attributes for sufficient statistics not initialized."
        assert isinstance(data,DataSet)
        assert isinstance(dist,ProductDistribution)
        assert dist.dist_nr == len(self.priorList)
        
        for i in range(dist.dist_nr):    
            if isinstance(dist.distList[i],MixtureModel):
             # XXX use of isinstance() should be removed, singleFeatureSubset(i) to replace getInternalFeature(i)  ?
                self.priorList[i].mapMStep(dist.distList[i],posterior,data.singleFeatureSubset(i), mix_pi, dist_ind)  
            else:    
                self.priorList[i].mapMStep(dist.distList[i],posterior,data.getInternalFeature(i), mix_pi, dist_ind)  


    def isValid(self,p):
        if not isinstance(p,ProductDistribution):
            raise InvalidDistributionInput,'Not a ProductDistribution.'
        if p.dist_nr != self.dist_nr:
            raise InvalidDistributionInput,'Different dimensions in ProductDistributionPrior and ProductDistribution: ' + str(p.dist_nr)+' != '+ str(self.dist_nr)
        for j in range(p.dist_nr):
            try:
                self[j].isValid(p[j]) 
            except InvalidDistributionInput,ex:
                ex.message += "\n\tin ProductDistributionPrior.priorList["+str(j)+"]"
                raise                       


class ProductDistribution(ProbDistribution):
    """ Class for joined distributions for a vector of random variables with (possibly) different
        types. We assume indepenence between the features.
        Implements the naive Bayes Model.

    """
    def __init__(self, distList):
        """
        Constructor

        @param distList: list of ProbDistribution objects
        """
        # initialize attributes
        self.distList = distList
        self.p = 0    
        self.freeParams = 0
        self.dataRange = []
        self.dist_nr = len(distList)

        # dimension and dataRange for sufficient statistics data
        self.suff_p = 0
        self.suff_dataRange = None
        for dist in distList:
            assert isinstance(dist,ProbDistribution)
            self.p += dist.p 
            self.dataRange.append(self.p)
            self.freeParams += dist.freeParams
            self.suff_p += dist.suff_p

        # initializing dimensions for sufficient statistics data
        self.update_suff_p()
 
    def __eq__(self,other):
        if other.p != self.p or other.dist_nr != self.dist_nr:
            return False
        for i in range(self.dist_nr):
            if not (other.distList[i] == self.distList[i]):
                return False
        return True
    
    def __copy__(self):
        copyList = []
        for i in range(len(self.distList)):
            copyList.append(copy.copy(self.distList[i]))
        
        copy_pd = ProductDistribution(copyList)
        copy_pd.suff_p = self.suff_p
        copy_pd.suff_dataRange = copy.copy(self.suff_dataRange)
        return copy_pd

    def __str__(self):
        outstr = "ProductDist: \n"
        for dist in self.distList:
            outstr += "  " + str(dist) + "\n"
        return outstr

    def __getitem__(self,ind):
        if ind < 0 or ind > self.dist_nr-1:
            raise IndexError
        else:
            return self.distList[ind]    

    def __setitem__(self, ind, value):
        if ind < 0 or ind > self.dist_nr-1:
            raise IndexError
        else:
            self.distList[ind] = value

    def __len__(self):
        return self.dist_nr

    def pdf(self,data):
        assert self.suff_dataRange and self.suff_p, "Attributes for sufficient statistics not initialized."
        if isinstance(data,DataSet):
            res = numpy.zeros(data.N, dtype='Float64')
            for i in range(self.dist_nr):
                if isinstance(self.distList[i], MixtureModel): # XXX only necessary for mixtures of mixtures
                    res += self.distList[i].pdf(data.singleFeatureSubset(i))
                else:
                    res += self.distList[i].pdf(data.getInternalFeature(i))
            return res 
        else:
            raise TypeError, 'DataSet object required, got '+ str(type(data))
        
    def sample(self):
        ls = []
        for i in range(len(self.distList)):    
            s = self.distList[i].sample()
            if type(s) != list:
                ls.append(s)
            else:
                ls += s
        return ls

    def sampleSet(self,nr):
        res = []
        for i in range(nr):
            res.append(self.sample() )
        return res

    def sampleDataSet(self, nr):
        """
        Returns a DataSet object of size 'nr'.
        
        @param nr: size of DataSet to be sampled
        
        @return: DataSet object
        """
        ls = []
        for i in range(nr):
            ls.append(self.sample())
        
        data = DataSet()
        data.dataMatrix= ls
        data.N = nr
        data.p = self.p
        data.sampleIDs = []
        
        for i in range(data.N):
            data.sampleIDs.append("sample"+str(i))    
        
        for h in range(data.p):
            data.headers.append("X_"+str(h))
                  
        data.internalInit(self)
        return data

    def MStep(self,posterior,data,mix_pi=None):

        assert self.suff_dataRange and self.suff_p, "Attributes for sufficient statistics not initialized."
        assert isinstance(data,DataSet),'DataSet required, got '+str(type(data))+'.'
        
        for i in range(self.dist_nr):    
            if isinstance(self.distList[i],MixtureModel):
                self.distList[i].MStep(posterior,data.singleFeatureSubset(i),mix_pi)  
            else:    
                self.distList[i].MStep(posterior,data.getInternalFeature(i))  

    def formatData(self,x):
        res = []
        last_index = 0
        for i in range(len(self.distList)):
            
            # XXX HACK: if distList[i] is an HMM feature there is nothing to be done
            # since all the HMM code was moved to mixtureHMM.py we check whether self.distList[i]
            # is an HMM by string matching  __class__ (for now).
            if self.distList[i].p == 1:
                strg = str(self.distList[i].__class__)
                if strg == 'mixtureHMM.HMM':
                    continue
            
            if self.dist_nr == 1:
                [new_p,dat] = self.distList[i].formatData(x)
                res += dat
            else:    
                if self.distList[i].suff_p == 1:
                    [new_p,dat] = self.distList[i].formatData(x[ self.dataRange[i] -1])
                    res += dat
                else:
                    [new_p,dat] = self.distList[i].formatData(x[last_index:self.dataRange[i]])
                    res += dat
                last_index = self.dataRange[i]        
        return [self.suff_p,res]

    def isValid(self,x):
        last_index = 0
        for i in range(len(self.distList)):
            if self.distList[i].p == 1:
                try:
                    self.distList[i].isValid(x[self.dataRange[i]-1])
                except InvalidDistributionInput,ex:
                    ex.message += "\n\tin ProductDistribution.distList["+str(i)+"]"
                    raise                       
            else:
                try: 
                    self.distList[i].isValid(x[last_index:self.dataRange[i]]) 
                except InvalidDistributionInput,ex:
                    ex.message += "\n\tin ProductDistribution.distList["+str(i)+"]"
                    raise                       
            last_index = self.dataRange[i] 
            
    def flatStr(self,offset):
        offset +=1
        s = "\t"*offset + ";Prod;"+str(self.p)+"\n"
        for d in self.distList:
            s+= d.flatStr(offset)
        return s

    def posteriorTraceback(self,x):
        res = []
        last_index = 0
        assert len(x) == self.suff_p, "Different number of dimensions in data and model."
        
        for i in range(len(self.distList)):
            if self.distList[i].suff_p == 1:
                res += self.distList[i].posteriorTraceback(x[:,self.suff_dataRange[i]-1])
            else:    
                res += self.distList[i].posteriorTraceback(x[:,last_index:self.suff_dataRange[i] ])
                    
            last_index = self.suff_dataRange[i]
        return res 

    def update_suff_p(self):
        old_suff_p = None
        # in case suff_ variables have already been initialized,
        # store old values and compare as consistency check
        if self.suff_dataRange is not None and self.suff_p:
            old_suff_dataRange = self.suff_dataRange
            old_suff_p = self.suff_p
            
        self.suff_dataRange = []
        self.suff_p = 0
        for dist in self.distList:
            self.suff_p += dist.update_suff_p() 
            self.suff_dataRange.append(self.suff_p)
 
        if old_suff_p:
            assert self.suff_p == old_suff_p,str(self.suff_p)+" != "+ str(old_suff_p)
            assert old_suff_dataRange == self.suff_dataRange,str(old_suff_dataRange)+" != "+ str(self.suff_dataRange)
        
        return self.suff_p    
        
#--------------------------------------------------------------------------------------------


class MixtureModelPrior(PriorDistribution):
    """
    Mixture model prior.
    """
    def __init__(self, structPrior, nrCompPrior,piPrior, compPrior):
        """
        Constructor
        
        @param structPrior: hyperparameter over structure complexity (0< structPrior < 1), stored on log scale internally
        @param nrCompPrior:  hyperparameter over number of components (0< nrCompPrior < 1), stored on log scale internally
        @param piPrior: DirichletPrior object
        @param compPrior: list of PriorDistribution objects
        """
        assert isinstance(piPrior,DirichletPrior)
        self.piPrior = piPrior

        for p in compPrior:
            assert isinstance(p,PriorDistribution)
    
        self.compPrior = ProductDistributionPrior(compPrior)  
        self.dist_nr = len(compPrior)
        self.structPrior = numpy.log(structPrior)
        self.nrCompPrior = numpy.log(nrCompPrior)
        
        self.constant_hyperparams = 1
        self.hp_update_indices = []
        for i in range(self.dist_nr):
            if self.compPrior.priorList[i].constant_hyperparams == 0:
                self.hp_update_indices.append(i)
            
        if len(self.hp_update_indices) > 0:
            self.constant_hyperparams = 0  # there is at least one prior which requires hyper parameter updates



    def __str__(self):
        outstr = "MixtureModelPrior: \n"
        outstr += "dist_nr="+str(self.dist_nr)+"\n"
        outstr += "structPrior ="+str(self.structPrior) +"\n"
        outstr += "nrCompPrior ="+str(self.nrCompPrior) +"\n"
        outstr += "  piPrior = "+str(self.piPrior)+"\n"
        for dist in self.compPrior:
            outstr += "    " + str(dist) + "\n"
        
        return outstr


    def __eq__(self,other):
        if not isinstance(other,MixtureModelPrior):
            return False

        if self.structPrior != other.structPrior or self.nrCompPrior != other.nrCompPrior:
            return False
        if not self.piPrior == other.piPrior:
            return False 

        if not self.compPrior == other.compPrior:
            return False    

        return True    
        

    def __copy__(self):
        cp_pi = copy.copy(self.piPrior)
        cp_comp = []
        for i in range(self.dist_nr):
            cp_comp.append( copy.copy(self.compPrior[i] ))
            
        # initialise copy with dummy values for .structPrior, .nrCompPrior   
        cp_pr = MixtureModelPrior( 0.0, 0.0, cp_pi, cp_comp )    
        # set values of hyperparameters
        cp_pr.structPrior = self.structPrior
        cp_pr.nrCompPrior = self.nrCompPrior

        return cp_pr
        

    def pdf(self, mix): 
        #assert isinstance(mix,MixtureModel), str(mix.__class__)
        #assert len(self.compPrior) == mix.components[0].dist_nr
        temp = DiscreteDistribution(mix.G,mix.pi)
        res = self.piPrior.pdf(temp)
        
        # XXX fixed components do not contribute to the prior (HACK)
        # this is needed if we use mixtures of mixtures to model missing data XXX
        if sum(mix.compFix) > 0: 
            for j in range(mix.components[0].dist_nr):
                for l in range(mix.G):
                    if mix.compFix[l] != 2:
                        p = self.compPrior[j].pdf( mix.components[l][j] )
                        res += p
        else:        
            # operate column wise on mix.components
            for j in range(self.dist_nr):

                if not isinstance(mix.components[0].distList[j], MixtureModel):
                    d_j = [mix.components[i].distList[j] for i in range(mix.G)]
                    res += numpy.sum(self.compPrior[j].pdf(d_j))
                else:
                    for i in range(mix.G):
                        res += self.compPrior[j].pdf(mix.components[i][j])
                    
                
        # prior over number of components
        res +=  self.nrCompPrior * mix.G 
        
        # prior over number of distinct groups
        if mix.struct:
            for j in range(mix.components[0].dist_nr):
                res +=  self.structPrior * len(mix.leaders[j])  
        else:
            for j in range(mix.components[0].dist_nr):
                res += self.structPrior * mix.G

        if numpy.isnan(res):
            # uncomment code below for detailed information where the nan value came from (DEBUG)        
#            print '--------------------------------'
#            print 'MixtureModelPrior.pdf ',res
#            temp = DiscreteDistribution(mix.G,mix.pi)
#            print '   MixPrior.pdf.pi:',temp,self.piPrior.pdf(temp)
#            if sum(mix.compFix) > 0: 
#                for j in range(mix.components[0].dist_nr):
#                    for l in range(mix.G):
#                        if mix.compFix[l] != 2:
#                            p = self.compPrior[j].pdf( mix.components[l][j] )
#                            print l,j,p
#            else:        
#                for l in range(mix.G):
#                    #print     mix.components[l]
#                    print '    comp ',l,':',self.compPrior.pdf(mix.components[l])
#                    for o,d in enumerate(self.compPrior):
#                        print '       dist',o,mix.components[l][o],':',d.pdf(mix.components[l][o])
#                        
#            print 'nrCompPrior=',  self.nrCompPrior * mix.G 
#
#            # prior over number of distinct groups
#            if mix.struct:
#                for j in range(mix.components[0].dist_nr):
#                    print '    struct:',self.structPrior * len(mix.leaders[j]) 
#            else:
#                for j in range(mix.components[0].dist_nr):
#                  print '    struct:', self.structPrior * mix.G
            raise ValueError, 'nan result in MixtureModelPrior.pdf'            
        return res
    


    # mapMStep is used for parameter estimation of lower hierarchy mixtures
    def mapMStep(self, dist, posterior, data, mix_pi=None, dist_ind = None):
        dist.mapEM(data, self,1,0.1,silent =True,mix_pi=mix_pi,mix_posterior=posterior)

    def updateHyperparameters(self, dists, posterior, data): 
        
        assert self.constant_hyperparams == 0
        assert isinstance(dists, MixtureModel) # XXX debug
        
        for j in self.hp_update_indices:
            d_j = [dists.components[i].distList[j] for i in range(dists.G)]
            self.compPrior[j].updateHyperparameters(d_j,posterior,data.getInternalFeature(j) )



    def flatStr(self,offset):
        offset +=1
        s = "\t"*offset + ";MixPrior;"+str(self.dist_nr)+";"+str(numpy.exp(self.structPrior))+";"+str(numpy.exp(self.nrCompPrior))+"\n"

        s+= self.piPrior.flatStr(offset)
        for d in self.compPrior:
            s+= d.flatStr(offset)
        
        return s
             
    def posterior(self, dist):
        raise NotImplementedError
  
    def isValid(self, m):
        if not isinstance(m,MixtureModel):
            raise InvalidDistributionInput, "MixtureModelPrior: " + str(m)
        else:
            if self.piPrior.M != m.G:
                raise InvalidDistributionInput, "MixtureModelPrior: invalid size of piPrior."
            
            try:
                # check validity of each component
                for i in range(m.G):
                    self.compPrior.isValid(m.components[i]) 
            except InvalidDistributionInput,ex:
                ex.message += "\n\tin MixtureModelPrior for component "+str(i)
                raise                       

    def structPriorHeuristic(self, delta, N):
        """
        Heuristic for setting the structure prior hyper-parameter 'self.structPrior', depending
        on the size of a data set 'N' and parameter 'delta'.
        """
        self.structPrior =  - numpy.log(1+delta)*N


    def mapMStepMerge(self, group_list):
        new_dist = copy.copy(group_list[0].dist)
        new_req_stat = copy.deepcopy(group_list[0].req_stat)
        
        assert new_dist.dist_nr == 1 # XXX
        assert new_dist.G == 2

        for i in range(new_dist.G):
            if new_dist.compFix[i] == 2:
                continue
            sub_group_list = []                
            for r in range(len(group_list)):
                sub_group_list.append( CandidateGroup( group_list[r].dist.components[i][0], group_list[r].post_sum, group_list[r].pi_sum,  group_list[r].req_stat[i] ) )

            d_i = self.compPrior[0].mapMStepMerge(sub_group_list)
            new_dist.components[i][0] = d_i.dist
            new_req_stat[i] = d_i.req_stat

        return CandidateGroup(new_dist, d_i.post_sum, d_i.pi_sum, new_req_stat )




class MixtureModel(ProbDistribution):
    """
    Class for a context-specific independence (CSI) mixture models.
    The components are naive Bayes models (i.e. ProductDistribution objects). 
    """
    def __init__(self,G, pi, components,compFix=None,struct=0, identifiable = 1):
        """
        Constructor

        @param G: number of components
        @param pi: mixture weights
        @param components: list of ProductDistribution objects, each entry is one component
        @param compFix: list of optional flags for fixing components in the reestimation
                         the following values are supported: 1 distribution parameters are fixed, 2 distribution 
                         parameters and mixture coefficients are fixed
        @param struct: Flag for CSI structure, 0 = no CSI structure, 1 = CSI structure
        """
        assert len(pi) == len(components) == G, str(len(pi)) +', '+ str(len(components))+', '+str(G)
        assert abs((1.0 - sum(pi))) < 1e-12, "sum(pi) = " + str(sum(pi)) +", "+str(abs((1.0 - sum(pi))))

        self.freeParams = 0
        self.p = components[0].p

        # Internally components must be a list of ProductDistribution objects. In case the input is a list of
        # ProbDistributions we convert components accordingly.
        if isinstance(components[0],ProbDistribution) and not isinstance(components[0],ProductDistribution):
            # make sure all elements of components are ProbDistribution of the same dimension
            for c in components:
                assert isinstance(c,ProbDistribution)
                assert c.p == self.p
            for i,c in enumerate(components)                :
                components[i] = ProductDistribution([c])

        self.dist_nr = components[0].dist_nr                

        # some checks to ensure model validity
        for c in components:
            # components have to be either ProductDistribution objects or a list of univariate ProbDistribution objects
            assert isinstance(c,ProductDistribution), "Got "+str(c.__class__)+" as component."
            assert self.p == c.p,str(self.p)+" != "+str(c.p)
            assert self.dist_nr == c.dist_nr
            self.freeParams += c.freeParams


        self.freeParams += G-1  # free parameters of the mixture coefficients

        self.G = G  # Number of components
        self.pi = numpy.array(pi,dtype='Float64')  # vector of mixture weights

        # XXX  check numpy capabilities for arrays of objects 
        self.components = components   # list of distribution objects 

        self.suff_p = None  # dimension of sufficient statistic data

        self.nr_tilt_steps = 10 # number of steps for deterministic annealing in the EM, 10 is default
        self.heat = 0.5  # initial heat parameter for deterministic annealing
        
        # compFix contains flags for each component, which determine whether the distribution parameter in a
        # component will be skipped in the reestimation, 
        if compFix:
            assert len(compFix) == self.G, str(len(compFix)) +" != "+ str(self.G)
            self.compFix = compFix
        else:
            self.compFix = [0] * self.G    

        # initializing dimensions for sufficient statistics data
        self.update_suff_p()

        self.struct = struct        
        if self.struct:
            self.initStructure()
            
        # flag that determines whether identifiability is enforced in training
        self.identFlag = identifiable    
        
        # error tolerance for the objective function in the parameter training
        self.err_tol = 1e-6
        
        # minmal mixture coefficient value
        self.min_pi = 0.05 # TEST

    def __eq__(self,other):
        res = False
        if isinstance(other,MixtureModel):
            if  numpy.allclose(other.pi,self.pi) and other.G == self.G:
                res = True
                for i in range(self.G):
                    if not (other.components[i] == self.components[i]):
                        #print other.components[i] ,"!=", self.components[i]
                        return False

        return res        

    def __copy__(self):
        copy_components = []
        copy_pi = copy.deepcopy(self.pi)  
        copy_compFix = copy.deepcopy(self.compFix)
        for i in range(self.G):
            copy_components.append(copy.deepcopy(self.components[i]))
        
        copy_model = MixtureModel(self.G, copy_pi, copy_components,compFix = copy_compFix)
        copy_model.nr_tilt_steps = self.nr_tilt_steps    
        copy_model.suff_p = self.suff_p
        copy.identFlag = self.identFlag
        
        if self.struct:
            copy_model.initStructure()
            
            copy_leaders = copy.deepcopy(self.leaders)
            copy_groups = copy.deepcopy(self.groups)
        
            copy_model.leaders = copy_leaders
            copy_model.groups = copy_groups
        
        return copy_model
        
    def __str__(self):
        s = "G = " + str(self.G)
        s += "\np = " + str(self.p)
        s += "\npi =" + str(self.pi)+"\n"
        s += "compFix = " + str(self.compFix)+"\n"
        for i in range(self.G):
            s += "Component " + str(i)+ ":\n"
            s += "  "+str(self.components[i])+"\n"
        
        if self.struct:
            s+= "\nCSI structure:\n"
            s += "leaders:"+str(self.leaders)+"\n"
            s += "groups:"+str(self.groups)+"\n"
        
        return s

    def initStructure(self):
        """
        Initializes the CSI structure.
        """
        
        self.struct = 1
        # for a model with group structure we have to check for valid model topology
        if self.struct:
            # ensuring identical number of distributions in all components
            nr = self.components[0].dist_nr
            for i in range(1,self.G):
                assert isinstance(self.components[i],ProductDistribution)
                assert self.components[i].dist_nr == nr
                # checking for consistent dimensionality of elementar distributions among components
                for j in range(nr):
                    assert self.components[i][j].p == self.components[0][j].p
                    assert self.components[i][j].freeParams == self.components[0][j].freeParams

            # if there is already a CSI structure in the model, components within the same group
            # share a reference to the same distribution object. For the new structure we need to make copies.
            if hasattr(self,'groups'):
                for j in range(nr):
                    for i in self.groups[j]:
                        for r in self.groups[j][i]:
                            self.components[r][j] = copy.copy(self.components[r][j])
                            
            # Variables for model structure.
            # leaders holds for each dimension the indixes of the group representing components. Initially each
            # component is the representative of itself, e.g. there are no groups.
            self.leaders = []
        
            # groups is a list of dictionaries, one for each dimension p. The group members of a leader are hashed with the
            # leaders index as key. Initially all components are inserted with an empty list.
            self.groups = []

            for i in range(nr):
               self.leaders.append(range(self.G))
               d= {}
               for j in range(self.G):
                    d[j] = []
               self.groups.append(d)

    def modelInitialization(self,data,rtype=1, missing_value = None):
        """
        Perform model initialization given a random assigment of the
        data to the components.
            
        @param data: DataSet object
        @param rtype: type of random assignments.
        0 = fuzzy assingment
        1 = hard assingment
        @param missing_value: missing symbol to be ignored in parameter estimation (if applicable)
        
        @return: posterior assigments
        """
        if not isinstance(data,DataSet):
            raise TypeError, "DataSet object required, got"+ str(data.__class__)
        else:
            if data.internalData is None:
                data.internalInit(self)
        
        # reset structure if applicable
        if self.struct:
            self.initStructure()
            
        l = numpy.zeros((self.G,len(data)),dtype='Float64')
        for i in range(len(data)):
            if rtype == 0:
              for j in range(self.G):
                 l[j,i] = random.uniform(0.1,1)
              s = sum(l[:,i])
              for j in range(self.G):
                  l[j,i] /= s
            else:
               l[random.randint(0,self.G-1),i] = 1 
       
       
        # do one M Step
        fix_pi = 1.0
        unfix_pi = 0.0
        fix_flag = 0   # flag for fixed mixture components
        for i in range(self.G):
            # setting values for pi
            if self.compFix[i] == 2:
                fix_pi -= self.pi[i]
                fix_flag = 1
            else:
                self.pi[i] =  l[i,:].sum() / len(data)                
                unfix_pi += self.pi[i]

            if self.compFix[i] == 1 or self.compFix[i] == 2:
                # fixed component
                continue
            else:    
                # components are product distributions that may contain mixtures
                last_index = 0

                for j in range(self.components[i].dist_nr):    
                    if isinstance(self.components[i][j],MixtureModel):
                        dat_j = data.singleFeatureSubset(j)
                        self.components[i][j].modelInitialization(dat_j,rtype=rtype,missing_value=missing_value)      
                    else:    
                        loc_l = l[i,:]
                        # masking missing values from parameter estimation
                        if data.missingSymbols.has_key(j):
                            ind_miss = data.getMissingIndices(j)
                            for k in ind_miss:
                                loc_l[k] = 0.0
                        
                        self.components[i][j].MStep(loc_l,data.getInternalFeature(j))      

        # renormalizing mixing proportions in case of fixed components
        if fix_flag:
            if unfix_pi == 0.0:
                #print "----\n",self,"----\n"
                 print "unfix_pi = ", unfix_pi
                 print "fix_pi = ", fix_pi
                 print "pi = ", self.pi
                 raise RuntimeError, "unfix_pi = 0.0"  
                
            for i in range(self.G):
                if self.compFix[i] == 0:
                    self.pi[i] = (self.pi[i] * fix_pi) / unfix_pi
        return l


    def pdf(self,x):
        logp_list = numpy.zeros((self.G,len(x)),dtype='Float64')
        for i in range(self.G):
            if self.pi[i] == 0.0:
                log_pi = float('-inf')
            else:
                log_pi =numpy.log(self.pi[i]) 
            logp_list[i] = log_pi + self.components[i].pdf(x)

        p = numpy.zeros(len(x),dtype='Float64')
        for j in range(len(x)):
            p[j] = sumlogs(logp_list[:,j])
        return p    

    def sample(self):
        sum = 0.0
        p = random.random()
        for k in range(self.G):
            sum += self.pi[k]
            if sum >= p:
                break
        return self.components[k].sample()
    
    def sampleSet(self, nr):
        ls = []
        for i in range(nr): 
            sum = 0.0
            p = random.random()
            for k in range(self.G):
                sum += self.pi[k]
                if sum >= p:
                    break
            ls.append(self.components[k].sample())
        return ls


    def sampleDataSet(self, nr):
        """
        Returns a DataSet object of size 'nr'.
        
        @param nr: size of DataSet to be sampled
        
        @return: DataSet object
        """
        ls = self.sampleSet(nr)
        data = DataSet()
        data.dataMatrix= ls
        data.N = nr
        data.p = self.p
        data.sampleIDs = []
        
        for i in range(data.N):
            data.sampleIDs.append("sample"+str(i))    
        
        for h in range(data.p):
            data.headers.append("X_"+str(h))
                  
        data.internalInit(self)
        return data

    def sampleDataSetLabels(self, nr):
        """
        Samples a DataSet of size 'nr' and returns the DataSet and the true
        component labels
        
        @param nr: size of DataSet to be sampled
        
        @return: tuple of DataSet object and list of labels
        """
        [c,ls] = self.sampleSetLabels(nr)
        
        data = DataSet()
        data.dataMatrix= ls
        data.N = nr
        data.p = self.p
        data.sampleIDs = []
        
        for i in range(data.N):
            data.sampleIDs.append("sample"+str(i))    

        for h in range(data.p):
            data.headers.append("X_"+str(h))
                  
        data.internalInit(self)
        return [data,c]
        
    def sampleSetLabels(self, nr):
        """ Same as sample but the component labels are returned as well. 
            Useful for testing purposes mostly.
        
        """
        ls = []
        label = []
        for i in range(nr): 
            sum = 0.0
            p = random.random()
            for k in range(self.G):
                sum += self.pi[k]
                if sum >= p:
                    break
            label.append(k)
            ls.append(self.components[k].sample())
        
        return [numpy.array(label), ls]

    def EM(self, data, max_iter, delta,silent = False, mix_pi=None, mix_posterior= None, tilt = 0, EStep=None, EStepParam = None ):
        """
        Reestimation of mixture parameters using the EM algorithm.
        
        @param data: DataSet object
        @param max_iter: maximum number of iterations
        @param delta: minimal difference in likelihood between two iterations before
        convergence is assumed.
        @param silent: 0/1 flag, toggles verbose output
        @param mix_pi: [internal use only] necessary for the reestimation of
        mixtures as components
        @param mix_posterior:[internal use only] necessary for the reestimation of
        mixtures as components
        @param tilt: 0/1 flag, toggles the use of a deterministic annealing in the training
        @param EStep: function implementing the EStep, by default self.EStep
        @param EStepParam: additional paramenters for more complex EStep implementations
        
        @return: tuple of posterior matrix and log-likelihood from the last iteration
        """
        if isinstance(data, numpy.ndarray):
            raise TypeError, "DataSet object required."
        elif  isinstance(data, DataSet):
            if data.internalData is None:
                if not silent:
                    sys.stdout.write("Parsing data set...")
                    sys.stdout.flush()
                data.internalInit(self)
                if not silent:
                    sys.stdout.write("done\n")
                    sys.stdout.flush()
        else:    
            raise ValueError, "Unknown input type format: " + str(data.__class__)

        log_p_old = -1.0
        step = 0

        if EStep == None:
            EStep = self.EStep

        # if deterministic annealing is activated, increase number of steps by self.nr_tilt_steps 
        if tilt:
            if not silent:
                sys.stdout.write("Running EM with "+ str(self.nr_tilt_steps) +" steps of deterministic annealing.\n" )
            max_iter += self.nr_tilt_steps

        log_p = 0.0            
        while 1:
            [log_l,log_p] = EStep(data,mix_posterior,mix_pi,EStepParam)
            if log_p_old != -1.0 and not silent and step != 0:
                if tilt and step <= self.nr_tilt_steps:
                    sys.stdout.write("TILT Step "+str(step)+": log likelihood: "+str(log_p)+"\n")
                else:
                    sys.stdout.write("Step "+str(step)+": log likelihood: "+str(log_p_old)+"   (diff="+str(diff)+")\n")

            # checking for convergence 
            diff = (log_p - log_p_old)
            
            if diff < 0.0 and step > 1 and abs(diff / log_p_old) > self.err_tol:
                print log_p,log_p_old, diff,step,abs(diff / log_p_old)
                print "WARNING: EM divergent."
                raise ConvergenceFailureEM,"Convergence failed."

            if numpy.isnan(log_p):
                print "WARNING: One sample was not assigned."
                raise ConvergenceFailureEM,"Non assigned element."               
                
            if (not tilt or (tilt and step+1 >= self.nr_tilt_steps)) and diff >= 0.0 and delta >= abs(diff) and max_iter != 1:
                if not silent:
                    sys.stdout.write("Step "+str(step)+": log likelihood: "+str(log_p)+"   (diff="+str(diff)+")\n")
                    sys.stdout.write("Convergence reached with log_p "+ str(log_p)+ " after "+str(step)+" steps.\n")
                if self.identFlag:
                    self.identifiable()
                return (log_l,log_p)           

            if step == max_iter:
                if not silent:
                    sys.stdout.write("Max_iter "+str(max_iter)+" reached -> stopping\n")
                if self.identFlag:
                    self.identifiable()
                return (log_l,log_p)  

            log_p_old = log_p            

            # compute posterior likelihood matrix from log posterior
            l = numpy.exp(log_l)

            # deterministic annealing, shifting posterior toward uniform distribution.
            if tilt and step+1 <= self.nr_tilt_steps and mix_posterior is None:
                h = self.heat - (step * (self.heat/(self.nr_tilt_steps))  )
                for j in range(data.N):
                    uni = 1.0 / self.G
                    tilt_l = (uni - l[:,j]) * h
                    
                    l[:,j] += tilt_l
                    #print l[:,j]

            # variables for component fixing
            fix_pi = 1.0
            unfix_pi = 0.0
            fix_flag = 0   # flag for fixed mixture components
            
            # update component parameters and mixture weights
            for i in range(self.G):
                if self.compFix[i] == 2:   # pi[i] is fixed
                    fix_pi -= self.pi[i]
                    fix_flag = 1
                    continue
                    
                else:
                    # for mixtures of mixtures we need to multiply in the mix_pi[i]s                 
                    if mix_pi is not None:
                        self.pi[i] =  ( l[i,:].sum() / (data.N * mix_pi))  
                    else:
                        self.pi[i] =  ( l[i,:].sum() / (data.N ) )
                    
                    unfix_pi += self.pi[i]

                if self.compFix[i] == 1 or self.compFix[i] == 2:
                    #print "  Component ",i," is skipped from reestimation."
                    continue                
                else:
                    # Check for model structure
                    if not self.struct:
                        # there might be mixtures down in the hierarchy, so pi[i] is passed to MStep
                        self.components[i].MStep(l[i,:], data, self.pi[i])
                    
            # if there is a model structure we update the leader distributions only 
            if self.struct:
                datRange = self.components[0].suff_dataRange
                        
                for j in range(self.dist_nr):
                    for k in self.leaders[j]:
                        if j == 0:
                            prev = 0
                        else:
                            prev = datRange[j-1]     
                        
                        # compute group posterior
                        g_post = numpy.array(l[k,:].tolist(),dtype='Float64')

                        for memb in self.groups[j][k]:
                            g_post += l[memb,:]
                        
                        if isinstance(self.components[k][j],MixtureModel):
                            self.components[k][j].MStep(g_post,data.singleFeatureSubset(j),self.pi[k] )
                        else:    
                            self.components[k][j].MStep(g_post,data.getInternalFeature(j))

            # renormalizing mixing proportions in case of fixed components
            if fix_flag:
                if unfix_pi == 0.0:
                    #print "----\n",self,"----\n"
                    print "unfix_pi = ", unfix_pi
                    print "fix_pi = ", fix_pi
                    print "pi = ", self.pi
                    raise RuntimeError, "unfix_pi = 0.0"  
                
                for i in range(self.G):
                    if self.compFix[i] == 0:
                        self.pi[i] = (self.pi[i] * fix_pi) / unfix_pi
            
            sys.stdout.flush()
            step += 1              

    # XXX all old code should be collected together (-> remove ?)
    def EStep_old(self,data,mix_posterior=None,mix_pi=None,EStepParam=None):
        """ [Old implementation, kept around for regression testing]
        
        Reestimation of mixture parameters using the EM algorithm.
                
        @param data: DataSet object
        @param mix_posterior:[internal use only] necessary for the reestimation of
        mixtures as components
        @param mix_pi: [internal use only] necessary for the reestimation of
        mixtures as components
        @param EStepParam: additional paramenters for more complex EStep implementations, in
        this implementaion it is ignored

        @return: tuple of log likelihood matrices and sum of log-likelihood of components
               
        """ 
        log_l = numpy.zeros((self.G,data.N),dtype='Float64')
        log_col_sum = numpy.zeros(data.N,dtype='Float64')  # array of column sums of log_l
        log_pi = numpy.log(self.pi)  # array of log mixture coefficients

        # compute log of mix_posterior (if present)
        if mix_posterior is not None:
            log_mix_posterior = numpy.log(mix_posterior)

        # computing log posterior distribution 
        for i in range(self.G):            
            #print i,self.components[i].pdf(data).tolist()
            
            # XXX cache redundant pdfs for models with CSI structure 
            log_l[i] = log_pi[i] + self.components[i].pdf(data)

        for j in range(data.N):
            log_col_sum[j] = sumlogs(log_l[:,j]) # sum over jth column of log_l    

            # if posterior is invalid, check for model validity
            if log_col_sum[j] == float('-inf'): 

                # if self is at the top of hierarchy, the model is unable to produce the
                # sequence and an exception is raised.
                if mix_posterior is None and not mix_pi:
                    #print "\n---- Invalid -----\n",self,"\n----------"
                    #print "\n---------- Invalid ---------------"
                    #print "mix_pi = ", mix_pi
                    #print "x[",j,"] = ", data.getInternalFeature(j)
                    #print "l[:,",j,"] = ", log_l[:,j] 
                    #print 'data[',j,'] = ',data.dataMatrix[j]
                    
                    raise InvalidPosteriorDistribution, "Invalid posterior distribution."

            # for valid posterior, normalize and go on    
            else:
                # normalizing log posterior
                log_l[:,j] = log_l[:,j] - log_col_sum[j]
                # adjusting posterior for lower hierarchy mixtures
                if mix_posterior is not None:
                    # multiplying in the posterior of upper hierarchy mixture
                    log_l[:,j] = log_l[:,j] + log_mix_posterior[j]
      
        # computing data log likelihood as criteria of convergence
        log_p = numpy.sum(log_col_sum)
        
        return log_l, log_p


    def EStep(self,data,mix_posterior=None,mix_pi=None,EStepParam=None):
        """Reestimation of mixture parameters using the EM algorithm.
                
        @param data: DataSet object
        @param mix_posterior:[internal use only] necessary for the reestimation of
        mixtures as components
        @param mix_pi: [internal use only] necessary for the reestimation of
        mixtures as components
        @param EStepParam: additional paramenters for more complex EStep implementations, in
        this implementaion it is ignored

        @return: tuple of log likelihood matrices and sum of log-likelihood of components
               
        """ 
        log_l = numpy.zeros((self.G,data.N),dtype='Float64')
        log_col_sum = numpy.zeros(data.N,dtype='Float64')  # array of column sums of log_l
        log_pi = numpy.log(self.pi)  # array of log mixture coefficients

        # compute log of mix_posterior (if present)
        if mix_posterior is not None:
            log_mix_posterior = numpy.log(mix_posterior)

        # computing log posterior distribution 
        for i in range(self.G):            
            #print i,self.components[i].pdf(data).tolist()
            
            # XXX cache redundant pdfs for models with CSI structure 
            log_l[i] = log_pi[i] + self.components[i].pdf(data)

        # log_l is normalized in-place 
        log_p = _C_mixextend.get_normalized_posterior_matrix(log_l)

        if log_p == float('-inf'):
            raise InvalidPosteriorDistribution, "Invalid posterior distribution."

        if mix_posterior is not None:
            # multiplying in the posterior of upper hierarchy mixture
            log_l = log_l + log_mix_posterior
           
        return log_l, log_p



    def randMaxEM(self,data,nr_runs,nr_steps,delta,tilt=0,silent=False):
        """
        Performs `nr_runs` normal EM runs with random initial parameters
        and returns the model which yields the maximum likelihood.
        
        @param data: DataSet object
        @param nr_runs: number of repeated random initializations
        @param nr_steps: maximum number of steps in each run
        @param delta: minimim difference in log-likelihood before convergence
        @param tilt: 0/1 flag, toggles the use of a deterministic annealing in the training
        @param silent:0/1 flag, toggles verbose output
        
        @return: log-likelihood of winning model
        """
        if isinstance(data, numpy.ndarray):
            raise TypeError, "DataSet object required."
        elif  isinstance(data, DataSet):
            if data.internalData is None:
                if not silent:
                    sys.stdout.write("Parsing data set...")
                    sys.stdout.flush()
                data.internalInit(self)
                if not silent:
                    sys.stdout.write("done\n")
                    sys.stdout.flush()
        else:    
            raise ValueError, "Unknown input type format: " + str(data.__class__)

        logp_list = []
        best_logp = float('-inf')
        best_model = None
        candidate_model = copy.copy(self)  # copying the model parameters 

        for i in range(nr_runs):
            # we do repeated random intializations until we get a model with valid posteriors in all components
            init = 0
            while not init:
                try:        
                    candidate_model.modelInitialization(data)    # randomizing parameters of the model copy
                except InvalidPosteriorDistribution:
                    pass
                else:
                    init = 1     
            
            try:        
                (l,log_p) = candidate_model.EM(data,nr_steps,delta,silent=silent,tilt=tilt)  # running EM 
            except ConvergenceFailureEM:
                sys.stdout.write("Run "+str(i)+" produced invalid model, omitted.\n")
            except InvalidPosteriorDistribution:
                sys.stdout.write("Run "+str(i)+" produced invalid model, omitted.\n")
            else:
                logp_list.append(log_p)

                # check whether current model is better than previous optimum
                if log_p > best_logp:
                    best_model = copy.copy(candidate_model)
                    best_logp = log_p
      
        if not silent:
            print "\nBest model likelihood over ",nr_runs,"random initializations:"
            print "Model likelihoods:",logp_list
            print "Average logp: ", sum(logp_list)/float(nr_runs)," SD:",numpy.array(logp_list).std()
            print "Best logp:",best_logp

        # check whether at least one run was sucessfully completed        
        if best_model == None:
            raise ConvergenceFailureEM, 'All '+ str(nr_runs)+' runs have failed.'
        
        self.components = best_model.components  # assign best parameter set to model 'self'
        self.pi = best_model.pi

        return best_logp  # return final data log likelihood


    def structureEM(self,data,nr_repeats,nr_runs,nr_steps,delta,tilt=0,silent=False):
        """ 
        EM training for models with CSI structure.
        First a candidate model is generated by using the randMaxEM procedure,
        then the structure is trained. 
            
        @param data: DataSet object
        @param nr_repeats: number of candidate models to be generated
        @param nr_runs: number of repeated random initializations
        @param nr_steps: maximum number of steps for the long training run
        @param delta: minimim difference in log-likelihood before convergence
        @param tilt: 0/1 flag, toggles the use of deterministic annealing in the training
        @param silent:0/1 flag, toggles verbose output
        
        @return: log-likelihood of winning model
        """
        if isinstance(data, numpy.ndarray):
            raise TypeError, "DataSet object required."
        elif  isinstance(data, DataSet):
            if data.internalData is None:
                sys.stdout.write("Parsing data set...")
                sys.stdout.flush()
                data.internalInit(self)
                sys.stdout.write("done\n")
                sys.stdout.flush()
        else:    
            raise ValueError, "Unknown input type format: " + str(data.__class__)
    
        assert self.struct 
        best_logp = None
        best_candidate = None 
        candidate_model = copy.copy(self)  # copying the model parameters         
        for r in range(nr_repeats):
            error = 0
            candidate_model.modelInitialization(data)    # randomizing parameters of the model copy
            log_p = candidate_model.randMaxEM(data,nr_runs,nr_steps,delta,tilt=tilt,silent=silent)
            ch = candidate_model.updateStructureGlobal(data,silent= silent) 
            if not silent:
                print "Changes = ",ch
            while ch != 0:
                try:
                    candidate_model.EM(data,30,0.01,silent=1,tilt=0)
                    ch = candidate_model.updateStructureGlobal(data,silent= silent)
                    if not silent:
                        print "Changes = ",ch
                except ConvergenceFailureEM:
                    error = 1
                    break
            
            if not error:
                (l,log_p) = candidate_model.EM(data,30,0.01,silent=1,tilt=0)
                if r == 0 or log_p > best_logp:
                    best_logp = log_p
                    best_candidate = copy.copy(candidate_model)
            else:
                continue

        self.components = best_candidate.components  # assign best parameter set to model 'self'
        self.pi = best_candidate.pi
        self.groups = best_candidate.groups
        self.leaders = best_candidate.leaders
        self.freeParams = best_candidate.freeParams
        return best_logp

    # MStep is used for parameter estimation of lower hierarchy mixtures
    def MStep(self,posterior,data,mix_pi=None):       
        self.EM(data,1,0.1,silent =True,mix_pi=mix_pi,mix_posterior=posterior)

    def mapEM(self,data, prior, max_iter, delta,silent = False, mix_pi=None, mix_posterior= None, tilt = 0):
        """
        Reestimation of maximum a posteriori (MAP) mixture parameters using the EM algorithm.
        
        @param data: DataSet object
        @param max_iter: maximum number of iterations
        @param prior: an appropriate MixtureModelPrior object
        @param delta: minimal difference in likelihood between two iterations before
        convergence is assumed.
        @param silent: 0/1 flag, toggles verbose output
        @param mix_pi: [internal use only] necessary for the reestimation of
        mixtures as components
        @param mix_posterior:[internal use only] necessary for the reestimation of
        mixtures as components
        @param tilt: 0/1 flag, toggles the use of a deterministic annealing in the training
        
        @return: tuple of posterior matrix and log-likelihood from the last iteration
        """

        if isinstance(data, numpy.ndarray):
            raise TypeError, "DataSet object required."
        elif  isinstance(data, DataSet):
            if data.internalData is None:
                if not silent:
                    sys.stdout.write("Parsing data set...")
                    sys.stdout.flush()
                data.internalInit(self) 
                if not silent:
                    sys.stdout.write("done\n")
                    sys.stdout.flush()
        else:    
            raise ValueError, "Unknown input type format: " + str(data.__class__)

        log_p_old = float('-inf')
        step = 0

        # if deterministic annealing is activated, increase number of steps by self.nr_tilt_steps 
        if tilt:
            if not silent:
                sys.stdout.write("Running EM with "+ str(self.nr_tilt_steps) +" steps of deterministic annealing.\n" )
            max_iter += self.nr_tilt_steps
            
        # for lower hierarchy mixture we need the log of mix_posterior
        if mix_posterior is not None:
            log_mix_posterior = numpy.log(mix_posterior)

        while 1:
            log_p = 0.0
            # matrix of log posterior probs: components# * (sequence positions)
            log_l = numpy.zeros((self.G,data.N),dtype='Float64')
            #log_col_sum = numpy.zeros(data.N,dtype='Float64')  # array of column sums of log_l
            log_pi = numpy.log(self.pi)  # array of log mixture coefficients
            
            # computing log posterior distribution 
            for i in range(self.G):            
                log_l[i] = log_pi[i] + self.components[i].pdf(data) 
            
            # computing data log likelihood as criteria of convergence
            # log_l is normalized in-place and likelihood is returned as log_p
            log_p = _C_mixextend.get_normalized_posterior_matrix(log_l)

            # adjusting posterior for lower hierarchy mixtures
            if mix_posterior is not None:
                # multiplying in the posterior of upper hierarchy mixture
                log_l = log_l + log_mix_posterior

            # compute posterior likelihood matrix from log posterior
            l = numpy.exp(log_l)

            # update prior hyper parametes in an empirical Bayes fashion, if appropriate
            if prior.constant_hyperparams == 0:
                prior.updateHyperparameters(self, l, data)

            
            # we have to take the parameter prior into account to form the objective function
            # Since we assume independence between parameters in different components, the prior
            # contribution is given by a product over the individual component and structure priors
            log_prior = prior.pdf(self)

            # calculate objective function
            log_p += log_prior

            # checking for convergence 
            # XXX special case for the parameter udpate of lower hierarchy mixtures
            if max_iter == 1 and mix_posterior != None and log_p == float('-inf'):
                diff = -1.0  # dummy value for diff 
            else:
                diff = (log_p - log_p_old) 

            if log_p_old != -1.0 and not silent and step > 0:
                if tilt and step <= self.nr_tilt_steps:
                    sys.stdout.write("TILT Step "+str(step)+": log posterior: "+str(log_p)+"\n")
                else:
                    sys.stdout.write("Step "+str(step)+": log posterior: "+str(log_p)+ "   (diff="+str(diff)+")\n")

            if diff < 0.0 and step > 1 and abs(diff / log_p_old) > self.err_tol:
                #print log_p,log_p_old, diff,step,abs(diff / log_p_old)
                #print "WARNING: EM divergent."
                raise ConvergenceFailureEM,"Convergence failed, EM divergent: "  
                                
            if (not tilt or (tilt and step+1 >= self.nr_tilt_steps)) and delta >= diff and max_iter != 1: 
                if not silent:
                    sys.stdout.write("Convergence reached with log_p "+ str(log_p)+ " after "+str(step)+" steps.\n")
                if self.identFlag:
                    self.identifiable()
                return (log_l,log_p)  

            log_p_old = log_p
            if step == max_iter:
                if not silent:
                    sys.stdout.write("Max_iter "+str(max_iter)+" reached -> stopping\n")

                if self.identFlag:
                    self.identifiable()
                return (log_l,log_p) 

            # deterministic annealing, shifting posterior toward uniform distribution.
            if tilt and step+1 <= self.nr_tilt_steps and mix_posterior is None:
                h = self.heat - (step * (self.heat/(self.nr_tilt_steps))  )
                for j in range(data.N):
                    uni = 1.0 / self.G
                    tilt_l = (uni - l[:,j]) * h
                    l[:,j] += tilt_l

            # variables for component fixing
            fix_pi = 1.0
            unfix_pi = 0.0
            fix_flag = 0   # flag for fixed mixture components
            
            # update component parameters and mixture weights
            for i in range(self.G):
                if self.compFix[i] & 2:   # pi[i] is fixed
                    fix_pi -= self.pi[i]
                    fix_flag = 1
                    continue
                else:
                    # for mixtures of mixtures we need to multiply in the mix_pi[i]s                 
                    if mix_pi is not None:                                                  
                        self.pi[i] =  ( l[i,:].sum() + prior.piPrior.alpha[i] -1.0 ) / ((data.N * mix_pi) + prior.piPrior.alpha_sum - self.G ) 
                        #print i, ( l[i,:].sum() + prior.piPrior.alpha[i] -1.0 ),((data.N * mix_pi) + prior.piPrior.alpha_sum - self.G )
                    else:
                        self.pi[i] =  ( l[i,:].sum() + prior.piPrior.alpha[i] -1.0 ) / (data.N + ( prior.piPrior.alpha_sum - self.G) ) 
                    
                    unfix_pi += self.pi[i]

                if self.compFix[i] & 1:
                    continue                
                else:
                    # Check for model structure
                    if not self.struct:
                        prior.compPrior.mapMStep(self.components[i],l[i,:],data,self.pi[i],i)
                    
            # if there is a model structure we update the leader distributions only 
            if self.struct:
                for j in range(self.dist_nr):
                    for k in self.leaders[j]:
                        # compute group posterior
                        # XXX extension function for pooled posterior ?
                        g_post = copy.deepcopy(l[k,:])  
                        g_pi = self.pi[k]
                        for memb in self.groups[j][k]:
                            g_post += l[memb,:]
                            g_pi += self.pi[memb]
                        
                        if isinstance(self.components[k][j],MixtureModel): 
                            prior.compPrior[j].mapMStep(self.components[k][j],g_post,data.singleFeatureSubset(j), g_pi,k)
                        else:    
                            try:
                                prior.compPrior[j].mapMStep(self.components[k][j],g_post,data.getInternalFeature(j),g_pi,k )
                            except InvalidPosteriorDistribution:
                                raise    
            
            # renormalizing mixing proportions in case of fixed components
            if fix_flag:
                if unfix_pi == 0.0:
                    #print "----\n",self,"----\n"
                    #print "unfix_pi = ", unfix_pi
                    #print "fix_pi = ", fix_pi
                    #print "pi = ", self.pi
                    #print self
                    raise ValueError, "unfix_pi = 0.0"  
                for i in range(self.G):
                    if self.compFix[i] == 0:
                        self.pi[i] = (self.pi[i] * fix_pi) / unfix_pi

            sys.stdout.flush()
            step += 1               

    def classify(self, data, labels = None, entropy_cutoff = None, silent = 0, EStep=None, EStepParam = None ):
        """
        Classification of input 'data'.
        Assignment to mixture components by maximum likelihood over
        the component membership posterior. No parameter reestimation.
        
        @param data: DataSet object
        @param labels: optional sample IDs
        @param entropy_cutoff: entropy threshold for the posterior distribution. Samples which fall
        above the threshold will remain unassigned
        @param silent: 0/1 flag, toggles verbose output        
        @param EStep: function implementing the EStep, by default self.EStep
        @param EStepParam: additional paramenters for more complex EStep implementations
        
        @return: list of class labels
        """
        if  isinstance(data, DataSet):
            if data.internalData is None:
                if not silent:
                    sys.stdout.write("Parsing data set...")
                    sys.stdout.flush()
                data.internalInit(self)
                if not silent:
                    sys.stdout.write("done\n")
                    sys.stdout.flush()
        else:    
            raise ValueError, "Invalid input type format: " + str(data.__class__)+", DataSet required."
       
        if EStep == None:
            EStep = self.EStep

        labels = data.sampleIDs
        
        # compute posterior distribution of component membership for cluster assignment
        [l,log_l] = EStep(data,EStepParam = EStepParam)
        
        if not silent:
            print "classify loglikelihood: "+str(log_l)+".\n"
 
        # cluster assingments initialised with -1
        z = numpy.ones(data.N,dtype='Int32') * -1       

        entropy_list = numpy.zeros(data.N,dtype='Float64')
        max_entropy = math.log(self.G,2)

        # compute posterior entropies 
        for i in range(data.N):
            exp_l = numpy.exp(l[:,i])
            if self.G == 1:
                entropy_list[i]  = entropy(exp_l)
            else:
                entropy_list[i]  = entropy(exp_l)/max_entropy

            if not entropy_cutoff:
                entropy_cutoff = float('inf')

#            if not silent:
#                print 'sample', data.sampleIDs[i],':',entropy_list[i]

            # apply entropy cutoff
            if entropy_list[i] < entropy_cutoff:
                # cluster assignment by maximum likelihood over the component membership posterior
                z[i] = numpy.argmax(l[:,i])
        
        if not silent:
            # printing out the clusters
            cluster = {}
            en = {}
            for j in range(-1,self.G,1):  
                cluster[j]=[]
                en[j]=[]

            for i in range(data.N):
                cluster[z[i]].append(labels[i])
                en[z[i]].append(entropy_list[i])
        
            print "\n** Clustering **"
            for j in range(self.G):
                print "Cluster ",j,', size',len(cluster[j])
                print cluster[j], "\n" 
        
            print "Unassigend due to entropy cutoff:"
            print cluster[-1], "\n"

        return z


    def isValid(self,x):
        """
        Exhaustive check whether a given DataSet is compatible with the model.
        If self is a lower hierarchy mixture 'x' is a single data sample in external representation.
        """
        if isinstance(x,DataSet):
            # self is at top level of the hierarchy
            for i in range(self.G):
                for j in range(x.N):
                    try:
                        self.components[i].isValid(x.dataMatrix[j])
                    except InvalidDistributionInput, ex:
                        ex.message += "\n\tin MixtureModel.components["+str(i)+"] for DataSet.dataMatrix["+str(j)+"]."
                        raise
        else:            
            for i in range(self.G):
                try:
                    self.components[i].isValid(x)
                except InvalidDistributionInput, ex:
                    ex.message += "\n\tMixtureModel.components["+str(i)+"]."
                    raise

    def formatData(self,x):
        [new_p, res] = self.components[0].formatData(x)
        return [new_p, res]

    def reorderComponents(self,order):
        """
        Reorder components into a new order
        
        @param order: list of indices giving the new order
        """
        # reordering components
        components_order = []
        pi_order = []
        compFix_order= []
        
        index_map = {}   # maps old indices to new indices
        for k,i in enumerate(order):    
            index_map[i] = k
            pi_order.append(self.pi[i])
            components_order.append( self.components[i]  )
            compFix_order.append(self.compFix[i])

        # assigning ordered parameters
        self.pi = numpy.array(pi_order,dtype='Float64')
        self.components = components_order
        self.compFix = compFix_order

        order_leaders = []
        order_groups = []
        # updating structure if necessary
        if self.struct:
            f = lambda x: index_map[x]
            for j in range(self.dist_nr):
                new_l = map(f,self.leaders[j])
                order_leaders.append(new_l)
                d = {}
                for h in self.leaders[j]:
                    new_g = map(f,self.groups[j][h])
                    d[index_map[h]] = new_g
                order_groups.append(d)    

            # reformat leaders and groups such that the minimal 
            # index of a group is used as the leader
            for j in range(self.dist_nr):
                for lind,lead in enumerate(order_leaders[j]):
                    tg = [lead] + order_groups[j][lead]

                    tg.sort()
                    nl = tg.pop(0)
                    order_leaders[j][lind] = nl
                    order_groups[j].pop(lead)
                    
                    order_groups[j][nl] = tg
                order_leaders[j].sort()

        self.leaders = order_leaders            
        self.groups = order_groups

    def identifiable(self):
        """ To provide identifiability the components are ordered by the mixture coefficient in 
            ascending order.
        """
        indices = {}
        # storing indices for sorting
        for i in range(self.G):
          indices[i] = self.pi[i]

        # determine new order of components by ascending mixture weight
        items = [(v, k) for k, v in indices.items()] 
        items.sort()
        
        order = []
        for it in items:
            order.append(it[1])
        self.reorderComponents(order)

    def flatStr(self,offset):
        offset +=1
        s = "\t"*offset+str(';Mix') +";"+str(self.G) + ";" + str(self.pi.tolist())+";"+str(self.compFix)+"\n"
        for c in self.components:   
            s += c.flatStr(offset)
        
        return s

    def printClusterEntropy(self,data):
        """
        Print out cluster stability measured by the entropy of the component membership posterior.
        
        @param data: DataSet object
        """
        if isinstance(data, numpy.ndarray):
            sequence = data
            seqLen = len(sequence)
                
        elif  isinstance(data, DataSet):
            sequence = data.internalData
       
        print "-------- getClusterEntropy ------------"
        post_entropy = []
        log_l = self.getPosterior(sequence)      
        l = numpy.exp(log_l)
        for i in range(data.N):
            post_entropy.append(entropy(l[:,i]))
        
        max_entropy = entropy([1.0/self.G]*self.G)
        print "Max entropy for G=",self.G,":",max_entropy
        
        print "--------\nPosterior distribuion: % max entropy"
        for i in range(data.N):
            print data.sampleIDs[i],": ",post_entropy[i] ," -> ",post_entropy[i] / max_entropy ,"%"
        print "--------\n"
        
    def posteriorTraceback(self,x):
        return self.pdf(x)[0] 

    def printTraceback(self,data,z,en_cut=1.01):
        """
        Prints out the posterior traceback, i.e. a detailed account of the
        contribution to the component membership posterior of each sample in each feature ordered
        by a clustering.
        
        @param data: DataSet object
        @param z: class labels
        @param en_cut: entropy threshold
        """
        if isinstance(data, numpy.ndarray):
            sequence = data
            seqLen = len(sequence)
        elif  isinstance(data, DataSet):
            templist = []
            for i in range(data.N):
                [t,dat] = self.components[0].formatData(data.dataMatrix[i])
                templist.append(dat)
            sequence = numpy.array(templist)
            labels = data.sampleIDs
            seqLen = len(sequence)

        l,log_p = self.EStep(data)

        print "seqLen = ", data.p
        print "pi = ",self.pi
        max_en = entropy([1.0/self.G]*self.G)
        for c in range(self.G):
            temp = numpy.where(z == c)
            c_index = temp[0]
            print "\n---------------------------------------------------"
            print "Cluster = ",c,": ",c_index
            for j in c_index:
                t = self.pdf(numpy.array([sequence[j]]))[0]
                print "\nj = ",j,", id =",data.sampleIDs[j],", log_l = ",t," -> ",numpy.exp(t)
                print "posterior = ",numpy.exp(l[:,j]).tolist(),"\n"
                tb = []
                for g in range(self.G):
                    ll = self.components[g].posteriorTraceback(data.internalData[j] )
                    tb.append(ll)
                tb_arr = numpy.array(tb,dtype='Float64')
                for i in range(len(tb[0])):
                    s = sumlogs(tb_arr[:,i])                
                    tb_arr[:,i] = tb_arr[:,i] - s
                exp_arr = numpy.exp(tb_arr)
                exp_tb = exp_arr.tolist()
                max_comp = numpy.zeros(len(tb[0]))                
                en_percent = numpy.zeros(len(tb[0]),dtype='Float64')
                for i in range(len(tb[0])):
                    max_comp[i] = numpy.argmax(tb_arr[:,i])
                    en_percent[i] = entropy(exp_arr[:,i]) / max_en
                print "     ",
                for q in range(len(tb[0])):
                    if en_percent[q] < en_cut:
                        head_len = len(data.headers[q])
                        print " " * (16-head_len) + str(data.headers[q]),
                print 
                print "     ",
                for q in range(len(tb[0])):
                   if en_percent[q] < en_cut:
                        x_len = len(str(data.dataMatrix[j][q]))
                        print " " * (16-x_len) + str(data.dataMatrix[j][q]),
                print 
                for e in range(self.G):
                    print e,": [",
                    if e != z[j]:
                        for k in range(len(exp_tb[e])):
                            if en_percent[k] < en_cut:
                                print "%16.10f" % (exp_tb[e][k],),
                    else:
                        for k in range(len(exp_tb[e])):
                            if en_percent[k] < en_cut:
                                print " *  %12.10f" % (exp_tb[e][k],),
                    print "]" 
                print "max  ",
                for e in range(len(data.headers)):
                    if en_percent[e] < en_cut:
                        print " "*15+str(max_comp[e]),
                print
                print "%EN  ",
                for e in range(len(data.headers)):
                    if en_percent[e] < en_cut:
                        print "%16.4f" % (en_percent[e],),
                print


    def update_suff_p(self):
        suff_p = self.components[0].update_suff_p()
        for i in range(1,self.G):
            sp =  self.components[i].update_suff_p() 
            assert sp == suff_p
        self.suff_p = suff_p
        return self.suff_p

    def updateStructureGlobal(self, data, silent=1): 
        """
        Updating CSI structure by chosing smallest KL distance merging, optimizing the AIC score.
        This was the first approach implemented for the CSI structure learning and using the Bayesian approach instead
        is stronly recommended.
        
        @param data: DataSet object
        @param silent: verbosity flag
        
        @return: number of structure changes 
        """
        assert self.struct == 1, "No structure in model."

        datRange = self.components[0].suff_dataRange
        new_leaders = []
        new_groups = []
        change = 0
        
        # building posterior factor matrix for the current group structure      
        l = numpy.zeros( (self.G,data.N,self.dist_nr ),dtype='Float64' )
        for j in range(self.dist_nr):
            if j == 0:
                prev = 0
            else:
                prev = datRange[j-1]     
            for lead_j in self.leaders[j]:
                if self.components[lead_j][j].suff_p == 1:
                    l_row = self.components[lead_j][j].pdf(data.internalData[:,datRange[j]-1] )
                else:
                    l_row = self.components[lead_j][j].pdf(data.internalData[:,prev:datRange[j]] )                
                l[lead_j, :, j] = l_row
                for v in self.groups[j][lead_j]:
                    l[v,:,j] = l_row
                                        
        g = numpy.sum(l,2) 
        for k in range(self.G):
            g[k,:] += numpy.log(self.pi[k])
        sum_logs = numpy.zeros(data.N,dtype='Float64')    
        for n in range(data.N):
            sum_logs[n] = sumlogs(g[:,n])
        lk = sum(sum_logs)
        for j in range(self.dist_nr):
            # initialize free parameters
            full_fp_0 = self.freeParams
            if not silent:
                print "\n************* j = ",j,"*****************\n"
            term = 0            
            while not term:
                nr_lead = len(self.leaders[j])
                if nr_lead == 1:
                    break
                min_dist = float('inf')
                merge_cand1 = -1
                merge_cand2 = -1
                # compute symmetric KL distances between all groups
                for i in range(nr_lead):
                    for k in range(i+1,nr_lead):
                        d = sym_dist(self.components[self.leaders[j][i]][j],self.components[self.leaders[j][k]][j])
                        if d < min_dist:
                            min_dist = d
                            merge_cand1 = self.leaders[j][i]
                            merge_cand2 = self.leaders[j][k]
                if not silent:
                    print "-------------------"
                    print merge_cand1," -> ",merge_cand2," = ",min_dist
                    print self.components[merge_cand1][j]
                    print self.components[merge_cand2][j]
                        
                full_BIC_0 = -2*lk + (full_fp_0  * numpy.log(data.N))
                # compute merged distribution of candidates with minimal KL distance
                candidate_dist = copy.copy(self.components[merge_cand1][j])
                merge_list = [ self.components[merge_cand2][j] ]
               
                # computing total weights for the two groups to be merged
                pi_list = [self.pi[merge_cand1], self.pi[merge_cand2]]
                for m in self.groups[j][merge_cand1]:
                    pi_list[0] += self.pi[m]
                for m in self.groups[j][merge_cand2]:
                    pi_list[1] += self.pi[m]
                
                # computing candidate leader distribution for the merged group
                candidate_dist.merge(merge_list, pi_list)

                if not silent:
                    print "candidate:", candidate_dist

                # computing the new reduced model complexity
                full_fp_1 = full_fp_0 - self.components[merge_cand1][j].freeParams

                # initialising new group structure with copies of the current structure
                new_leaders = copy.deepcopy(self.leaders)
                new_groups =  copy.deepcopy(self.groups)

                # removing merged leader from self.leaders
                ind = new_leaders[j].index(merge_cand2)
                new_leaders[j].pop(ind)
            
                # joining merged groups and removing old group entry
                new_groups[j][merge_cand1] += [merge_cand2] + new_groups[j][merge_cand2]
                new_groups[j].pop(merge_cand2)             
            
                if not silent:
                    print "\ncandidate model structure:"
                    print "lead = ",new_leaders[j]
                    print "groups = ",new_groups[j],"\n"
		
        		# updating likelihood matrix
                l_1 = copy.copy(l)
                if j == 0:
                    prev = 0
                else:
                    prev = datRange[j-1]     
                if candidate_dist.suff_p == 1:
                    l_row = candidate_dist.pdf(data.internalData[:,datRange[j]-1] )
                else:
                    l_row = candidate_dist.pdf(data.internalData[:,prev:datRange[j]] )                
                l_1[merge_cand1, :, j] = l_row
                for v in new_groups[j][merge_cand1]:
                    l_1[v,:,j] = l_row
                g = numpy.sum(l_1,2) 
                for k in range(self.G):
                    g[k,:] += numpy.log(self.pi[k])
            
                sum_logs = numpy.zeros(data.N,dtype='Float64')    
                for n in range(data.N):
                    sum_logs[n] = sumlogs(g[:,n])
                lk_1 = sum(sum_logs)
                full_BIC_1 = -2*lk_1 + (full_fp_1  * numpy.log(data.N))
                AIC_0 = -2*lk + ( 2 * full_fp_0 )
                AIC_1 = -2*lk_1 + ( 2 * full_fp_1 )

                if not silent:                
                    print "LK_0: ", lk
                    print "LK_1: ", lk_1

                    print "full_fp_0 =",full_fp_0
                    print "full_fp_1 =",full_fp_1
                
                    print "\nfull_BIC_0 =",full_BIC_0
                    print "full_BIC_1 =",full_BIC_1
                
                    print "Vorher: AIC_0 =",AIC_0
                    print "Nachher: AIC_1 =",AIC_1
            
                    #if  AIC_1 < AIC_0:
                    #    print "Merge accepted according to AIC"
                    #else:
                    #    print "Merge rejected according to AIC"
                
                if AIC_1 < AIC_0:
                    if not silent:
                        print "\n*** Merge accepted !"
                    change += 1

                    # new_model.components[merge_cand1][j]
                    # assigning leader distribution to group members
                    self.components[merge_cand1][j] = copy.copy(candidate_dist)
                    for g in new_groups[j][merge_cand1]:
                        self.components[g][j] = self.components[merge_cand1][j]
                    
                    self.leaders = new_leaders
                    self.groups = new_groups
                    
                    lk = lk_1
                    l = l_1
                    
                    full_fp_0 = full_fp_1

                    self.freeParams = full_fp_1                    

                    # if only one group is left terminate, update free parameters and continue with next variable
                    if len(self.leaders[j]) == 1:
                        if not silent:
                            print "*** Fully merged !"
                        term = 1
                else:
                    if not silent:
                        print "\n*** Merge rejected: Abort !"
                    # finished with this variable, update free parameters and go on
                    term = 1
        
       # reassinging groups and leaders 
        for j in range(self.dist_nr):
           for l in self.leaders[j]:
               for g in self.groups[j][l]:
                   self.components[g][j]  = self.components[l][j]
        return change

    def minimalStructure(self):
        """ Finds redundant components in the model structure and collapses the 
            structure to a minimal representation.
        """
        assert self.struct == 1, "No structure in model."
        
        distNr = self.components[0].dist_nr
        # features with only one group can be excluded
        exclude = []
        for i in range(distNr):
            if len(self.leaders[i]) == 1:
                exclude.append(i)
        # get first feature with more than one group
        first = -1
        for i in range(distNr):
            if i not in exclude:
                first = i
                break
        # initialising group dictionaries for first non-trivial group structure
        firstgroup_dicts = []   
        for j in range(len(self.leaders[first])):
            d = {}
            for k in [self.leaders[first][j]] + self.groups[first][self.leaders[first][j]]:
                d[k] = 1
            firstgroup_dicts.append(d)    

        # initialising group dictionaries for remaining features
        allgroups_dicts = []        
        for i in range(first+1,distNr,1):
            if i in exclude:
                continue
            gdicts_i = []
            for l in self.leaders[i]:
                d = {}
                for k in [l] + self.groups[i][l]:
                    d[k]  = 1
                gdicts_i.append(d)    
            allgroups_dicts.append(gdicts_i)

        toMerge = []
        # for each group in first non-trivial structure
        for g, fdict in enumerate(firstgroup_dicts):
            candidate_dicts = [fdict]
            # for each of the other non-trivial features
            for i, dict_list in enumerate(allgroups_dicts):
                new_candidate_dicts = [] 
                # for each group in the i-th feature
                for adict in dict_list:
                    # for each candidate group
                    for j,cand_dict in enumerate(candidate_dicts):
                        # find intersection
                        inter_d = dict_intersection(cand_dict, adict)       
                        if len(inter_d) >= 2:
                            new_candidate_dicts.append(inter_d)
                candidate_dicts = new_candidate_dicts
                # check whether any valid candidates are left
                if len(candidate_dicts) == 0:
                    break
            if len(candidate_dicts) > 0:
                for c in candidate_dicts:
                    toMerge.append(c.keys())

        if len(toMerge) > 0:  # postprocess toMerge to remove single entry sets
            for i in range(len(toMerge)-1,-1,-1):
                if len(toMerge[i]) == 1:
                    toMerge.pop(i)

        d_merge = None
        if len(toMerge) > 0:
            d_merge = {}  # map from indices to be merged to respective leader index
            for m in range(len(toMerge)):
                tm = toMerge[m]
                tm.sort()
                lead = tm.pop(0)
                for x in tm:
                    d_merge[x] = lead
            new_pi = self.pi.tolist()
            new_compFix = copy.copy(self.compFix)
            
            l = d_merge.keys()
            l.sort()
            l.reverse()
            for j in l:
                
                # update new_pi
                pi_j = new_pi.pop(j)
                new_pi[d_merge[j]] += pi_j
                self.components.pop(j)  # remove component
                
                # update compFix
                cf_j = new_compFix.pop(j)
                new_compFix[d_merge[j]] = new_compFix[d_merge[j]] or cf_j
                
                # update leaders
                for l1 in range(len(self.leaders)):
                    for l2 in range(len(self.leaders[l1])):
                        if self.leaders[l1][l2] == j:
                            self.leaders[l1][l2] = d_merge[j]
                        elif self.leaders[l1][l2] > j:   
                            self.leaders[l1][l2] -= 1

                # update component indices in groups
                for g_j in range(len(self.groups)):
                    for g in self.groups[g_j].keys():
                        if g == j:
                            tmp = self.groups[g_j][g]
                            self.groups[g_j].pop(j)          
                            self.groups[g_j][d_merge[j]] = tmp
                        elif g > j:
                            tmp = self.groups[g_j][g]
                            self.groups[g_j].pop(g)          
                            self.groups[g_j][g-1] = tmp

                # remove merged component from groups
                for g_j in range(len(self.groups)):
                    for g in self.groups[g_j].keys():
                        for gm in range(len(self.groups[g_j][g])-1,-1,-1):
                            if self.groups[g_j][g][gm] == j:
                                self.groups[g_j][g].pop(gm)
                            elif self.groups[g_j][g][gm] > j:
                                self.groups[g_j][g][gm] -= 1

            self.G = self.G - len(l)  # update number of components
            self.pi = numpy.array(new_pi,dtype='Float64')  # set new pi in model
            self.compFix = new_compFix

        self.updateFreeParams()
        if self.identFlag:
            self.identifiable()

        return d_merge
 
    def removeComponent(self,ind):
        """
        Deletes a component from the model.

        @param ind: ind of component to be removed
        """
        self.G = self.G - 1  # update number of components
        tmp = self.pi.tolist() # update pi
        tmp.pop(ind)
        tmp = map(lambda x: x / sum(tmp),tmp) # renormalize pi
        self.pi = numpy.array(tmp,dtype='Float64')  # set new pi in model
        self.components.pop(ind)  # remove component
        if self.compFix:
            self.compFix.pop(ind)

        # update CSI structure if necessary
        if self.struct:
            #update leaders
            for k,ll in enumerate(self.leaders):
                try:  # remove ind from leader lists
                    r = ll.index(ind)
                    self.leaders[k].pop(r)
                except ValueError:
                    pass
                for i,l in enumerate(self.leaders[k]): # update to new indices
                    if l > ind:
                        self.leaders[k][i] -= 1

            new_groups = []
            # update groups
            for i,dg in enumerate(self.groups):
                new_groups.append({})
                for k in dg.keys():
                    if k == ind:
                        # case ind is leader: remove ind and select new leader
                        gr = self.groups[i].pop(k)
                        if len(gr) > 0:  # need to re-enter with new leader
                            new_l = gr.pop(0)
                            if new_l > ind: new_l -= 1
                            for r in range(len(gr)):
                                if gr[r] > ind:
                                    gr[r] -=1
                            new_groups[i][new_l] = gr
                            self.leaders[i].append(new_l)
                    else:
                        # case ind is not leader but might be in the group
                        gr = self.groups[i].pop(k)
                        if ind in gr:
                            gr.remove(ind)
                        for r in range(len(gr)):
                            if gr[r] > ind:
                                gr[r] -=1
                        if k > ind: 
                            new_groups[i][k-1] = gr  # need to change key                        
                        else:
                            new_groups[i][k] = gr
        self.groups = new_groups
        self.updateFreeParams()
        if self.identFlag:
            self.identifiable()

    def merge(self,dlist, weights):
        raise DeprecationWarning, 'Part of the outdated structure learning implementation.'
        coeff = 1.0 / ( len(dlist) +1 )
        m_pi = self.pi * coeff
        for i in range(len(dlist)):
            assert isinstance(dlist[i],MixtureModel)
            assert dlist[i].G == self.G
        for j in range(self.G):
            group = []
            for i in range(len(dlist)):
                group.append( dlist[i].components[j] )
                m_pi[j] += coeff * dlist[i].pi[j]
            self.components[j].merge(group,weights)
        self.pi = m_pi

    def printStructure(self,data= None):
        """
        Pretty print of the model structure
        """
        assert self.struct == 1, "No model structure."
        if data:
            headers = data.headers
        else:
            headers = range(self.dist_nr)    
        for i in range(self.dist_nr):
            print "Feature "+str(i)+": " + str(headers[i])
            for j,l in enumerate(self.leaders[i]):
                if self.groups[i][l] == []:
                    print "\tGroup "+str(j)+": "+"("+str(l)+")"
                else:
                    print "\tGroup "+str(j)+": "+str(tuple([l]+self.groups[i][l]))
                print "\t  ",self.components[l][i],"\n"
    

    def updateFreeParams(self):
        """
        Updates the number of free parameters for the current group structure
        """
        if self.struct == 0:
            self.freeParams = (self.components[0].freeParams * self.G) + self.G-1
        else:
            fp = 0            
            for i in range(self.dist_nr):
                for l in self.leaders[i]:
                    fp += self.components[l][i].freeParams
            fp += self.G -1
            self.freeParams = fp
            

    def validStructure(self):
        """
        Checks whether the CSI structure is syntactically correct. Mostly for debugging.
        """
        if self.struct == 0:
            return True
        r = range(self.G)
        try:
            # check valid entries in group and leader
            for j in range(self.dist_nr):
                for l in self.leaders[j]:
                    assert l in r
                    for g in self.groups[j][l]:
                        assert g in r
            # check completeness of structure
            for j in range(self.dist_nr):
                tmp = copy.copy(self.leaders[j])
                for g in self.groups[j].keys():
                    tmp += copy.copy(self.groups[j][g])
                tmp.sort()
                assert tmp == r
        except AssertionError:
            print 'Invalid structure:',j
            print self.leaders
            print self.groups
            raise

    def sufficientStatistics(self, posterior, data):
        """
        Returns sufficient statistics for a given data set and posterior.
        """
        assert self.dist_nr == 1
        sub_post = get_posterior(self, data,logreturn=True)
        
        suff_stat = []
        dat = data.getInternalFeature(0)
        for i in range(self.G):
            if self.compFix[i] == 2:
                suff_stat.append([float('-inf'),float('inf')])
                continue

            np =  sub_post[i]+posterior
            inds = numpy.where(np != float('-inf'))
            suff_stat.append( self.components[i][0].sufficientStatistics( np[inds], dat[inds]) )
        
        return suff_stat



#--------------------------- TEST ------------------------------------------------
class CandidateGroupHISTORY:  # XXX reomve ? ...
    def __init__(self, l, dist_prior,dist ):
        #self.indices = indices  # leader indices for this merge
        #self.post = post   # total posterior of this merge
        self.l = l       # vector of likelihoods of the merge for each sample in a single feature
        self.dist_prior = dist_prior  # prior density of candidate distribution
        self.dist = dist  # candidate distribution

        #self.lead = lead  # candidate leaders
        #self.groups = groups  # candidate groups
       
        #self.l_j_1 = l_j_1
        #self.log_prior_list_j = log_prior_list_j


class CandidateGroup:  
    """
    CandidateGroup is a simple container class. 
    It holds the parameters and sufficient statistics for a candidate grouping
    in the CSI structure. It is used as part of the structure learning.
    """

    def __init__(self, dist, post_sum, pi_sum , req_stat, l=None, dist_prior=None):
        """
        Constructor
        
        @param dist: candidate distribution
        @param post_sum: sum over component membership posterior for candidate distribution
        @param pi_sum:  sum of pi's corresponding to candidate distribution
        @param req_stat:  additional statistics required  for paramter updates
        @param l: vector of likelihoods induced by the candidate distribution for each sample in a single feature
        @param dist_prior:  prior density over candidate distribution
        """
        
        self.dist = dist  # candidate distribution
        self.post_sum = post_sum  # sum over posterior for candidate merge
        self.pi_sum = pi_sum   # sum of pi's corresponding to candidate merge
        self.req_stat = req_stat  # additional required statistics for paramter updates by merge
        
        self.l = l       # vector of likelihoods of the merge for each sample in a single feature
        self.dist_prior = dist_prior  # prior density of candidate distribution
        
#----------------------------------------------------------------------------------

            
class BayesMixtureModel(MixtureModel):
    """
    Bayesian mixture models
    """
    def __init__(self,G, pi, components, prior,compFix=None,struct=0,identifiable = 1):
        """
        Constructor

        @param G: number of components
        @param pi: mixture weights
        @param components: list of ProductDistribution objects, each entry is one component
        @param prior: MixtureModelPrior object
        @param compFix: list of optional flags for fixing components in the reestimation
                         the following values are supported:
                         1 distribution parameters are fixed, 
                         2 distribution parameters and mixture coefficients are fixed
        @param struct: Flag for CSI structure, 
            0 = no CSI structure
            1 = CSI structure
        """
        MixtureModel.__init__(self,G, pi, components, compFix=compFix, struct=struct, identifiable = identifiable)

        # check and set model prior
        self.prior = None
        prior.isValid(self)
        self.prior = prior

    def __str__(self):
        s = MixtureModel.__str__(self)
        s += "\n" + str(self.prior)
        return s

    def __eq__(self,other):
        if not isinstance(other, BayesMixtureModel):
            return False
        res = MixtureModel.__eq__(self,other)
        if res == False:
            return res
        else:
            res = self.prior.__eq__(other.prior)    
            return res

    def __copy__(self):
        copy_components = []
        copy_pi = copy.deepcopy(self.pi) 
        copy_compFix = copy.deepcopy(self.compFix)
        for i in range(self.G):
            copy_components.append(copy.deepcopy(self.components[i]))
        copy_prior = copy.copy(self.prior)
        copy_model = BayesMixtureModel(self.G, copy_pi, copy_components, copy_prior, compFix = copy_compFix)
        copy_model.nr_tilt_steps = self.nr_tilt_steps    
        copy_model.suff_p = self.suff_p
        copy_model.identFlag = self.identFlag
        
        if self.struct:
            copy_model.initStructure()
            copy_leaders = copy.deepcopy(self.leaders)
            copy_groups = copy.deepcopy(self.groups)
            copy_model.leaders = copy_leaders
            copy_model.groups = copy_groups
        return copy_model
    
    def modelInitialization(self,data,rtype=1):
        """
        Perform model initialization given a random assigment of the
        data to the models.
            
        @param data: DataSet object
        @param rtype: type of random assignments.
        0 = fuzzy assignment
        1 = hard assignment

        @return: posterior assigments
        """
        if not isinstance(data,DataSet):
            raise TypeError, "DataSet object required, got"+ str(data.__class__)
        else:
            if data.internalData is None:
                data.internalInit(self)
        # reset structure if applicable
        if self.struct:
            self.initStructure()
        # generate 'random posteriors'
        l = numpy.zeros((self.G,len(data)),dtype='Float64')
        for i in range(len(data)):
            if rtype == 0:
              for j in range(self.G):
                 l[j,i] = random.uniform(0.1,1)
              s = sum(l[:,i])
              for j in range(self.G):
                  l[j,i] /= s
            else:
               l[random.randint(0,self.G-1),i] = 1 

        # do one M Step
        fix_pi = 1.0
        unfix_pi = 0.0
        fix_flag = 0   # flag for fixed mixture components
        for i in range(self.G):
                # setting values for pi
                if self.compFix[i] == 2:
                    fix_pi -= self.pi[i]
                    fix_flag = 1
                else:
                    #self.pi[i] =  l[i,:].sum() / len(data)                
                    self.pi[i] =  ( l[i,:].sum() + self.prior.piPrior.alpha[i] -1.0 ) / (len(data) + ( self.prior.piPrior.alpha_sum - self.G) ) 
                    unfix_pi += self.pi[i]
                if self.compFix[i] == 1 or self.compFix[i] == 2:
                    # fixed component
                    continue
                else:    
                    # components are product distributions that may contain mixtures
                    if isinstance(self.components[i], ProductDistribution):
                        last_index = 0
                        for j in range(self.components[i].dist_nr):    
                            if isinstance(self.components[i][j],MixtureModel):
                                dat_j = data.singleFeatureSubset(j)
                                self.components[i][j].modelInitialization(dat_j,rtype=rtype)      
                            else:    
                                loc_l = l[i,:]
                                # masking missing values from parameter estimation
                                if data.missingSymbols.has_key(j):
                                    ind_miss = data.getMissingIndices(j)
                                    for k in ind_miss:
                                        loc_l[k] = 0.0
                                self.prior.compPrior[j].mapMStep(self.components[i][j], loc_l,data.getInternalFeature(j), dist_ind = i )  
                    else:  # components are not ProductDistributions -> invalid    
                        raise TypeError

        # renormalizing mixing proportions in case of fixed components
        if fix_flag:
            if unfix_pi == 0.0:
                #print "----\n",self,"----\n"
                 print "unfix_pi = ", unfix_pi
                 print "fix_pi = ", fix_pi
                 print "pi = ", self.pi
                 raise RuntimeError, "unfix_pi = 0.0"  
            for i in range(self.G):
                if self.compFix[i] == 0:
                    self.pi[i] = (self.pi[i] * fix_pi) / unfix_pi

        # updating hyperparameters in prior where apprpriate
        if self.prior.constant_hyperparams != 1:
            self.prior.updateHyperparameters(self, l, data)

        return l

    def mapEM(self,data, max_iter, delta,silent = False, mix_pi=None, mix_posterior= None, tilt = 0):
        return MixtureModel.mapEM(self,data, self.prior, max_iter, delta,silent = silent, mix_pi=mix_pi, mix_posterior= mix_posterior, tilt = tilt)
 
    def removeComponent(self,ind):
        """
        Deletes a component from the model.

        @param ind: index of component to be removed 
        """
        MixtureModel.removeComponent(self,ind)
        # update component prior
        alpha = self.prior.piPrior.alpha.tolist()
        alpha.pop(ind)
        n_pipr = DirichletPrior(self.G,alpha)
        self.prior.piPrior = n_pipr


    # XXX CODE REVIEW


    def updateStructureBayesian(self,data,objFunction='MAP',silent=1):
        """
        Updating structure by chosing optimal local merge with respect to the posterior.

        
        Features: - store merges in a learning history to prevent recomputation of merge parameters
                  - compute parameters of candidate structures from parameters of groups to be merged
        
       
        @param data: DataSet object
        @param silent: verbosity flag

        @return: number of structure changes 
        """
        assert self.struct == 1, "No structure in model."
        assert objFunction in ['MAP'] # for now only MAP estimation

        new_leaders = []
        new_groups = []
        change = 0
        # building data likelihood factor matrix for the current group structure      
        l = numpy.zeros( (self.dist_nr, self.G, data.N),dtype='Float64' )
        for j in range(self.dist_nr):
            # extracting current feature from the DataSet
            if isinstance(self.components[0][j], MixtureModel): # XXX
                data_j = data.singleFeatureSubset(j)
            else:
                data_j = data.getInternalFeature(j)
            
            for lead_j in self.leaders[j]:
                l_row = self.components[lead_j][j].pdf(data_j)  
                l[j,lead_j,:] = l_row
                for v in self.groups[j][lead_j]:
                    l[j,v,:] = l_row

        # g is the matrix of log posterior probabilities of the components given the data
        g = numpy.sum(l, axis=0) 
        for k in range(self.G):
            g[k,:] += numpy.log(self.pi[k])


        sum_logs = matrixSumlogs(g)
        
        try:
            g_norm = g - sum_logs
        except FloatingPointError:
            print sum_logs
            raise
            
        tau = numpy.exp(g_norm)

        if not silent:
            print "\ntau="
            for tt in tau:
                print tt.tolist()
            print

        # computing posterior as model selection criterion
        temp = DiscreteDistribution(self.G,self.pi)
        pi_prior = self.prior.piPrior.pdf(temp)
        log_prior = pi_prior
        log_prior_list = [0.0] * self.dist_nr
        for j in range(self.dist_nr):
            for r in range(self.G):
                   log_prior_list[j] += self.prior.compPrior[j].pdf( self.components[r][j] )
        log_prior += sum(log_prior_list)

        # prior over number of components
        log_prior += self.prior.nrCompPrior * self.G 
        # prior over number of distinct groups
        for j in range(self.dist_nr):
            log_prior += self.prior.structPrior * len(self.leaders[j])   

        # get posterior
        lk = numpy.sum(sum_logs) 
        post = lk + log_prior
        if not silent:
            print "0: ",  lk ,"+", log_prior,"=", post
            print log_prior_list

        changes = 0        
        g_wo_j = numpy.zeros((self.G, data.N),dtype='Float64')

        # initialising temporary group structure with copies of the current structure
        temp_leaders = copy.deepcopy(self.leaders)
        temp_groups =  copy.deepcopy(self.groups)
        for j in range(self.dist_nr):
            L = {}  # initialize merge history

            if not silent:
                print "\n************* j = ",j,"*****************\n"

            # unnormalized posterior matrix without the contribution of the jth feature
            try:
                g_wo_j  = g - l[j]
            except FloatingPointError:

                # if there was an exception we have to compute each
                # entry in g_wo_j seperately to set -inf - -inf = -inf
                g_wo_j = _C_mixextend.substract_matrix(g,l[j])

            # checking whether feature j is already fully merged 
            nr_lead = len(self.leaders[j])
            if nr_lead == 1:
                continue  # nothing to be done...

            term = 0            
            if not silent:
                print self.leaders
                print self.groups

            # extracting current feature from the DataSet
            if isinstance(self.components[0][j], MixtureModel): # XXX
                data_j = data.singleFeatureSubset(j)
            else:
                data_j = data.getInternalFeature(j)

            # initialize merge history
            tau_pool = numpy.zeros(data.N, dtype='Float64')

            for lead in self.leaders[j]:
                #el_dist = copy.copy(self.components[lead][j])
                
                # NOT a copy, changes in el_dist changes self !
                el_dist = self.components[lead][j]
                tau_pool = copy.copy(tau[lead, :])
                pi_pool = self.pi[lead]
                for z in self.groups[j][lead]:
                    tau_pool += tau[z, :]                    
                    pi_pool += self.pi[z]

                if objFunction == 'MAP':
                    self.prior.compPrior[j].mapMStep(el_dist, tau_pool, data_j,  pi_pool)  
                else:
                    # should never get here...
                    raise TypeError

                stat = el_dist.sufficientStatistics(tau_pool, data_j)
                M = CandidateGroup(el_dist, numpy.sum(tau_pool), pi_pool, stat)
                l_row = el_dist.pdf(data_j)  
                cdist_prior = self.prior.compPrior[j].pdf( el_dist ) 
                M.l = l_row
                M.dist_prior = cdist_prior

                L[(lead,)+tuple(self.groups[j][lead])] = M

                # update likelihood matrix for initial model
                for ll in [lead] + self.groups[j][lead]:
                    l[j,ll,:] = l_row

            while not term:
                best_dist = None   # accepted candidate distributions
                best_post = float('-inf')   # corresponding posteriors
                best_indices = None
                best_l_j = l[j] 
                for mc1 in range(len(temp_leaders[j])):
                    merge_cand1 = temp_leaders[j][mc1]
                    for mc2 in range(mc1+1,len(temp_leaders[j])):
                        merge_cand2 = temp_leaders[j][mc2]
                        if not silent:
                            print "-------------------"
                            print merge_cand1," -> ",merge_cand2
                            print self.components[merge_cand1][j], '( sum(tau) = ',sum(tau[merge_cand1,:]),')'
                            print self.components[merge_cand2][j], '( sum(tau) = ',sum(tau[merge_cand2,:]),')'

                        nr_leaders_j = len(temp_leaders[j])-1
                        cand_group_j = temp_groups[j][merge_cand1] + [merge_cand2] + temp_groups[j][merge_cand2]

                        hist_ind_part1 = (merge_cand1,)+tuple(temp_groups[j][merge_cand1])
                        hist_ind_part2 = (merge_cand2,)+tuple(temp_groups[j][merge_cand2])
                        hist_ind_complete = hist_ind_part1 + hist_ind_part2
                        
                        recomp = 0 
                        
                        if L.has_key(hist_ind_complete):
                            recomp = 1
                        
                        if not silent:
                            print "\ncandidate model structure: "
                            print 'merge:',hist_ind_part1, hist_ind_part2, '->' ,hist_ind_complete
                            #print "lead = ",leaders_j
                            #print "groups = ",groups_j
                            #print "others = ",others,"\n"

                        if not recomp:
                            assert L.has_key( hist_ind_part1),str(hist_ind_part1)+' missing.'
                            assert L.has_key( hist_ind_part2),str(hist_ind_part2)+' missing.'
                            
                            M = self.prior.compPrior[j].mapMStepMerge([L[hist_ind_part1], L[hist_ind_part2]])
                            candidate_dist = M.dist

                            if not silent:
                                print "candidate:", candidate_dist

                            l_row = candidate_dist.pdf(data_j)  
                            cdist_prior = self.prior.compPrior[j].pdf( candidate_dist ) 
                           
                            M.l = l_row
                            M.dist_prior = cdist_prior
                            L[hist_ind_complete] = M
                        else:
                            # retrieve merge data from history
                            candidate_dist = L[hist_ind_complete].dist                            

                            if not silent:
                                print "candidate:", candidate_dist

                            l_row = L[hist_ind_complete].l
                            cdist_prior = L[hist_ind_complete].dist_prior

                	    # computing change in likelihood matrix for this step
                        l_j_1 = copy.copy(l[j]) 

                        # updating l_j_1 with the new candidate distribution
                        l_j_1[merge_cand1, :] = l_row
                        for v in cand_group_j:
                            l_j_1[v,:] = l_row

                        # get updated unnormalized posterior matrix
                        _C_mixextend.add_matrix(g, g_wo_j, l_j_1)
                        sum_logs = matrixSumlogs(g)
                        lk_1 = numpy.sum(sum_logs)

                        # computing posterior as model selection criterion
                        log_prior_1 = pi_prior

                        # compute parameter prior for the candidate merge parameters
                        log_prior_list_j = 0.0
                        for r in range(self.G):
                            if r in [merge_cand1] + cand_group_j:
                                log_prior_list_j += cdist_prior
                            else:
                                log_prior_list_j += self.prior.compPrior[j].pdf( self.components[r][j] )

                        log_prior_1 += sum(log_prior_list)
                        log_prior_1 -= log_prior_list[j]
                        log_prior_1 += log_prior_list_j

                        # prior over number of components
                        log_prior_1 += self.prior.nrCompPrior * self.G 
                        # prior over number of distinct groups
                        for z in range(self.dist_nr):
                            if z == j:
                                log_prior_1 += self.prior.structPrior * nr_leaders_j  
                            else:
                                log_prior_1 += self.prior.structPrior * len(temp_leaders[z])   # XXX len could be cached ?

                        post_1 = lk_1 + log_prior_1

                        if not silent:
                            print '\nPosterior:',post_1 ,'=', lk_1 ,'+', log_prior_1

                        if post_1 >= post:
                            if not silent:
                                print "*** Merge accepted", post_1 ,">=", post

                            if post_1 > best_post:  # current merge is better than previous best
                                best_dist = candidate_dist
                                best_post = post_1
                                best_indices = [merge_cand1, merge_cand2]
                                best_l_j = l_j_1
                                best_log_prior_list_j = log_prior_list_j
                        else:
                            if not silent:
                                print "*** Merge rejected:", post_1 ,"!>", post

                # if there is no possible merge that increases the score we are done
                if best_post == float('-inf'):
                    if not silent:
                        print "*** Finished !"
                    # setting updated structure in model
                    self.leaders[j] = temp_leaders[j]
                    self.groups[j] = temp_groups[j]

                    # reset posterior matrix to the last accepted merge
                    g = g_wo_j + best_l_j
                    term = 1

                # otherwise we update the model with the best merge found
                else:  
                    if not silent:
                        print "\n--- Winner ---"
                        print "post:", best_post
                        print "indices:",best_indices
                        print "dist: ",best_dist
                        #print "lead:",best_leaders
                        #print "group:", best_groups

                    post = best_post
                    l[j] = best_l_j
                    g = g_wo_j + best_l_j  # posterior matrix for the next iteration
                    log_prior_list[j] = best_log_prior_list_j

                    # updating model
                    # removing merged leader from new_leaders
                    ind = temp_leaders[j].index(best_indices[1])
                    temp_leaders[j].pop(ind)

                    # joining merged groups and removing old group entry
                    temp_groups[j][best_indices[0]] += [best_indices[1]] + temp_groups[j][best_indices[1]]
                    temp_groups[j].pop(best_indices[1]) 

                    # assigning distributions according to new structure
                    self.components[best_indices[0]][j] = best_dist
                    for d in temp_groups[j][best_indices[0]]:
                        self.components[d][j] = self.components[best_indices[0]][j]
                    change += 1

        return change 

#---------------------------------------------------------------------------------------------------------------

    def updateStructureBayesianFullEnumerationFixedOrder(self,data,objFunction='MAP',silent=1):
        """
        Updating structure by chosing optimal local merge with respect to the posterior.

        Enumerates all possible structures for each feature seperately, i.e. returns the optimal structure
        for the given ordering of features. For a complete enumeration of the structure space, including permutation of
        feature order use the functions in fullEnumerationExhaustive.py.
       
        @param data: DataSet object
        @param silent: verbosity flag

        @return: number of structure changes 
        """
        assert self.struct == 1, "No structure in model."
        assert objFunction in ['MAP'] # for now only MAP estimation

        new_leaders = []
        new_groups = []
        change = 0
        # building data likelihood factor matrix for the current group structure      
        l = numpy.zeros( (self.dist_nr, self.G, data.N),dtype='Float64' )
        for j in range(self.dist_nr):
            for lead_j in self.leaders[j]:
                l_row = self.components[lead_j][j].pdf(data.getInternalFeature(j) )  
                l[j,lead_j,:] = l_row
                for v in self.groups[j][lead_j]:
                    l[j,v,:] = l_row

        # g is the matrix of log posterior probabilities of the components given the data
        g = numpy.sum(l, axis=0) 
        for k in range(self.G):
            g[k,:] += numpy.log(self.pi[k])

        sum_logs = matrixSumlogs(g)
        g_norm = g - sum_logs
        tau = numpy.exp(g_norm)
        
        if not silent:
            print "\ntau="
            for tt in tau:
                print tt.tolist()
            print

        # computing posterior as model selection criterion
        temp = DiscreteDistribution(self.G,self.pi)
        pi_prior = self.prior.piPrior.pdf(temp)
        log_prior = pi_prior
        log_prior_list = [0.0] * self.dist_nr
        for j in range(self.dist_nr):
            for r in range(self.G):
                   log_prior_list[j] += self.prior.compPrior[j].pdf( self.components[r][j] )
        log_prior += sum(log_prior_list)

        # prior over number of components
        log_prior += self.prior.nrCompPrior * self.G 
        # prior over number of distinct groups
        for j in range(self.dist_nr):
            log_prior += self.prior.structPrior * len(self.leaders[j])   

        # get posterior
        lk = numpy.sum(sum_logs) 
        best_post = lk + log_prior
        if not silent:
            print "0: ",  lk ,"+", log_prior,"=", best_post
            print log_prior_list

        changes = 0        
        g_wo_j = numpy.zeros((self.G, data.N),dtype='Float64')

        # initialising temporary group structure with copies of the current structure
        #temp_leaders = copy.deepcopy(self.leaders)
        #temp_groups =  copy.deepcopy(self.groups)
        for j in range(self.dist_nr):
            L = {}  # initialize merge history

            if not silent:
                print "\n************* j = ",j,"*****************\n"

            # unnormalized posterior matrix without the contribution of the jth feature
            try:
                g_wo_j  = g - l[j]
            except FloatingPointError:
                # if there was an exception we have to compute each
                # entry in g_wo_j seperately to set -inf - -inf = -inf
                g_wo_j = _C_mixextend.substract_matrix(g,l[j])


            # checking whether feature j is already fully merged 
            nr_lead = len(self.leaders[j])
            if nr_lead == 1:
                continue  # nothing to be done...

            term = 0            
            if not silent:
                print self.leaders
                print self.groups,'\n'
                
            # extracting current feature from the DataSet
            if isinstance(self.components[0][j], MixtureModel): # XXX
                data_j = data.singleFeatureSubset(j)
            else:
                data_j = data.getInternalFeature(j)

            for lead in self.leaders[j]:
                el_dist = copy.copy(self.components[lead][j])

                tau_pool = copy.copy(tau[lead, :])
                pi_pool = self.pi[lead]
                for z in self.groups[j][lead]:
                    tau_pool += tau[z, :]                    
                    pi_pool += self.pi[z]

                if objFunction == 'MAP':
                    self.prior.compPrior[j].mapMStep(el_dist, tau_pool, data_j,  pi_pool)  
                else:
                    # should never get here...
                    raise TypeError
           
                stat = el_dist.sufficientStatistics(tau_pool, data_j)
                M = CandidateGroup(el_dist, numpy.sum(tau_pool), pi_pool, stat)

                l_row = el_dist.pdf(data_j)  
                cdist_prior = self.prior.compPrior[j].pdf( el_dist ) 

                M.l = l_row
                M.dist_prior = cdist_prior
                L[(lead,)+tuple(self.groups[j][lead])] = M
                
            # first partition is full structure matrix
            kappa, max_kappa = init_last(self.G)  # XXX first partition posterior is computed twice ! 

            best_dist = None   # accepted candidate distributions
            best_indices = None
            best_l_j = l[j] 
            best_log_prior_list_j = log_prior_list[j]
            best_partition = [ (ll,)+tuple(self.groups[j][ll]) for ll in self.leaders[j] ]
            best_post = float('-inf')
            
            while 1:
                curr_part = setPartitions.decode_partition(numpy.arange(self.G) ,kappa,max_kappa)
                if not silent:
                    print "\n-------------------"
                    #print 'History:', L.keys()
                    print 'Current structure:',kappa,'->',j ,curr_part
                    
                    fullstruct = []
                    for jj in range(self.dist_nr):
                        if jj != j:
                            fullstruct.append( [tuple( [ll] + self.groups[jj][ll] ) for ll in self.leaders[jj]])
                        else:   
                            fullstruct.append(curr_part)
                    print 'Full:' ,fullstruct

                # computing change in likelihood matrix for this step
                #l_j_1 = copy.copy(l[j]) 
                l_j_1 = numpy.zeros( (self.G, data.N ) )  # XXX needs only be done once

                for group in curr_part:
                    if L.has_key(group):
                        # retrieve merge data from history
                        candidate_dist = L[group].dist                            

                        if not silent:
                            print "  candidate:", group, candidate_dist

                        l_row = L[group].l
                        cdist_prior = L[group].dist_prior
                    else:
                        M = self.prior.compPrior[j].mapMStepMerge( [L[(c,)]  for c in group]  )
                        candidate_dist = M.dist

                        if not silent:
                            print "  candidate:", group,candidate_dist

                        l_row = candidate_dist.pdf(data_j)  
                        cdist_prior = self.prior.compPrior[j].pdf( candidate_dist ) 

                        M.l = l_row
                        M.dist_prior = cdist_prior

                        L[group] = M
 
                    for g in group:
                        l_j_1[g,:] = l_row
                                
                # get updated unnormalized posterior matrix
                g = g_wo_j + l_j_1
                sum_logs = matrixSumlogs(g)
                lk_1 = numpy.sum(sum_logs)

                # computing posterior as model selection criterion
                log_prior_1 = pi_prior

                # compute parameter prior for the candidate merge parameters
                log_prior_list_j = 0.0
                for r in curr_part:
                    log_prior_list_j += L[r].dist_prior * len(r)

                log_prior_1 += sum(log_prior_list)
                log_prior_1 -= log_prior_list[j]
                log_prior_1 += log_prior_list_j

                # prior over number of components
                log_prior_1 += self.prior.nrCompPrior * self.G 
                # prior over number of distinct groups
                for z in range(self.dist_nr):
                    if z == j:
                        log_prior_1 += self.prior.structPrior * len(curr_part) 
                    else:
                        log_prior_1 += self.prior.structPrior * len(self.leaders[z])   # XXX len could be cached ?
                post_1 = lk_1 + log_prior_1

                if not silent:
                    print '\nPosterior:',post_1 ,'=', lk_1 ,'+', log_prior_1

                if post_1 >= best_post: # current candidate structure is better than previous best
                    if not silent:
                        print "*** New best candidate", post_1 ,">=", best_post
                    best_post = post_1
                    best_partition = curr_part # XXX
                    best_l_j = l_j_1
                    best_log_prior_list_j = log_prior_list_j

                else:
                    if not silent:
                        print "*** worse than previous best",best_partition,'(', post_1 ,"!>", best_post,')'
                
                ret = prev_partition(kappa, max_kappa)
                if ret == None:  # all candidate partitions have been scored
                    if not silent:
                        print "*** Finished with post=",best_post
                    # setting updated structure in model
                    lead = []
                    groups = {}
                    for gr in best_partition:
                        gr_list = list(gr)
                        gr_lead = gr_list.pop(0)
                        lead.append(gr_lead)
                        groups[gr_lead] = gr_list
                        
                        # assigning distributions according to new structure
                        self.components[gr_lead][j] = L[gr].dist
                        for d in gr_list:
                            self.components[d][j] = self.components[gr_lead][j]
                    
                    self.leaders[j] = lead
                    self.groups[j] = groups

                    # reset posterior matrix to the last accepted merge
                    g = g_wo_j + best_l_j
                    log_prior_list[j] = best_log_prior_list_j
                    break
                
                kappa,max_kappa  = ret
                
                
#---------------------------------------------------------------------------------------------------------------

    def updateStructureBayesianBottomUp(self,data,objFunction='MAP',silent=1):
        """
        Updating structure by chosing optimal local merge with respect to the posterior. The procedure
        starts with the minimally complex structure, i.e. every feature has a single distribution.

       
        @param data: DataSet object
        @param silent: verbosity flag

        @return: number of structure changes 
        """
        assert self.struct == 1, "No structure in model."
        assert objFunction in ['MAP'] # for now only MAP estimation

        new_leaders = []
        new_groups = []
        change = 0
        # building data likelihood factor matrix for the current group structure      
        l = numpy.zeros( (self.dist_nr, self.G, data.N),dtype='Float64' )
        for j in range(self.dist_nr):
            for lead_j in self.leaders[j]:
                l_row = self.components[lead_j][j].pdf(data.getInternalFeature(j) )  
                l[j,lead_j,:] = l_row
                for v in self.groups[j][lead_j]:
                    l[j,v,:] = l_row

        # g is the matrix of log posterior probabilities of the components given the data
        g = numpy.sum(l, axis=0) 
        for k in range(self.G):
            g[k,:] += numpy.log(self.pi[k])

        sum_logs = matrixSumlogs(g)
        g_norm = g - sum_logs
        tau = numpy.exp(g_norm)
        
        if not silent:
            print "\ntau="
            for tt in tau:
                print tt.tolist()
            print

        # computing posterior as model selection criterion
        temp = DiscreteDistribution(self.G,self.pi)
        pi_prior = self.prior.piPrior.pdf(temp)
        
        # compute feature wise parameter prior contributions 
        log_prior_list = [0.0] * self.dist_nr
        for j in range(self.dist_nr):
            for r in range(self.G):
                   log_prior_list[j] += self.prior.compPrior[j].pdf( self.components[r][j] )

        changes = 0        
        g_wo_j = numpy.zeros((self.G, data.N),dtype='Float64')

        # initialising starting group structure
        temp_leaders = copy.copy(self.leaders)
        temp_groups = copy.copy(self.groups)
        for j in range(self.dist_nr):
            L = {}  # initialize merge history

            if not silent:
                print "\n************* j = ",j,"*****************\n"

            # initialising starting group structure for feature j
            temp_leaders[j] =  [0]
            temp_groups[j] = {0: range(1,self.G)}

            # unnormalized posterior matrix without the contribution of the jth feature
            try:
                g_wo_j  = g - l[j]
            except FloatingPointError:
                # if there was an exception we have to compute each
                # entry in g_wo_j seperately to set -inf - -inf = -inf
                g_wo_j = _C_mixextend.substract_matrix(g,l[j])

            term = 0            

            if not silent:
                print temp_leaders
                print temp_groups

            # extracting current feature from the DataSet
            if isinstance(self.components[0][j], MixtureModel): # XXX
                data_j = data.singleFeatureSubset(j)
            else:
                data_j = data.getInternalFeature(j)

            # initial model structure
            tau_pool = numpy.ones(data.N, dtype='Float64')
            pi_pool = 1.0
            el_dist = copy.copy(self.components[0][j])

            if objFunction == 'MAP':
                self.prior.compPrior[j].mapMStep(el_dist, tau_pool, data_j,  pi_pool)  
            else:
                # should never get here...
                raise TypeError
            
            # set initial distribution in model
            for i in range(self.G):
                self.components[i][j] = el_dist
            
            stat = el_dist.sufficientStatistics(tau_pool, data_j)
            M = CandidateGroup(el_dist, numpy.sum(tau_pool), pi_pool, stat)

            l_row = el_dist.pdf(data_j)  
            cdist_prior = self.prior.compPrior[j].pdf( el_dist ) 

            M.l = l_row
            M.dist_prior = cdist_prior
            L[tuple(range(self.G))] = M

            sum_logs = matrixSumlogs( g_wo_j + l_row)
            
            temp = copy.copy(l_row)
            temp = temp.reshape(1,data.N)
            l[j] = temp.repeat(self.G,axis=0)

            # get likelihood
            lk = numpy.sum( sum_logs )

            log_prior = pi_prior
            log_prior_list[j] = self.G * self.prior.compPrior[j].pdf( el_dist)
            log_prior += numpy.sum(log_prior_list)

            # prior over number of components
            log_prior += self.prior.nrCompPrior * self.G 
            # prior over number of distinct groups
            for jj in range(self.dist_nr):
                log_prior += self.prior.structPrior * len(temp_leaders[jj])   
            
            post = lk + log_prior
            if not silent:
                print "0: ",  lk ,"+", log_prior,"=", post

            split_dist = copy.copy(self.components[0][j])  
            while not term:
                best_split_dist = None   # accepted split distributions
                best_remainder_dist = None   # accepted remainder distributions
                best_log_prior_list_j = log_prior_list[j]                
                best_post = float('-inf')   # corresponding posteriors
                best_leaders = None
                best_groups = None
                #best_indices = None
                best_l_j = l[j] 
                for mc1 in range(len(temp_leaders[j])):
                    merge_cand1 = temp_leaders[j][mc1]

                    if len(temp_groups[j][merge_cand1]) == 0:
                        continue    # nothing to split

                    for mc1_grp in temp_groups[j][merge_cand1]:
                        if not silent:
                            print "-------------------"
                            print '*** leader '+str(merge_cand1)+': split',mc1_grp,'from',temp_groups[j][merge_cand1]
                            print self.components[merge_cand1][j], '( sum(tau) = ',sum(tau[merge_cand1,:]),')'
                            
                        # initialising candidate group structure with copies of the temporary structure
                        leaders_j = copy.copy(temp_leaders[j])
                        groups_j =  copy.deepcopy(temp_groups[j])

                        hist_ind_presplit = (merge_cand1,)+tuple(groups_j[merge_cand1])

                        # adding new leader created by split to leaders_j
                        leaders_j.append( mc1_grp)

                        # removing new leader from merge_cand1 group
                        ind = groups_j[merge_cand1].index(mc1_grp)
                        groups_j[merge_cand1].pop(ind)
                        
                        # adding new group 
                        groups_j[mc1_grp] = []
                        
                        if not silent:
                            print "\ncandidate model structure:"
                            print "lead = ",leaders_j
                            print 'groups =',groups_j
                            print j,[ (ll,)+tuple(groups_j[ll]) for ll in leaders_j ]
                            
                        nr_leaders_j = len(temp_leaders[j])-1
                        
                        hist_ind_remainder = (merge_cand1,)+tuple(groups_j[merge_cand1])
                        hist_ind_split = (mc1_grp,)

                	    # computing change in likelihood matrix for this step
                        l_j_1 = copy.copy(l[j]) 
                        
                        # computing parameters for new single component group
                        if objFunction == 'MAP':
                            self.prior.compPrior[j].mapMStep(split_dist, tau[mc1_grp,:], data_j,  self.pi[mc1_grp])  
                        else:
                            # should never get here...
                            raise TypeError
                        if not silent:
                            print "split dist:", split_dist
                        
                        # updating l_j_1 with the new split distribution
                        l_row = split_dist.pdf(data_j)    
                        l_j_1[mc1_grp, :] = l_row
                        
                        # add candidategroup 
                        stat = split_dist.sufficientStatistics(tau[mc1_grp,:], data_j)
                        M = CandidateGroup(split_dist, numpy.sum(tau[mc1_grp,:]), self.pi[mc1_grp], stat)
                        L[hist_ind_split] = M

                        split_dist_prior = self.prior.compPrior[j].pdf( split_dist ) 

                        M.l = l_row
                        M.dist_prior = split_dist_prior

                        # computing parameters for group which has been split
                        recomp = 0
                        
                        if L.has_key(hist_ind_remainder):
                            recomp = 1

                        if not recomp:
                            assert L.has_key( hist_ind_presplit),str(hist_ind_presplit)+' missing.'
                            assert L.has_key( hist_ind_split),str(hist_ind_split)+' missing.'
                            M = self.prior.compPrior[j].mapMStepSplit(L[hist_ind_presplit], L[hist_ind_split])
                            remainder_dist = M.dist

                            if not silent:
                                print "remainder dist:", remainder_dist

                            l_row = remainder_dist.pdf(data_j)  
                            remainder_dist_prior = self.prior.compPrior[j].pdf( remainder_dist ) 
                           
                            M.l = l_row
                            M.dist_prior = remainder_dist_prior
                            
                            L[hist_ind_remainder] = M
                        else:
                            # retrieve merge data from history
                            remainder_dist = L[hist_ind_remainder].dist                            

                            if not silent:
                                print "remainder dist:", remainder_dist

                            l_row = L[hist_ind_remainder].l
                            remainder_dist_prior = L[hist_ind_remainder].dist_prior

                        # updating l_j_1 with the new remainder distribution
                        l_j_1[merge_cand1, :] = l_row
                        for v in groups_j[merge_cand1]:
                            l_j_1[v,:] = l_row

                        sum_logs = matrixSumlogs(g)
                        lk_1 = numpy.sum(sum_logs)

                        # computing posterior as model selection criterion
                        log_prior_1 = pi_prior

                        # compute parameter prior for the candidate merge parameters
                        log_prior_list_j = 0.0
                        for r in range(self.G):
                            if r in [merge_cand1] + groups_j[merge_cand1]:
                                log_prior_list_j += remainder_dist_prior
                            elif r == mc1_grp:
                                log_prior_list_j += split_dist_prior
                            else:
                                log_prior_list_j += self.prior.compPrior[j].pdf( self.components[r][j] )

                        log_prior_1 += sum(log_prior_list)
                        log_prior_1 -= log_prior_list[j]
                        log_prior_1 += log_prior_list_j

                        # prior over number of components
                        log_prior_1 += self.prior.nrCompPrior * self.G 
                        # prior over number of distinct groups
                        for z in range(self.dist_nr):  
                            if z == j:
                                log_prior_1 += self.prior.structPrior * (len(temp_leaders[z])+1)
                            else:
                                log_prior_1 += self.prior.structPrior * len(temp_leaders[z])   # XXX len could be cached ?

                        post_1 = lk_1 + log_prior_1
                        if not silent:
                            print '\nPosterior:',post_1 ,'=', lk_1 ,'+', log_prior_1

                        if post_1 >= post:
                            if not silent:
                                print "*** Split accepted", post_1 ,">=", post

                            if post_1 > best_post:  # current merge is better than previous best
                                best_split_dist = copy.copy(split_dist)
                                best_remainder_dist = copy.copy(remainder_dist)
                                best_post = post_1
                                best_leaders = leaders_j
                                best_groups = groups_j

                                best_indices = [merge_cand1, mc1_grp]
                                #best_l_j = copy.copy(l_j_1)
                                best_l_j = l_j_1
                                best_log_prior_list_j = log_prior_list_j
                        else:
                            if not silent:
                                print "*** Split rejected:", post_1 ,"!>", post

                # if there is no possible split that increases the score we are done
                if best_post == float('-inf'):
                    if not silent:
                        print "*** Finished with post", post
                    # setting updated structure in model
                    self.leaders[j] = temp_leaders[j]
                    self.groups[j] = temp_groups[j]

                    # reset posterior matrix to the last accepted merge
                    l[j] = best_l_j
                    g = g_wo_j + best_l_j
                    term = 1

                # otherwise we update the model with the best merge found
                else:  
                    if not silent:
                        print "\n--- Winner ---"

                        print "post:", best_post
                        print "indices:",best_indices
                        print "remainder dist: ",best_remainder_dist
                        print "split dist: ",best_split_dist

                        print "lead:",best_leaders
                        print "group:", best_groups

                    post = best_post
                    l[j] = best_l_j
                    g = g_wo_j + best_l_j  # posterior matrix for the next iteration
                    log_prior_list[j] = best_log_prior_list_j

                    # updating model
                    temp_leaders[j] = best_leaders
                    temp_groups[j] = best_groups

                    # assigning distributions according to new structure
                    self.components[best_indices[0]][j] = copy.copy(best_remainder_dist)
                    
                    for d in temp_groups[j][best_indices[0]]:
                        self.components[d][j] = self.components[best_indices[0]][j]
                        
                    self.components[best_indices[1]][j] = copy.copy(best_split_dist)
                    change += 1

        return change 

# --------------------------------------------------------------------------------------------------------------------

    def KLFeatureRanks(self,data, comps, silent=False):
        """
        Ranks the features by the symmetric relative entropy between the parameters induced in a structure
        where all components in 'comps' are merged and the parameters induced by all components (i.e. the uninformative
        structure).

        This gives a ranking of the relevance,i.e. discriminatory information of features for distinguishing a
        subset of components.
        
        @param data: DataSet object
        @param comps: list of component indices
        
        @return: list of tuples of (feature index, score) pairs in descending score order 
        """
        assert type(comps) == list  # XXX for debugging mostly
        
        # some initial checks
        for c in comps:
            assert c in range(self.G) # check for valid entries 
        assert len(comps) < self.G  # all components doesn`t make sense
        
        comps.sort()

        others = range(0,self.G) 
        for c in comps:
            others.remove(c)        

        if not silent:
            print 'comps',comps            
            print 'others:',others
        
        # building data likelihood factor matrix for the current group structure      
        l = numpy.zeros( (self.dist_nr, self.G, data.N),dtype='Float64' )
        for j in range(self.dist_nr):
            if isinstance(self.components[0][j], MixtureModel):
                data_j = data.singleFeatureSubset(j)
            else:
                data_j = data.getInternalFeature(j)              
            for lead_j in self.leaders[j]:
                l_row = self.components[lead_j][j].pdf(data_j )  
                l[j,lead_j,:] = l_row
                for v in self.groups[j][lead_j]:
                    l[j,v,:] = l_row

        # g is the matrix of log posterior probabilities of the components given the data
        g = numpy.sum(l,axis=0) 

        for k in range(self.G):
            g[k,:] += numpy.log(self.pi[k])

        sum_logs = numpy.zeros(data.N,dtype='Float64')    
        g_norm = numpy.zeros((self.G, data.N),dtype='Float64')
        for n in range(data.N):
            sum_logs[n] = sumlogs(g[:,n])
            # normalizing log posterior
            g_norm[:,n] = g[:,n] - sum_logs[n]

        tau = numpy.exp(g_norm)
        model = copy.copy(self)
        score = []
        for j in range(model.dist_nr):
            if len(model.leaders[j]) == 1:
                if not silent:
                    print 'Feature '+str(data.headers[j])+' uninformative.'
                # if the feature is uninformative already the score is set to zero    
                score.append( (0.0,j) )
                continue
            else:
                # this whole section is more general than needed, might get optimized as some point
                
                if not silent:
                    print '\nFeature '+str(data.headers[j])+' ( index '+str(j)+' ) usefull. '
            
                # data for the jth feature
                if isinstance(model.components[0][j], MixtureModel):
                    data_j = data.singleFeatureSubset(j)
                else:
                    data_j = data.getInternalFeature(j)              

                # updating groups and leaders:
                # backup leaders and groups of feature j
                backup_leaders = copy.copy(model.leaders[j])
                backup_groups = copy.copy(model.groups[j])

                # clear structure for j
                model.groups[j] = {}

                # assign structure for the ranking
                model.leaders[j] = [ comps[0],others[0] ]
                model.groups[j][comps[0]] = comps[1:]
                model.groups[j][others[0]] = others[1:]

                # saving distribution parameters for feature j
                temp = []            
                for i in range(model.G):
                    temp.append(model.components[i][j])

                # estimating distribution for component to be ranked
                comps_dist_j = copy.copy(model.components[0][j])
                tau_pool = copy.copy(tau[comps[0],:])
                pi_pool = self.pi[comps[0]]
                for z in model.groups[j][comps[0]]:
                    tau_pool += tau[z, :]                    
                    pi_pool += self.pi[z]
                model.prior.compPrior[j].mapMStep(comps_dist_j, tau_pool, data_j,  pi_pool)  

                # estimating distribution for all components except the one to be ranked
                others_dist_j = copy.copy(model.components[0][j])
                tau_pool = copy.copy(tau[others[0],:])
                pi_pool = self.pi[others[0]]
                for z in model.groups[j][others[0]]:
                    tau_pool += tau[z, :]                    
                    pi_pool += self.pi[z]
                model.prior.compPrior[j].mapMStep(others_dist_j, tau_pool, data_j,  pi_pool)  

                # compute feature score
                score.append( ( sym_kl_dist(comps_dist_j, others_dist_j), j) )

        score.sort()
        score.reverse()
        return score         


    def randMaxTraining(self, data,nr_runs,nr_steps,delta,tilt=0,objFunction='MAP',rtype=1,silent=False):
        """
        Performs `nr_runs` MAP EM runs with random initial parameters
        and returns the model which yields the maximum likelihood.
        
        @param data: DataSet object
        @param nr_runs: number of repeated random initializations
        @param nr_steps: maximum number of steps in each run
        @param delta: minimim difference in log-likelihood before convergence
        @param tilt: 0/1 flag, toggles the use of a deterministic annealing in the training
        @param silent:0/1 flag, toggles verbose output
        
        @return: log-likelihood of winning model
        """
        if isinstance(data, numpy.ndarray):
            raise TypeError, "DataSet object required."
        elif  isinstance(data, DataSet):
            if data.internalData is None:
                sys.stdout.write("Parsing data set...")
                sys.stdout.flush()
                data.internalInit(self)
                sys.stdout.write("done\n")
                sys.stdout.flush()
        else:    
            raise ValueError, "Unknown input type format: " + str(data.__class__)

        assert objFunction in ['MAP']  # only MAP estimation for now

        best_logp = float('-inf')
        best_model = None
        logp_list = []
        for i in range(nr_runs):
            # copying inside of the loop is necessary for the dirichlet 
            # models to get independent parameter initializations
            candidate_model = copy.copy(self)  # copying the model parameters 
            # we do repeated random intializations until we get a model with valid posteriors in all components
            init = 0
            while not init:
                try:        
                    candidate_model.modelInitialization(data,rtype=rtype)    # randomizing parameters of the model copy
                except InvalidPosteriorDistribution:
                    pass
                else:
                    init = 1     
            try:        
                if objFunction == 'MAP': # maximum a posteriori estimation
                    (l,log_p) = candidate_model.mapEM(data, nr_steps,delta,silent=silent,tilt=tilt)  # running EM 
                else:
                    # should never get here... 
                    raise TypeError
                logp_list.append(log_p)

            except ConvergenceFailureEM:
                sys.stdout.write("Run "+str(i)+" produced invalid model, omitted.\n")
            except ValueError:
                sys.stdout.write("Run "+str(i)+" produced invalid model, omitted.\n")
            except InvalidPosteriorDistribution:
                sys.stdout.write("Run "+str(i)+" produced invalid model, omitted.\n")
            else:
                # current model is better than previous optimum
                if log_p > best_logp:
                    best_model = copy.copy(candidate_model)
                    best_logp = log_p

        if not silent:
            print "\nBest model likelihood over ",nr_runs,"random initializations ( "+str(nr_runs - len(logp_list))+" runs failed):"
            if len(logp_list) > 0:
                print "Model likelihoods:",logp_list
                print "Average logp: ", sum(logp_list)/len(logp_list)," SD:",numpy.array(logp_list).std()
                print "Best logp:",best_logp

        self.components = best_model.components  # assign best parameter set to model 'self'
        self.pi = best_model.pi
        return best_logp  # return final data log likelihood
                
    def shortInitMAP(self,data,nr_runs,nr_init,nr_steps,delta,tilt=0,nr_seed_steps=5,silent=False):
        """ 
        EM strategy:
            - 'nr_init' random initial models
            - short EM runs with each start model
            - take model with best likelihood for long EM run.
            
        The short runs are set to 5 iterations by default    
        
        @param data: DataSet object
        @param nr_runs: number of repeated random initializations
        @param nr_init: number of random models for each run
        @param nr_steps: maximum number of steps for the long training run
        @param delta: minimim difference in log-likelihood before convergence
        @param tilt: 0/1 flag, toggles the use of a deterministic annealing in the training
        @param nr_seed_steps: number of EM steps for each seed model, default is 5
        @param silent:0/1 flag, toggles verbose output
        
        @return: log-likelihood of winning model
        """
        if isinstance(data, numpy.ndarray):
            raise TypeError, "DataSet object required."
        elif  isinstance(data, DataSet):
            if data.internalData is None:
                sys.stdout.write("Parsing data set...")
                sys.stdout.flush()
                data.internalInit(self)
                sys.stdout.write("done\n")
                sys.stdout.flush()
        else:    
            raise ValueError, "Unknown input type format: " + str(data.__class__)
            
        seed_model = copy.copy(self)  # copying the model parameters 
        best_model = None
        best_model_logp = float('-inf')
        model_logp = []        
        for i in range(nr_runs):
            print "i = " , i
            best_seed_model = None
            best_seed_logp = float('-inf')
            seed_logp = []

            for j in range(nr_init):
                seed_model.modelInitialization(data)    # randomizing parameters of the model copy
                invalid_model = 0  # flag for models that produced an exception in EM
                try:
                    (l,logp_i) = seed_model.mapEM(data,nr_seed_steps,0.1,silent=silent)    
                except InvalidPosteriorDistribution:
                    print "Invalid seed model discarded: ",j
                    invalid_model = 1
                
                # only valid models are considered in the maximization
                if not invalid_model:     
                    if logp_i > best_seed_logp:
                        #print 'better',j
                        best_seed_logp = logp_i
                        best_seed_model = copy.copy(seed_model)
                    seed_logp.append(logp_i)
 
            sys.stdout.write("--- picking best seed model ---\n")
            try:
                (l,logp) = best_seed_model.mapEM(data,nr_steps,delta,silent=silent,tilt=tilt)
            except InvalidPosteriorDistribution:
                sys.stdout.write("***** Run "+str(i)+" produced invalid model.\n")
            except ZeroDivisionError:
                sys.stdout.write("***** Run "+str(i)+" is numerically instable.\n")  
            else:
                if logp > best_model_logp:
                    best_model_logp = logp
                    best_model = copy.copy(best_seed_model)
                model_logp.append(logp)

        if not silent:
            print "\nBest model likelihood over ",nr_runs," repeats ( "+str(nr_runs - len(model_logp))+" runs failed):"
            print "Model likelihoods:",model_logp
            print "Average logp: ", sum(model_logp)/len(model_logp)," SD:",numpy.array(model_logp).std()
            print "Best logp:",best_model_logp

        final_logp = numpy.array(model_logp,dtype='Float64')
        # assign best parameter set to model 'self' 
        self.components = best_model.components
        self.pi = best_model.pi
        return best_modellogp

    def bayesStructureEM(self,data,nr_repeats,nr_runs,nr_steps,delta,tilt=0,objFunction='MAP',silent=False,min_struct=1,rtype=1):
        """ 
        EM training for models with CSI structure.
        First a candidate model is generated by using the randMaxMAP procedure,
        then the structure is trained. 
            
        @param data: DataSet object
        @param nr_repeats: number of candidate models to be generated
        @param nr_runs: number of repeated random initializations
        @param nr_steps: maximum number of steps for the long training run
        @param delta: minimim difference in log-likelihood before convergence
        @param tilt: 0/1 flag, toggles the use of a deterministic annealing in the training
        @param silent:0/1 flag, toggles verbose output
        @param min_struct: 0/1 flag, toggles merging of components with identical paramters
        
        @return: log-likelihood of winning model
        """
        if isinstance(data, numpy.ndarray):
            raise TypeError, "DataSet object required."
        elif  isinstance(data, DataSet):
            if data.internalData is None:
                sys.stdout.write("Parsing data set...")
                sys.stdout.flush()
                data.internalInit(self)
                sys.stdout.write("done\n")
                sys.stdout.flush()
        else:    
            raise ValueError, "Unknown input type format: " + str(data.__class__)

        assert objFunction in ['MAP']  # only MAP estimation for now
        assert self.struct 
        logp_list = []
        best_logp = None
        best_candidate = None
        for r in range(nr_repeats):
            error = 0
            candidate_model = copy.copy(self)  # copying the model parameters 
            # we do repeated random intializations until we get a model with valid posteriors in all components
            init = 0
            while not init:
                try: 
                    #print "bayesStructureEM: try init   " 
                    candidate_model.modelInitialization(data, rtype=rtype)    # randomizing parameters of the model copy
                except InvalidPosteriorDistribution:
                    pass
                else:
                    init = 1     
            
            log_p = candidate_model.randMaxTraining(data,nr_runs,nr_steps,delta,tilt=tilt,objFunction=objFunction,rtype=rtype,silent=silent)
            
            try:
                ch = candidate_model.updateStructureBayesian(data,objFunction=objFunction,silent=1)  
            except ValueError:
                error = 1
                print 'ERROR: failed structure lerarning.'
                continue
                    
            if not silent:
                print "Changes = ",ch
            while ch != 0:
                try:
                    if objFunction == 'MAP':
                        candidate_model.mapEM(data, nr_steps,delta,silent=silent,tilt=0)
                    else:  
                        # should never get here...
                        raise TypeError 
                    ch = candidate_model.updateStructureBayesian(data,objFunction=objFunction,silent=1) 
                    if not silent:
                        print "Changes = ",ch
                except ConvergenceFailureEM:
                    print 'FAILURE:'
                    error = 1
                    break
                        
            try:
                # DEBUG: check structure validity
                self.validStructure()
            except AssertionError:
                print 'ERROR: Produced invalid structure'
                error = 1                
            if not error:
                try: 
                    if objFunction == 'MAP':
                        (l,log_p) =candidate_model.mapEM(data, nr_steps,delta,silent=silent,tilt=0)
                    else:  
                        # should never get here...
                        raise TypeError 
                except ConvergenceFailureEM:
                    continue
                logp_list.append(log_p)
                if r == 0 or log_p > best_logp:
                    best_logp = log_p
                    best_candidate = candidate_model
            else:
                continue

        self.components = best_candidate.components  # assign best parameter set to model 'self'
        self.pi = best_candidate.pi
        self.groups = best_candidate.groups
        self.leaders = best_candidate.leaders
        if min_struct:
            # remove redundant components
            self.minimalStructure()

        # update free parameters 
        self.updateFreeParams()
        if not silent:
            print 'Structural EM (',nr_repeats,' runs over',nr_runs,'random inits each):'
            print 'logp:',logp_list
            print "Average logp: ", sum(logp_list)/len(logp_list)," SD:",numpy.array(logp_list).std()
            print "Best logp:",best_logp
        return best_logp            



    def minimalStructure(self):
        d = MixtureModel.minimalStructure(self)
        # if there have been changes to the number of
        # components the prior needs to be updated
        if d is not None:
            l = d.keys()
            l.sort()
            l.reverse()
            alpha = self.prior.piPrior.alpha.tolist()
            for i in l:
                a = alpha.pop(i)
                alpha[d[i]] += (a-1.0)
            
            #print "alpha =",alpha
            assert len(alpha) == self.G
            new_piPr = DirichletPrior(self.G,alpha)
            self.prior.piPrior = new_piPr

        

#---------------------------------- Partial Supervised Learning -------------------------------------------------     
            
class ConstrainedDataSet(DataSet):
    """
    Extension of the DataSet object that can hold pairwise or label constraints in the objects.
    This data set is required  for the semi-supervised learning EM 
    """   
    def __init__(self):
        DataSet.__init__(self)
        self.pairwisepositive = None
        self.pairwisenegative = None
        self.labels = []
        self.noLabels = 0

    def __copy__(self):
        cop = DataSet.__copy(self)
        cop.labels = copy.copy(self.labels)
        cop.pairwisepositive = copy.copy(self.pairwisepositive)
        cop.pairwisenegative = copy.copy(self.pairwisenegative)
        cop.noLabels = self.noLabels

    def setConstrainedLabels(self,labels):         
        """
        Sets labels for semi-supervised learning.

        @param labels: list of lists of sample indices. The index in 'labels' denotes the component the samples are
        assigned to. For instance labels = [[0,2],[4,6]] would mean samples 0 and 2 are labelled with component 0.
        """
        assert sum([len(i) for i in labels]) <= self.N, 'Label constraints must be within the number of observations'
        self.labels = labels
        self.noLabels = len(labels)       

    def setPairwiseConstraints(self,positive,negative):
        """
        Set pairwise constraints.
        
        XXX add params
        """
        if positive != None:
          assert len(positive) == self.N, 'Pairwise Constraints should cover the all observations'
        if negative != None:
           assert len(negative) == self.N, 'Pairwise Constraints should cover the all observations'       
        self.pairwisepositive = positive
        self.pairwisenegative = negative


def LMMfromMM(mm):
    """
    Convenience function. Takes a MixtureModel or and returns a LabeledMixtureModel with
    the same parameters.
    
    @param mm: MixtureModel object
    """
    return LabeledMixtureModel(mm.G,mm.pi,mm.components,mm.compFix,mm.struct)


class LabeledMixtureModel(MixtureModel):

    """ 
    Class for a mixture model containing the label constrained
    version of the E-Step See 
    A. Schliep, C. Steinhoff,
    A. A. Schonhuth Robust inference of groups in gene expression
    time-courses using mixtures of HMMs Bioinformatics. 2004 Aug 4;20
    Suppl 1:I283-I289 (Proceedings of the ISMB 2004).
    for details 
    """
    def __init__(self,G, pi, components,compFix=None,struct=0):
       MixtureModel.__init__(self,G, pi, components, compFix=compFix, struct=struct,identifiable=0)

    def EM(self, data, max_iter, delta,silent = False, mix_pi=None, mix_posterior= None, tilt = 0):
        """
        Reestimation of mixture parameters using the EM algorithm.
        This method do some initial checking and call the EM from MixtureModel with the constrained labels E step
        
        @param data: DataSet object
        @param max_iter: maximum number of iterations
        @param delta: minimal difference in likelihood between two iterations before
        convergence is assumed.
        @param silent: 0/1 flag, toggles verbose output
        @param mix_pi: [internal use only] necessary for the reestimation of
        mixtures as components
        @param mix_posterior:[internal use only] necessary for the reestimation of
        mixtures as components
        @param tilt: 0/1 flag, toggles the use of a deterministic annealing in the training
        
        @return: tuple of posterior matrix and log-likelihood from the last iteration
        """
        
        assert  isinstance(data,ConstrainedDataSet), 'Data set does not contain labels, Labeled EM can not be performed'
        assert data.labels != None, 'Data set does not contain labels, Labeled EM can not be performed'
        assert data.noLabels <= self.G, 'Number of components should be equal or greatedr then the number of label classes'

        return MixtureModel.EM(self,data, max_iter, delta,silent = silent, mix_pi=mix_pi,
                               mix_posterior= mix_posterior,tilt = tilt, EStep=self.EStep, EStepParam = None)

    def EStep(self,data,mix_posterior=None,mix_pi=None,EStepParam=None):
        """Reestimation of mixture parameters using the EM algorithm.
                
        @param data: DataSet object
        @param mix_pi: [internal use only] necessary for the reestimation of
        mixtures as components
        @param mix_posterior:[internal use only] necessary for the reestimation of
        mixtures as components
        @param EStepParam: additional paramenters for more complex EStep implementations, in
        this implementaion it is ignored

        @return: tuple of log likelihood matrices and sum of log-likelihood of components
        """
        # computing log posterior distribution
        #[log_l,log_p] = MixtureModel.EStep(self,data,mix_posterior,mix_pi,None)

        log_l = numpy.zeros((self.G,data.N),dtype='Float64')
        log_col_sum = numpy.zeros(data.N,dtype='Float64')  # array of column sums of log_l
        log_pi = numpy.log(self.pi)  # array of log mixture coefficients

        # compute log of mix_posterior (if present)
        if mix_posterior is not None:
            log_mix_posterior = numpy.log(mix_posterior)

        # computing log posterior distribution 
        for i in range(self.G):            
            log_l[i] = log_pi[i] + self.components[i].pdf(data)

    	# Partially supervised training
        # For sequences with partial information use a weight vector s.t.
        # P[model|seq] = 1 if model = label[seq] and 0 else
        for i,cl in enumerate(data.labels): # for each class
            for o in cl: # for each observation in a class
              v = log_l[i,o]
              p_vec = numpy.zeros(self.G, dtype='Float64')
              p_vec[:] = float('-inf')
              p_vec[i] = v
              log_l[:,o] = p_vec

        log_col_sum = numpy.zeros(data.N,dtype='Float64')  # array of column sums of log_l
        for j in range(data.N):
            log_col_sum[j] = sumlogs(log_l[:,j]) # sum over jth column of log_l    
            # if posterior is invalid, check for model validity
            if log_col_sum[j] == float('-inf'): 
                # if self is at the top of hierarchy, the model is unable to produce the
                # sequence and an exception is raised. Otherwise normalization is not necessary.
                if mix_posterior is None and not mix_pi:
                    #print "\n---- Invalid -----\n",self,"\n----------"
                    #print "\n---------- Invalid ---------------"
                    #print "mix_pi = ", mix_pi
                    #print "x[",j,"] = ", data.getInternalFeature(j)
                    #print "l[:,",j,"] = ", log_l[:,j] 
                    #print 'data[',j,'] = ',data.dataMatrix[j]
                    raise InvalidPosteriorDistribution, "Invalid posterior distribution."
            # for valid posterior, normalize and go on    
            else:
                # normalizing log posterior
                log_l[:,j] = log_l[:,j] - log_col_sum[j]
                # adjusting posterior for lower hierarchy mixtures
                if mix_posterior is not None:
                    # multiplying in the posterior of upper hierarch mixture
                    log_l[:,j] = log_l[:,j] + log_mix_posterior[j]
        return log_l, numpy.sum(log_col_sum)

    def modelInitialization(self,data,rtype=1, missing_value = None):
        """
        Perform model initialization given a random assigment of the
        data to the models.
            
        @param data: DataSet object
        @param rtype: type of random assignments.
        0 = fuzzy assingment
        1 = hard assingment
        @param missing_value: missing symbol to be ignored in parameter estimation (if applicable)
        """
        assert  isinstance(data,ConstrainedDataSet), 'Data set does not contain labels, Labeled EM can not be performed'
        assert data.labels != None, 'Data set does not contain labels, Labeled EM can not be performed'
        assert data.noLabels <= self.G, 'Number of components should be equal or greater than the number of label classes'

        if not isinstance(data,DataSet):
            raise TypeError, "DataSet object required, got"+ str(data.__class__)
        else:
            if data.internalData is None:
                data.internalInit(self)
        # reset structure if applicable
        if self.struct:
            self.initStructure()

        # generate 'random posteriors'
        l = numpy.zeros((self.G,len(data)),dtype='Float64')
        for i in range(len(data)):
            if rtype == 0:
              for j in range(self.G):
                 l[j,i] = random.uniform(0.1,1)
              s = sum(l[:,i])
              for j in range(self.G):
                  l[j,i] /= s
            else:
               l[random.randint(0,self.G-1),i] = 1

        # peform label assigments (non random!)
        for i,cl in enumerate(data.labels): # for each class
            for o in cl: # for each observation in a class
              p_vec = numpy.zeros(self.G, dtype='Float64')
              p_vec[i] = 1.0
              l[:,o] = p_vec
              
        # do one M Step
        fix_pi = 1.0
        unfix_pi = 0.0
        fix_flag = 0   # flag for fixed mixture components
        for i in range(self.G):
                # setting values for pi
                if self.compFix[i] == 2:
                    fix_pi -= self.pi[i]
                    fix_flag = 1
                else:
                    self.pi[i] =  l[i,:].sum() / len(data)                
                    unfix_pi += self.pi[i]
                if self.compFix[i] == 1 or self.compFix[i] == 2:
                    # fixed component
                    continue
                else:    
                    last_index = 0
                    for j in range(self.components[i].dist_nr):    
                        if isinstance(self.components[i][j],MixtureModel):
                            self.components[i][j].modelInitialization(data.getInternalFeature(j),rtype=rtype,missing_value=missing_value)      
                        else:    
                            loc_l = l[i,:]
                            if missing_value:
                                # masking missing values from parameter estimation
                                for k,d in enumerate(data.getInternalFeature(j)):
                                    if d == missing_value:
                                        loc_l[k] = 0.0
                            self.components[i][j].MStep(loc_l,data.getInternalFeature(j))      

        # renormalizing mixing proportions in case of fixed components
        if fix_flag:
            if unfix_pi == 0.0:
                #print "----\n",self,"----\n"
                 print "unfix_pi = ", unfix_pi
                 print "fix_pi = ", fix_pi
                 print "pi = ", self.pi
                 raise RuntimeError, "unfix_pi = 0.0"  

            for i in range(self.G):
                if self.compFix[i] == 0:
                    self.pi[i] = (self.pi[i] * fix_pi) / unfix_pi    

    def classify(self, data, labels = None, entropy_cutoff = None, silent = 0):
        """
        Classification of input 'data'.
        Assignment to mixture components by maximum likelihood over
        the component membership posterior. No parameter reestimation.
        
        @param data: DataSet object
        @param labels: optional sample IDs
        @param entropy_cutoff: entropy threshold for the posterior distribution. Samples which fall
        above the threshold will remain unassigned
        @param silent: 0/1 flag, toggles verbose output        

        @return: list of class labels
        """        
        return  MixtureModel.classify(self,data,labels = labels, entropy_cutoff = entropy_cutoff, silent = silent, EStep=self.EStep, EStepParam = None )


class labeledBayesMixtureModel(BayesMixtureModel):
    """

    Bayesian mixture models with labeled data.

    """
    def __init__(self,G, pi, components, prior,compFix=None,struct=0):
       BayesMixtureModel.__init__(self,G, pi, components, prior,compFix=compFix,struct=struct,identifiable=0)

    
    def __copy__(self):
        copy_components = []
        
        #copy_pi = copy.copy(self.pi)
        copy_pi = copy.deepcopy(self.pi) 
        copy_compFix = copy.deepcopy(self.compFix)
        
        for i in range(self.G):
            copy_components.append(copy.deepcopy(self.components[i]))
        
        copy_prior = copy.copy(self.prior)

        copy_model = labeledBayesMixtureModel(self.G, copy_pi, copy_components, copy_prior, compFix = copy_compFix)
        copy_model.nr_tilt_steps = self.nr_tilt_steps    
        copy_model.suff_p = self.suff_p
        
        copy_model.identFlag = self.identFlag
        
        if self.struct:
            copy_model.initStructure()
            
            copy_leaders = copy.deepcopy(self.leaders)
            copy_groups = copy.deepcopy(self.groups)
        
            copy_model.leaders = copy_leaders
            copy_model.groups = copy_groups
        
        return copy_model
        
    
    def modelInitialization(self,data,rtype=1):
        """
        Perform model initialization given a random assigment of the
        data to the models.
            
        @param data: DataSet object
        @param rtype: type of random assignments.
        0 = fuzzy assingment
        1 = hard assingment
        @return posterior assigments
        """
        assert self.G >= len(data.labels), 'Insufficent number of components for given labeling.'

        if not isinstance(data,ConstrainedDataSet):
            raise TypeError, "DataSet object required, got"+ str(data.__class__)

        else:
            if data.internalData is None:
                data.internalInit(self)

        # reset structure if applicable
        if self.struct:
            self.initStructure()

        # generate 'random posteriors'
        l = numpy.zeros((self.G,len(data)),dtype='Float64')
        for i in range(len(data)):
            if rtype == 0:
              for j in range(self.G):
                 l[j,i] = random.uniform(0.1,1)

              s = sum(l[:,i])
              for j in range(self.G):
                  l[j,i] /= s
            else:
               l[random.randint(0,self.G-1),i] = 1 

        # peform label assigments (non random!)
        for i,cl in enumerate(data.labels): # for each class
            for o in cl: # for each observation in a class
              p_vec = numpy.zeros(self.G, dtype='Float64')
              p_vec[i] = 1.0
              l[:,o] = p_vec



        #print 'Random init l:'
        #for u in range(self.G):
        #    print u,l[u,:].tolist()

        # do one M Step
        fix_pi = 1.0
        unfix_pi = 0.0
        fix_flag = 0   # flag for fixed mixture components
        for i in range(self.G):
    
                # setting values for pi
                if self.compFix[i] == 2:
                    fix_pi -= self.pi[i]
                    fix_flag = 1

                else:
                    self.pi[i] =  l[i,:].sum() / len(data)                
                    unfix_pi += self.pi[i]
                    
                if self.compFix[i] == 1 or self.compFix[i] == 2:
                    # fixed component
                    continue
                else:    
                    # components are product distributions that may contain mixtures
                    if isinstance(self.components[i], ProductDistribution):
                        last_index = 0
                        for j in range(self.components[i].dist_nr):    
                            if isinstance(self.components[i].distList[j],MixtureModel):
                                dat_j = data.singleFeatureSubset(j)
                                self.components[i].distList[j].modelInitialization(dat_j,rtype=rtype)      
                            else:    
                                loc_l = l[i,:]
                                # masking missing values from parameter estimation
                                if data.missingSymbols.has_key(j):
                                    ind_miss = data.getMissingIndices(j)
                                    for k in ind_miss:
                                        loc_l[k] = 0.0
                                #self.components[i].distList[j].mapMStep(loc_l,data.getInternalFeature(j),self.prior.compPrior[j] )      
                                self.prior.compPrior.priorList[j].mapMStep(self.components[i].distList[j], loc_l,data.getInternalFeature(j) )  
                    
                    else:  # components are not ProductDistributions -> invalid    
                        raise TypeError

        # renormalizing mixing proportions in case of fixed components
        if fix_flag:
            if unfix_pi == 0.0:
                #print "----\n",self,"----\n"
                 print "unfix_pi = ", unfix_pi
                 print "fix_pi = ", fix_pi
                 print "pi = ", self.pi
                 raise RuntimeError, "unfix_pi = 0.0"  

            for i in range(self.G):
                if self.compFix[i] == 0:
                    self.pi[i] = (self.pi[i] * fix_pi) / unfix_pi
        return l



    def mapEM(self,data, max_iter, delta,silent = False, mix_pi=None, mix_posterior= None, tilt = 0):
        """
        Reestimation of maximum a posteriori mixture parameters using the EM algorithm.
        
        @param data: DataSet object
        @param max_iter: maximum number of iterations
        @param delta: minimal difference in likelihood between two iterations before
        convergence is assumed.
        @param silent: 0/1 flag, toggles verbose output
        @param mix_pi: [internal use only] necessary for the reestimation of
        mixtures as components
        @param mix_posterior:[internal use only] necessary for the reestimation of
        mixtures as components
        @param tilt: 0/1 flag, toggles the use of a deterministic annealing in the training
        
        @return: tuple of posterior matrix and log-likelihood from the last iteration
        """
        assert self.G >= len(data.labels), 'Insufficent number of components for given labeling.'

        if isinstance(data, numpy.ndarray):
            raise TypeError, "DataSet object required."
        elif  isinstance(data, DataSet):
            if data.internalData is None:
                if not silent:
                    sys.stdout.write("Parsing data set...")
                    sys.stdout.flush()
                data.internalInit(self) 
                if not silent:
                    sys.stdout.write("done\n")
                    sys.stdout.flush()
        else:    
            raise ValueError, "Unknown input type format: " + str(data.__class__)

        log_p_old = float('-inf')
        step = 0

        # if deterministic annealing is activated, increase number of steps by self.nr_tilt_steps 
        if tilt:
            if not silent:
                sys.stdout.write("Running EM with "+ str(self.nr_tilt_steps) +" steps of deterministic annealing.\n" )
            max_iter += self.nr_tilt_steps
            
        # for lower hierarchy mixture we need the log of mix_posterior
        if mix_posterior is not None:
            log_mix_posterior = numpy.log(mix_posterior)

        while 1:
            log_p = 0.0
            # matrix of log posterior probs: components# * (sequence positions)
            log_l = numpy.zeros((self.G,data.N),dtype='Float64')
            #log_col_sum = numpy.zeros(data.N,dtype='Float64')  # array of column sums of log_l
            log_pi = numpy.log(self.pi)  # array of log mixture coefficients
            
            # computing log posterior distribution 
            for i in range(self.G):            
                log_l[i] = log_pi[i] + self.components[i].pdf(data) 

    	    # Partially supervised training
            # For sequences with partial information use a weight vector s.t.
            # P[model|seq] = 1 if model = label[seq] and 0 else
            for i,cl in enumerate(data.labels): # for each class
                for o in cl: # for each observation in a class
                  v = log_l[i,o]
                  p_vec = numpy.zeros(self.G, dtype='Float64')
                  p_vec[:] = float('-inf')
                  p_vec[i] = v
                  log_l[:,o] = p_vec

            
            # computing data log likelihood as criteria of convergence
            # log_l is normalized in-place and likelihood is returned as log_p
            log_p = _C_mixextend.get_normalized_posterior_matrix(log_l)

            # adjusting posterior for lower hierarchy mixtures
            if mix_posterior is not None:
                # multiplying in the posterior of upper hierarchy mixture
                log_l[:,j] = log_l[:,j] + log_mix_posterior[j]

            
            # we have to take the parameter prior into account to form the objective function
            # If we assume independence between parameters in different components, the prior
            # contribution is given by a product over the individual component and structure priors
            try:
                log_prior = self.prior.pdf(self)
            except ValueError:  # catch zero probability under prior
                
                raise ConvergenceFailureEM,"Zero probability under prior." 
    
            # calculate objective function
            log_p += log_prior
            # checking for convergence 
            diff = (log_p - log_p_old) 

            if log_p_old != -1.0 and not silent and step > 0:
                if tilt and step <= self.nr_tilt_steps:
                    sys.stdout.write("TILT Step "+str(step)+": log posterior: "+str(log_p)+"\n")
                else:
                    sys.stdout.write("Step "+str(step)+": log posterior: "+str(log_p)+ "   (diff="+str(diff)+")\n")

            if diff < 0.0 and step > 1 and abs(diff / log_p_old) > self.err_tol:
                #print log_p,log_p_old, diff,step,abs(diff / log_p_old)
                #print "WARNING: EM divergent."
                raise ConvergenceFailureEM,"Convergence failed, EM divergent: "  
                                
            if (not tilt or (tilt and step+1 >= self.nr_tilt_steps)) and delta >= diff and max_iter != 1: 
                if not silent:
                    sys.stdout.write("Convergence reached with log_p "+ str(log_p)+ " after "+str(step)+" steps.\n")
                if self.identFlag:
                    self.identifiable()
                return (log_l,log_p)  

            log_p_old = log_p
            if step == max_iter:
                if not silent:
                    sys.stdout.write("Max_iter "+str(max_iter)+" reached -> stopping\n")

                if self.identFlag:
                    self.identifiable()
                return (log_l,log_p) 

           
            # compute posterior likelihood matrix from log posterior
            l = numpy.exp(log_l)

            # deterministic annealing, shifting posterior toward uniform distribution.
            if tilt and step+1 <= self.nr_tilt_steps and mix_posterior is None:
                h = self.heat - (step * (self.heat/(self.nr_tilt_steps))  )
                for j in range(data.N):
                    uni = 1.0 / self.G
                    tilt_l = (uni - l[:,j]) * h
                    l[:,j] += tilt_l

            # variables for component fixing
            fix_pi = 1.0
            unfix_pi = 0.0
            fix_flag = 0   # flag for fixed mixture components
            
            # update component parameters and mixture weights
            for i in range(self.G):
                if self.compFix[i] & 2:   # pi[i] is fixed
                    fix_pi -= self.pi[i]
                    fix_flag = 1
                    continue
                else:
                    # for mixtures of mixtures we need to multiply in the mix_pi[i]s                 
                    if mix_pi is not None:                                                  
                        self.pi[i] =  ( l[i,:].sum() + self.prior.piPrior.alpha[i] -1.0 ) / ((data.N * mix_pi) + self.prior.piPrior.alpha_sum - self.G ) 
                    else:
                        self.pi[i] =  ( l[i,:].sum() + self.prior.piPrior.alpha[i] -1.0 ) / (data.N + ( self.prior.piPrior.alpha_sum - self.G) ) 
                    
                    unfix_pi += self.pi[i]

                if self.compFix[i] & 1:
                    continue                
                else:
                    # Check for model structure
                    if not self.struct:
                        self.prior.compPrior.mapMStep(self.components[i],l[i,:],data,self.pi[i])
                    
            # if there is a model structure we update the leader distributions only 
            if self.struct:
                for j in range(self.dist_nr):
                    for k in self.leaders[j]:
                        # compute group posterior
                        # XXX extension function for pooled posterior ?
                        g_post = copy.deepcopy(l[k,:])  
                        g_pi = self.pi[k]
                        for memb in self.groups[j][k]:
                            g_post += l[memb,:]
                            g_pi += self.pi[memb]
                        
                        if isinstance(self.components[k][j],MixtureModel): 
                            self.prior.compPrior[j].mapMStep(self.components[k][j],g_post,data.singleFeatureSubset(j), g_pi)
                        else:    
                            try:
                                self.prior.compPrior[j].mapMStep(self.components[k][j],g_post,data.getInternalFeature(j) )
                            except InvalidPosteriorDistribution:
                                raise    
            
            # renormalizing mixing proportions in case of fixed components
            if fix_flag:
                if unfix_pi == 0.0:
                    #print "----\n",self,"----\n"
                    print "unfix_pi = ", unfix_pi
                    print "fix_pi = ", fix_pi
                    print "pi = ", self.pi
                    raise RuntimeError, "unfix_pi = 0.0"  
                for i in range(self.G):
                    if self.compFix[i] == 0:
                        self.pi[i] = (self.pi[i] * fix_pi) / unfix_pi

            sys.stdout.flush()
            step += 1              


    def updateStructureBayesian(self,data,objFunction='MAP',silent=1):
        """
        Updating structure by chosing optimal local merge with respect to the posterior.

        
        Features: - store merges in a learning history to prevent recomputation of merge parameters
                  - compute parameters of candidate structures from expected sufficient statistics of groups to be merged
        
       
        @param data: DataSet object
        @param silent: verbosity flag

        @return: number of structure changes 
        """
        assert self.struct == 1, "No structure in model."
        assert objFunction in ['MAP'] # for now only MAP estimation

        new_leaders = []
        new_groups = []
        change = 0
        # building data likelihood factor matrix for the current group structure      
        l = numpy.zeros( (self.dist_nr, self.G, data.N),dtype='Float64' )
        for j in range(self.dist_nr):
            # extracting current feature from the DataSet
            if isinstance(self.components[0][j], MixtureModel): # XXX
                data_j = data.singleFeatureSubset(j)
            else:
                data_j = data.getInternalFeature(j)


            for lead_j in self.leaders[j]:
                l_row = self.components[lead_j][j].pdf(data_j )  
                l[j,lead_j,:] = l_row
                for v in self.groups[j][lead_j]:
                    l[j,v,:] = l_row

        # apply labels to posterior matrix 
        # XXX there should be a more efficient way to do this ... XXX 
        for i,cl in enumerate(data.labels): # for each class
            for o in cl: # for each observation in a class
                ind = range(self.G)
                ind.pop(i)
                for ng in ind:
                    l[:,ng,o] = float('-inf')

        # g is the matrix of log posterior probabilities of the components given the data
        g = numpy.sum(l, axis=0) 
        for k in range(self.G):
            g[k,:] += numpy.log(self.pi[k])

        sum_logs = matrixSumlogs(g)
        g_norm = g - sum_logs

        tau = numpy.exp(g_norm)
        if not silent:
            print "\ntau="
            for tt in tau:
                print tt.tolist()
            print

        # computing posterior as model selection criterion
        temp = DiscreteDistribution(self.G,self.pi)
        pi_prior = self.prior.piPrior.pdf(temp)
        log_prior = pi_prior
        log_prior_list = [0.0] * self.dist_nr
        for j in range(self.dist_nr):
            for r in range(self.G):
                   log_prior_list[j] += self.prior.compPrior[j].pdf( self.components[r][j] )
        log_prior += sum(log_prior_list)

        # prior over number of components
        log_prior += self.prior.nrCompPrior * self.G 
        # prior over number of distinct groups
        for j in range(self.dist_nr):
            log_prior += self.prior.structPrior * len(self.leaders[j])   

        # get posterior
        lk = numpy.sum(sum_logs) 
        post = lk + log_prior
        if not silent:
            print "0: ",  lk ,"+", log_prior,"=", post
            print log_prior_list

        changes = 0        
        g_wo_j = numpy.zeros((self.G, data.N),dtype='Float64')

        # initialising temporary group structure with copies of the current structure
        temp_leaders = copy.deepcopy(self.leaders)
        temp_groups =  copy.deepcopy(self.groups)
        for j in range(self.dist_nr):
            L = {}  # initialize merge history

            if not silent:
                print "\n************* j = ",j,"*****************\n"

            # unnormalized posterior matrix without the contribution of the jth feature
            try:
                g_wo_j  = g - l[j]
            except FloatingPointError:
                # if there was an exception we have to compute each
                # entry in g_wo_j seperately to set -inf - -inf = -inf
                g_wo_j = _C_mixextend.substract_matrix(g,l[j])

            # checking whether feature j is already fully merged 
            nr_lead = len(self.leaders[j])
            if nr_lead == 1:
                continue  # nothing to be done...

            term = 0            

            # extracting current feature from the DataSet
            if isinstance(self.components[0][j], MixtureModel): # XXX
                data_j = data.singleFeatureSubset(j)
            else:
                data_j = data.getInternalFeature(j)

            tau_pool = numpy.zeros(data.N, dtype='Float64')
            for lead in self.leaders[j]:
                el_dist = copy.copy(self.components[lead][j])

                tau_pool = copy.copy(tau[lead, :])
                pi_pool = self.pi[lead]
                for z in self.groups[j][lead]:
                    tau_pool += tau[z, :]                    
                    pi_pool += self.pi[z]
                
                stat = el_dist.sufficientStatistics(tau_pool, data_j)
                M = CandidateGroup(el_dist, numpy.sum(tau_pool), pi_pool, stat)
                L[(lead,)+tuple(self.groups[j][lead])] = M
                
            while not term:
                best_dist = None   # accepted candidate distributions
                best_post = float('-inf')   # corresponding posteriors
                best_indices = None
                best_l_j = l[j] 
                best_log_prior_list_j = log_prior_list[j]                
                
                for mc1 in range(len(temp_leaders[j])):
                    merge_cand1 = temp_leaders[j][mc1]
                    for mc2 in range(mc1+1,len(temp_leaders[j])):
                        merge_cand2 = temp_leaders[j][mc2]
                        if not silent:
                            print "-------------------"
                            print merge_cand1," -> ",merge_cand2
                            print self.components[merge_cand1][j]
                            print self.components[merge_cand2][j]

                        nr_leaders_j = len(temp_leaders[j])-1
                        cand_group_j = temp_groups[j][merge_cand1] + [merge_cand2] + temp_groups[j][merge_cand2]
                        
                        hist_ind_part1 = (merge_cand1,)+tuple(temp_groups[j][merge_cand1])
                        hist_ind_part2 = (merge_cand2,)+tuple(temp_groups[j][merge_cand2])
                        hist_ind_complete = hist_ind_part1 + hist_ind_part2


                        recomp = 0
                        if L.has_key(hist_ind_complete):
                            recomp = 1
                        if not silent:
                            print "\ncandidate model structure: XXX"
                            #print "lead = ",leaders_j
                            #print "groups = ",groups_j
                            #print "others = ",others,"\n"

                        if not recomp:
                            assert L.has_key( hist_ind_part1),str(hist_ind_part1)+' missing.'
                            assert L.has_key( hist_ind_part2),str(hist_ind_part2)+' missing.'
                            
                            M = self.prior.compPrior[j].mapMStepMerge([L[hist_ind_part1], L[hist_ind_part2]])

                            candidate_dist = M.dist

                            if not silent:
                                print "candidate:", candidate_dist

                            l_row = candidate_dist.pdf(data_j)  
                            cdist_prior = self.prior.compPrior[j].pdf( candidate_dist ) 
                           
                            M.l = l_row
                            M.dist_prior = cdist_prior
                            L[hist_ind_complete] = M
                        else:
                            # retrieve merge data from history
                            candidate_dist = L[hist_ind_complete].dist                            
                            if not silent:
                                print "candidate:", candidate_dist

                            l_row = L[hist_ind_complete].l
                            cdist_prior = L[hist_ind_complete].dist_prior

                	    # computing change in likelihood matrix for this step
                        l_j_1 = copy.copy(l[j]) 

                        # updating l_j_1 with the new candidate distribution
                        l_j_1[merge_cand1, :] = l_row
                        for v in cand_group_j:
                            l_j_1[v,:] = l_row

                        # applying labels 
                        for i,cl in enumerate(data.labels): # for each class
                            for o in cl: # for each observation in a class
                                ind = range(self.G)
                                ind.pop(i)
                                for ng in ind:
                                    l_j_1[ng,o] = float('-inf')

                        # get updated unnormalized posterior matrix
                        g = g_wo_j + l_j_1

                        sum_logs = matrixSumlogs(g)
                        lk_1 = numpy.sum(sum_logs)

                        # computing posterior as model selection criterion
                        log_prior_1 = pi_prior

                        # compute parameter prior for the candidate merge parameters
                        log_prior_list_j = 0.0
                        for r in range(self.G):
                            if r in [merge_cand1] + cand_group_j:
                                log_prior_list_j += cdist_prior
                            else:
                                log_prior_list_j += self.prior.compPrior[j].pdf( self.components[r][j] )

                        log_prior_1 += sum(log_prior_list)
                        log_prior_1 -= log_prior_list[j]
                        log_prior_1 += log_prior_list_j

                        # prior over number of components
                        log_prior_1 += self.prior.nrCompPrior * self.G 
                        # prior over number of distinct groups
                        for z in range(self.components[0].dist_nr):
                            if z == j:
                                log_prior_1 += self.prior.structPrior * nr_leaders_j  
                            else:
                                log_prior_1 += self.prior.structPrior * len(temp_leaders[z])   # XXX len could be cached ?

                        post_1 = lk_1 + log_prior_1

                        if not silent:
                            print 'Posterior:',post_1 ,'=', lk_1 ,'+', log_prior_1

                        if post_1 >= post:
                            if not silent:
                                print "*** Merge accepted", post_1 ,">=", post

                            if post_1 > best_post:  # current merge is better than previous best
                                best_dist = candidate_dist
                                best_post = post_1
                                best_indices = [merge_cand1, merge_cand2]
                                best_l_j = l_j_1
                                best_log_prior_list_j = log_prior_list_j
                        else:
                            if not silent:
                                print "*** Merge rejected:", post_1 ,"!>", post

                # if there is no possible merge that increases the score we are done
                if best_post == float('-inf'):
                    if not silent:
                        print "*** Finished !"
                    # setting updated structure in model
                    self.leaders[j] = temp_leaders[j]
                    self.groups[j] = temp_groups[j]

                    # reset posterior matrix to the last accepted merge
                    g = g_wo_j + best_l_j
                    term = 1
                # otherwise we update the model with the best merge found
                else:  
                    if not silent:
                        print "--- Winner ---"
                        print "post:", best_post
                        print "indices:",best_indices
                        print "dist: ",best_dist
                        #print "lead:",best_leaders
                        #print "group:", best_groups

                    post = best_post
                    l[j] = best_l_j
                    g = g_wo_j + best_l_j  # posterior matrix for the next iteration
                    log_prior_list[j] = best_log_prior_list_j

                    # updating model
                    # removing merged leader from new_leaders
                    ind = temp_leaders[j].index(best_indices[1])
                    temp_leaders[j].pop(ind)

                    # joining merged groups and removing old group entry
                    temp_groups[j][best_indices[0]] += [best_indices[1]] + temp_groups[j][best_indices[1]]
                    temp_groups[j].pop(best_indices[1]) 

                    # assigning distributions according to new structure
                    self.components[best_indices[0]][j] = best_dist
                    for d in temp_groups[j][best_indices[0]]:
                        self.components[d][j] = self.components[best_indices[0]][j]
                    change += 1
        return change 


    def classify(self, data, labels = None, entropy_cutoff = None, silent = 0, EStep=None, EStepParam = None ):
        """
        Classification of input 'data'.
        Assignment to mixture components by maximum likelihood over
        the component membership posterior. No parameter reestimation.
        
        Classification of labelled samples is fixed a priori and overrrides the assignment by
        maximum posterior.
        
        @param data: ConstrainedDataSet object
        @param labels: optional sample IDs
        @param entropy_cutoff: entropy threshold for the posterior distribution. Samples which fall
        above the threshold will remain unassigned
        @param silent: 0/1 flag, toggles verbose output        
        @param EStep: function implementing the EStep, by default self.EStep
        @param EStepParam: additional paramenters for more complex EStep implementations
     
        
        @return: list of class labels
        """
        
        assert isinstance(data,ConstrainedDataSet)
        
        # standard classification 
        c = BayesMixtureModel.classify(self,data,labels=labels,entropy_cutoff=entropy_cutoff,silent=1,EStep=EStep,EStepParam=EStepParam)
        
        # apply labels
        for i,cl in enumerate(data.labels): # for each class
            for o in cl: # for each observation in a class
                c[o] = i
        
        if not silent:
            # printing out the clusters
            cluster = {}
            en = {}
            for j in range(-1,self.G,1):  
                cluster[j]=[]

            for i in range(data.N):
                cluster[c[i]].append(data.sampleIDs[i])
        
            print "\n** Clustering **"
            for j in range(self.G):
                print "Cluster ",j,', size',len(cluster[j])
                print cluster[j], "\n" 
        
            print "Unassigend due to entropy cutoff:"
            print cluster[-1], "\n"


        return c

def CMMfromMM(mm):
    """
    Convenience function. Takes a MixtureModel or and returns a ConstrainedMixtureModel with
    the same parameters.
    
    @param mm: MixtureModel object
    """
    return ConstrainedMixtureModel(mm.G,mm.pi,mm.components,mm.compFix,mm.struct)
    

class ConstrainedMixtureModel(MixtureModel):
    """
    Class for a mixture model containing the pairwise constrained version of the E-Step
    """
    def __init__(self,G, pi, components,compFix=None,struct=0):
       MixtureModel.__init__(self,G, pi, components, compFix=compFix, struct=struct,identifiable=0)

    def __copy__(self):
        copy_model = MixtureModel.__copy__(self)
        return CMMfromMM(copy_model)
    
    def EM(self, data, max_iter, delta, prior_positive,
           prior_negative, previous_posterior, prior_type, normaliziering=False, silent = False,
           mix_pi=None, mix_posterior= None, tilt = 0):
        """ 
        Reestimation of mixture parameters using the EM algorithm.
        This method do some initial checking and call the EM from
        MixtureModel with the constrained labels E step
                
        @param data: DataSet object
        @param max_iter: maximum number of iterations
        @param delta: minimal difference in likelihood between two iterations before
        convergence is assumed.
        @param prior_positive: importance parameter for positive constraints
        @param prior_negative: importance parameter for negative constraints
        @param previous_posterior: matrix containing the posterior of the previous model assigments
        @param prior_type: 1 positive constr.
                           2 negative constr.
                           3 positive and negative constr.
        @param silent: 0/1 flag, toggles verbose output
        @param mix_pi: [internal use only] necessary for the reestimation of
        mixtures as components
        @param mix_posterior:[internal use only] necessary for the reestimation of
        mixtures as components
        @param tilt: 0/1 flag, toggles the use of a deterministic annealing in the training
        
        @return: tuple of posterior matrix and log-likelihood from the last iteration
        """
        assert  isinstance(data,ConstrainedDataSet), 'Data set does not contain labels, Labeled EM can not be performed'
        if data.pairwisepositive == None:
          assert (prior_positive == 0), 'Data set does not contain pairwise constraints, Labeled EM can not be performed'
        if data.pairwisenegative == None:
          assert (prior_negative == 0), 'Data set does not contain pairwise constraints, Labeled EM can not be performed'

        class EArgs:
             def __init__(self):
                 self.prior_positive = prior_positive
                 self.prior_negative = prior_negative
                 self.normaliziering = normaliziering
                 self.previous_posterior = previous_posterior
                 self.prior_type = prior_type

        return MixtureModel.EM(self,data, max_iter, delta,silent = silent, mix_pi=mix_pi, mix_posterior= mix_posterior,
                               tilt = tilt, EStep=self.EStep, EStepParam = EArgs())        
        

    def modelInitialization(self,data,prior_positive, prior_negative, prior_type, normaliziering=False,rtype=1):
        """
        Perform model initialization given a random assigment of the
        data to the models.
            
        @param data: DataSet object
        @param rtype: type of random assignments.
        0 = fuzzy assingment
        1 = hard assingment
        @return posterior assigments
        """

        if not isinstance(data,ConstrainedDataSet):
            raise TypeError, "DataSet object required, got"+ str(data.__class__)

        else:
            if data.internalData is None:
                data.internalInit(self)

        # reset structure if applicable
        if self.struct:
            self.initStructure()

        # generate 'random posteriors'
        l = numpy.zeros((self.G,len(data)),dtype='Float64')
        for i in range(len(data)):
            if rtype == 0:
              for j in range(self.G):
                 l[j,i] = random.uniform(0.1,1)

              s = sum(l[:,i])
              for j in range(self.G):
                  l[j,i] /= s
            else:
               l[random.randint(0,self.G-1),i] = 1 

        class EArgs:
             def __init__(self):
                 self.prior_positive = prior_positive
                 self.prior_negative = prior_negative
                 self.normaliziering = normaliziering
                 self.previous_posterior = l
                 self.prior_type = prior_type


        # peform constrain updates (non random!)
        self.EStep(data, EStepParam=EArgs())

        #print 'Random init l:'
        #for u in range(self.G):
        #    print u,l[u,:].tolist()

        # do one M Step
        fix_pi = 1.0
        unfix_pi = 0.0
        fix_flag = 0   # flag for fixed mixture components
        for i in range(self.G):
    
                # setting values for pi
                if self.compFix[i] == 2:
                    fix_pi -= self.pi[i]
                    fix_flag = 1

                else:
                    self.pi[i] =  l[i,:].sum() / len(data)                
                    unfix_pi += self.pi[i]
                    
                if self.compFix[i] == 1 or self.compFix[i] == 2:
                    # fixed component
                    continue
                else:    
                    # components are product distributions that may contain mixtures
                    if isinstance(self.components[i], ProductDistribution):
                        last_index = 0
                        for j in range(self.components[i].dist_nr):    
                            if isinstance(self.components[i].distList[j],MixtureModel):
                                dat_j = data.singleFeatureSubset(j)
                                self.components[i].distList[j].modelInitialization(dat_j,rtype=rtype)      
                            else:    
                                loc_l = l[i,:]
                                # masking missing values from parameter estimation
                                if data.missingSymbols.has_key(j):
                                    ind_miss = data.getMissingIndices(j)
                                    for k in ind_miss:
                                        loc_l[k] = 0.0
                                self.components[i][j].MStep(loc_l,data.getInternalFeature(j))      
                    
                    else:  # components are not ProductDistributions -> invalid    
                        raise TypeError

        # renormalizing mixing proportions in case of fixed components
        if fix_flag:
            if unfix_pi == 0.0:
                #print "----\n",self,"----\n"
                 print "unfix_pi = ", unfix_pi
                 print "fix_pi = ", fix_pi
                 print "pi = ", self.pi
                 raise RuntimeError, "unfix_pi = 0.0"  

            for i in range(self.G):
                if self.compFix[i] == 0:
                    self.pi[i] = (self.pi[i] * fix_pi) / unfix_pi
        return l


    def classify(self, data, prior_positive, prior_negative, previous_posterior,
                 prior_type, labels = None, entropy_cutoff = None,
                 silent = 0,  normaliziering=False):
        """
        Classification of input 'data'.  Assignment to mixture components by maximum likelihood
        over the component membership posterior. No parameter
        reestimation.
        
        @param data: DataSet object
        @param labels: optional sample IDs
        @param prior_positive: importance parameter for positive constraints
        @param prior_negative: importance parameter for negative constraints
        @param previous_posterior: matrix containing the posterior of the previous model assigments
        @param prior_type: 1 positive constr.
                           2 negative constr.
                           3 positive and negative constr.
        @param entropy_cutoff: entropy threshold for the posterior distribution. Samples which fall
        above the threshold will remain unassigned
        @param silent: 0/1 flag, toggles verbose output        

        @return: list of class labels
        """
        
        assert  isinstance(data,ConstrainedDataSet), 'Data set does not contain labels, Labeled EM can not be performed'
        if data.pairwisepositive == None:
          assert (prior_positive == 0), 'Data set does not contain pairwise constraints, Labeled EM can not be performed'
        if data.pairwisenegative == None:
          assert (prior_negative == 0), 'Data set does not contain pairwise constraints, Labeled EM can not be performed'

        class EArgs:
             def __init__(self):
                 self.prior_positive = prior_positive
                 self.prior_negative = prior_negative
                 self.normaliziering = normaliziering
                 self.previous_posterior = previous_posterior
                 self.prior_type = prior_type
                 
        return  MixtureModel.classify(self,data,labels = labels, entropy_cutoff = entropy_cutoff, 
                                      silent = silent, EStep=self.EStep, EStepParam = EArgs() )

    def EStep(self,data,mix_posterior=None,mix_pi=None,EStepParam=None):
        """
        Reestimation of mixture parameters using the EM algorithm.
                
        @param data: DataSet object
        @param mix_pi: [internal use only] necessary for the reestimation of
        mixtures as components
        @param mix_posterior:[internal use only] necessary for the reestimation of
        mixtures as components
        @param EStepParam: additional paramenters for more complex EStep implementations, in
        this implementaion it is ignored

        @return: tuple of log likelihood matrices and sum of log-likelihood of components
        """
        prior_pos = EStepParam.prior_positive
        prior_neg = EStepParam.prior_negative
        norm = EStepParam.normaliziering
        prior_type = EStepParam.prior_type
        # HACK - I am updating at each iteration the previous posterior,
        # since this parameter is not part of of regular EStep
        previous_posterior = EStepParam.previous_posterior
        positive_constraints = data.pairwisepositive
        negative_constraints = data.pairwisenegative

        log_l = numpy.zeros((self.G,data.N),dtype='Float64')
        log_col_sum_nopen =  numpy.zeros(data.N,dtype='Float64')  # array of column sums of log_l without penalty
        log_col_sum = numpy.zeros(data.N,dtype='Float64')  # array of column sums of log_l
        log_pi = numpy.log(self.pi)  # array of log mixture coefficient

        # computing log posterior distribution 
        for i in range(self.G):            
            log_l[i] = log_pi[i] + self.components[i].pdf(data)

        log_col_sum_nopen[i] = sumlogs(log_l[:,i]) # sum over jth column of log_l without penalization
      
        #changing order of indices assigments
        indices = range(data.N)

        random.shuffle(indices)
        
        pen = numpy.zeros(self.G,dtype='Float64')
        penn = numpy.zeros(self.G,dtype='Float64')     
        
        
        
        for x,i in enumerate(indices):

#          print '\n----- sample'+str(i)
#          print 'log_l=',log_l[:,i]
#          print 'norm l', numpy.exp( log_l[:,i] - sumlogs(log_l[:,i]) ) 

          # Leaves -Inf values unchanged
          pen[:] = 0.0
          penn[:] = 0.0 
          for y,j in enumerate(indices):
              
              #print x,y,'->',i,j
              
              # calculating penalities
              # in a Gibbs sampling manner (using either previous or current posteriors
              if y > x: # if posterior not yet calculated, use of previous one posterior
                coc =  numpy.multiply(previous_posterior[:,j],1.0) # posterior of y
                
                
                if prior_type == 1 or prior_type == 3:
                  if positive_constraints[i][j] > 0.0:
                    if norm:                    
                      pen += numpy.divide(numpy.multiply(1-coc, positive_constraints[i][j]),((1-self.pi)))
                    else:
                      pen += numpy.multiply(1-coc, positive_constraints[i][j])

#                      print '   +  1 - coc',j, 1 - coc
#                      print '   +',pen                      


                if prior_type == 2 or prior_type == 3:
                  if negative_constraints[i][j] > 0.0:
                    if norm:
                      penn += numpy.divide(numpy.multiply(coc, negative_constraints[i][j]),self.pi)
                    else:
                      penn += numpy.multiply(coc, negative_constraints[i][j])

#                      print '   - coc',j, coc
#                      print '   -',penn                      


              elif y < x:                      
                coc =  numpy.multiply(numpy.exp(log_l[:,j]),1)
                
                
                
                if prior_type == 1 or prior_type == 3:
                  if positive_constraints[i][j] > 0.0:
                    if norm:
                      pen += numpy.divide(numpy.multiply(1-coc, positive_constraints[i][j]),((1-self.pi)))
                    else:
                      pen += numpy.multiply(1-coc, positive_constraints[i][j])
                
#                      print '   + coc',j, coc
#                      print '   +',pen                      
                      
                if prior_type == 2 or prior_type == 3:
                    if negative_constraints[i][j] > 0.0:
                      if norm:
                        penn += numpy.divide(numpy.multiply(coc, negative_constraints[i][j]),self.pi)
                      else:
                        penn += numpy.multiply(coc, negative_constraints[i][j])
   
#                        print '   - coc',j, coc
#                        print '   -',penn                      

          #print '\n-------', i
          #print log_l[:,i]
#          print - numpy.multiply(pen,prior_pos)
#          print - numpy.multiply(penn,prior_neg)
        

              
          log_l[:,i] += (-numpy.multiply(pen,prior_pos) - numpy.multiply(penn,prior_neg))
          # l[k,i] = log( (a_k * * P[seq i| model k]) + P[W+|y] * P[W-|y] )

#          print '-> log_l=',log_l[:,i]
#          print '-> norm l=',numpy.exp( log_l[:,i] - sumlogs(log_l[:,i]) ) 

          #print '->',log_l[:,i]

          log_col_sum[i] = sumlogs(log_l[:,i]) # sum over jth column of log_l    
          # if posterior is invalid, check for model validity
          if log_col_sum[i] == float('-inf'): 

                # if self is at the top of hierarchy, the model is unable to produce the
                # sequence and an exception is raised. Otherwise normalization is not necessary.
                if mix_posterior is None and not mix_pi:
                    print "\n---- Invalid -----\n",self,"\n----------"
                    #print "\n---------- Invalid ---------------"
                    print "mix_pi = ", mix_pi
                    #print "x[",i,"] = ", sequence[j]
                    print "l[:,",i,"] = ", log_l[:,i] 
                    raise InvalidPosteriorDistribution, "Invalid posterior distribution."
          # for valid posterior, normalize and go on    
          else:
                # normalizing log posterior
                log_l[:,i] = log_l[:,i] - log_col_sum[i]

                # adjusting posterior for lower hierarchy mixtures
                if mix_posterior is not None:
                    log_mix_posterior = numpy.zeros(len(mix_posterior),dtype='Float64')
                    for k in range(len(mix_posterior)):
                        if mix_posterior[k] == 0.0:
                            log_mix_posterior[k] = float('-inf')
                        else:
                            log_mix_posterior[k] = numpy.log(mix_posterior[k])

                    # multiplying in the posterior of upper hierarch mixture
                    log_l[:,i] = log_l[:,i] + log_mix_posterior[i]

        # final penalty (P(W|Y))
        penalty = 0.0
        for x,i in enumerate(indices):
            for y,j in enumerate(indices):
                coc =  numpy.multiply(numpy.exp(log_l[:,j]),1)
                if prior_type == 1 or prior_type == 3:
                  if positive_constraints[i][j] > 0.0:
                    if norm:
                      penalty += numpy.sum(numpy.divide(numpy.multiply(1-coc, positive_constraints[i][j]),((1-self.pi))))
                    else:
                      penalty += numpy.sum(numpy.multiply(1-coc, positive_constraints[i][j]))
                if prior_type == 2 or prior_type == 3:
                    if negative_constraints[i][j] > 0.0:
                      if norm:
                        penalty += numpy.sum(numpy.divide(numpy.multiply(coc, negative_constraints[i][j]),self.pi))
                      else:
                        penalty += numpy.sum(numpy.multiply(coc, negative_constraints[i][j]))
                    
        # HACK - I am updating at each iteration the previous posterior,
        #since this parameter is not part of of regular EStep          
        EStepParam.previous_posterior = numpy.exp(log_l)
        # computing data log likelihood as criteria of convergence        
        log_p = numpy.sum(log_col_sum) + penalty
        return log_l, numpy.sum(log_col_sum) + penalty
       




    def randMaxEM(self,data,nr_runs,nr_steps, delta, prior_positive,
           prior_negative, prior_type, tilt=0,silent=False):
        """
        Performs `nr_runs` normal EM runs with random initial parameters
        and returns the model which yields the maximum likelihood.
        
        @param data: DataSet object
        @param nr_runs: number of repeated random initializations
        @param nr_steps: maximum number of steps in each run
        @param delta: minimim difference in log-likelihood before convergence
        @param tilt: 0/1 flag, toggles the use of a deterministic annealing in the training
        @param silent:0/1 flag, toggles verbose output
        
        @return: log-likelihood of winning model
        """
        if isinstance(data, numpy.ndarray):
            raise TypeError, "DataSet object required."
        elif  isinstance(data, DataSet):
            if data.internalData is None:
                if not silent:
                    sys.stdout.write("Parsing data set...")
                    sys.stdout.flush()
                data.internalInit(self)
                if not silent:
                    sys.stdout.write("done\n")
                    sys.stdout.flush()
        else:    
            raise ValueError, "Unknown input type format: " + str(data.__class__)

        logp_list = []
        best_logp = float('-inf')
        best_model = None
        candidate_model = copy.copy(self)  # copying the model parameters 

        for i in range(nr_runs):
            # we do repeated random intializations until we get a model with valid posteriors in all components
            init = 0
            while not init:
                try:        
                    previous_posterior = candidate_model.modelInitialization(data,prior_positive, prior_negative, prior_type,)    # randomizing parameters of the model copy
                except InvalidPosteriorDistribution:
                    pass
                else:
                    init = 1     


            
            try:        
#                                 def EM(self, data, max_iter, delta, prior_positive, prior_negative, previous_posterior, prior_type, normaliziering=False, silent = False,  mix_pi=None, mix_posterior= None, tilt = 0):
                (l,log_p) = candidate_model.EM(data,nr_steps, delta, prior_positive, prior_negative, previous_posterior, prior_type, silent=silent)  # running EM 
            except ConvergenceFailureEM:
                sys.stdout.write("Run "+str(i)+" produced invalid model, omitted.\n")
            except InvalidPosteriorDistribution:
                sys.stdout.write("Run "+str(i)+" produced invalid model, omitted.\n")
            else:
                logp_list.append(log_p)

                # check whether current model is better than previous optimum
                if log_p > best_logp:
                    best_model = copy.copy(candidate_model)
                    best_logp = log_p
                    best_l = copy.copy(l)
        if not silent:
            print "\nBest model likelihood over ",nr_runs,"random initializations:"
            print "Model likelihoods:",logp_list
            print "Average logp: ", sum(logp_list)/float(nr_runs)," SD:",numpy.array(logp_list).std()
            print "Best logp:",best_logp

        # check whether at least one run was sucessfully completed        
        if best_model == None:
            raise ConvergenceFailureEM, 'All '+ str(nr_runs)+' runs have failed.'
        
        self.components = best_model.components  # assign best parameter set to model 'self'
        self.pi = best_model.pi

        return best_l, best_logp  # return final data log likelihood



#---------------------------------- Miscellaneous -------------------------------------------------        

        
def structureAccuracy(true,m):
    """
    Returns the accuracy of two model structures with respect to
    the component partition they define.
    
    @param true: MixtureModel object with CSI structure
    @param m: MixtureModel object with CSI structure

    @return: agreement of the two structures as measure by the accuracy
    """

    tp = fn = tn = fp = 0
    for j in range(true.dist_nr):
                
        ltrue = setPartitions.encode_partition(true.leaders[j], true.groups[j], true.G)
        lm = setPartitions.encode_partition(m.leaders[j], m.groups[j], m.G)

        # For all unordered pairs
        for i in range(m.G):
            for j in range(i+1,m.G):
                
                if ltrue[i] == ltrue[j]: # (i,j) is a positive
                    if lm[i] == lm[j]:
                        tp += 1
                    else:
                        fn += 1

                else: # (i,j) is a negative
                    if lm[i] == lm[j]:
                        fp += 1
                    else:
                        tn += 1             

    acc = float(tp+tn) / (tp + fn + tn + fp)

    return acc



def modelSelection(data,models, silent=False):
    """
    Computes model selection criterias NEC, BIC and AIC for a list of models. 

    @param data: DataSet object
    @param models: list of MixtureModel objects order with ascending number of components.
    
    @return: list of optimal components number according to [NEC, BIC, AIC], in that order. 
    """
    assert models[0].G == 1, "One component model needed for NEC."
    
    m_1 = models[0]
    data.internalInit(m_1)
    L_1 = get_loglikelihood(m_1,data)
    #P_1 = L_1 + m_1.prior.pdf(m_1)

    G_list  = [1]
    NEC = [1.0]
    BIC = [-2*L_1 - (m_1.freeParams * numpy.log(data.N))]
    AIC = [-2*L_1 + ( 2 * m_1.freeParams )]
    #bBIC = [ -2*P_1 - (m_1.freeParams * numpy.log(data.N)) ]  # test: BIC with MAP instead of ML
    for i in range(1,len(models)):
        m_i = models[i]
        G_list.append(m_i.G)
        (log_l,L_G) = m_i.EStep(data) 
        l = numpy.exp(log_l)        

        if m_i.G == 1:
            NEC.append(1.0)  # if G=1, then NEC = 1.0 by definition
        else:

            # entropy term of the NEC
            E_g = 0
            for g in range(m_i.G):
                for n in range(data.N):
                    if log_l[g,n] != float('-inf'):
                        E_g += l[g,n]  * log_l[g,n]
            E_g = -E_g
            NEC_G = E_g / ( L_G - L_1 )
            NEC.append(NEC_G)
 
        BIC_G = -2*L_G + (m_i.freeParams * numpy.log(data.N))
        BIC.append(BIC_G)
        
        AIC_G =  -2*L_G + ( 2 * m_i.freeParams )
        AIC.append(AIC_G)
        
    if not silent:
        print "NEC = ", NEC
        print "BIC = ", BIC
        print "AIC = ", AIC


    NEC_min = numpy.argmin(numpy.array(NEC,dtype='Float64'))
    AIC_min = numpy.argmin(numpy.array(AIC,dtype='Float64'))
    BIC_min = numpy.argmin(numpy.array(BIC,dtype='Float64'))

    if not silent:
        print G_list
        print '**',NEC_min
        print "Minimal NEC at G = "+str(G_list[NEC_min])+" with "+ str(NEC[NEC_min])
        print "Minimal BIC at G = "+str(G_list[BIC_min])+" with "+ str(BIC[BIC_min])
        print "Minimal AIC at G = "+str(G_list[AIC_min])+" with "+ str(AIC[AIC_min])

    return (NEC,BIC,AIC)

def kl_dist(d1,d2):
    """
    Kullback-Leibler divergence for two distributions. Only accept MultinomialDistribution and
    NormalDistribution objects for now.
    
    @param d1: MultinomialDistribution or NormalDistribution instance
    @param d2: MultinomialDistribution or NormalDistribution instance

    @return: Kullback-Leibler divergence between input distributions
    """
    # Kullback = sum[1..P](ln(SIGMA2/SIGMA1)) 
    # + sum[1..P](SIGMA1^2 / (2*(SIGMA2^2))) 
    # + sum[1..P]((MU1-MU2)^2 / (2*(SIGMA2^2))) - P/2
    if isinstance(d1,NormalDistribution) and isinstance(d2,NormalDistribution):
        res = ( (0.5 * numpy.log(d2.sigma**2/d1.sigma**2)) - 0.5 + d1.sigma**2/(2*d2.sigma**2)
              + (abs(d2.mu-d1.mu)**2)/(2*d2.sigma**2) )
        return res
    elif isinstance(d1,MultinomialDistribution) and isinstance(d2,MultinomialDistribution):
        assert d1.M == d2.M
        en = 0
        for i in range(d1.M):
            en += d1.phi[i] * numpy.log(d1.phi[i]/d2.phi[i])
        return en    
    elif isinstance(d1, ProductDistribution) and isinstance(d2, ProductDistribution):
        assert d1.dist_nr == d2.dist_nr == 1
        return kl_dist(d1[0], d2[0])
    
    else:
        raise TypeError, "Type mismatch or distribution not yet supported by kl_dist: "+str(d1.__class__)+", "+str(d2.__class__)


def sym_kl_dist(d1,d2):
    """
    Symmetric Kullback-Leibler divergence for two distributions. Only accept MultinomialDistribution and
    NormalDistribution objects for now.
    
    @param d1: MultinomialDistribution or NormalDistribution instance
    @param d2: MultinomialDistribution or NormalDistribution instance

    @return: Symmetric Kullback-Leibler divergence between input distributions
    """
    if ( (isinstance(d1,NormalDistribution) and isinstance(d2,NormalDistribution)) 
        or (isinstance(d1,MultinomialDistribution) and isinstance(d2,MultinomialDistribution)) ):
        d12 = kl_dist(d1,d2)
        d21 = kl_dist(d2,d1)
        dist = (d12 + d21 )/2.0
    elif isinstance(d1,MixtureModel) and isinstance(d2,MixtureModel):
        assert d1.G == d2.G, "Unequal number of components"
        d12 = 0
        d21 = 0
        for i in range(d1.G):
            d12 += d1.pi[i] * kl_dist(d1.components[i],d2.components[i])
            d21 += d2.pi[i] * kl_dist(d2.components[i],d1.components[i])
        dist = (d12 + d21 ) / 2.0
    else:
        raise TypeError,str(d1.__class__)+" != "+str(d2.__class__)

    if dist < 0.0:
        #raise ValueError,"Negative distance in sym_kl_dist."
        #print 'WARNING: Negative distance in sym_kl_dist.'
        return 0.0
    else: 
        return dist
        
def computeErrors(classes, clusters):
    """
    For an array of class labels and an array of cluster labels
    compute true positives, false negatives, true negatives and
    false positives.

    Assumes identical order of objects.

    Class and cluster labels can be arbitrary data types supporting
    '==' operator.
   
    @param classes: list of class labels (true labels)
    @param clusters: list of cluster labels (predicted labels)
   
    @return: Ratios for true positives, false negatives, true 
    negatives, false postitives   (tp, fn, tn, fp) 
    """
    
    #print 'true:',classes
   # print 'pred:', clusters
    
    assert len(classes) == len(clusters)
    tp = fn = tn = fp = 0
    
    classList = []
    clustList = []
    # samples with cluster or class label -1 are excluded 
    for i in xrange(len(classes)):
        if clusters[i] != -1 and classes[i] != -1:
            classList.append(classes[i])
            clustList.append(clusters[i])

        #else:
        #    print i,'discarded:' ,  classes[i] , clusters[i]       
                    
    # For all unordered pairs
    for i in xrange(len(classList)):
        for j in xrange(i+1, len(classList)):
            if classList[i] == classList[j]: # (i,j) is a positive
                if clustList[i] == clustList[j]:
                    tp += 1
                else:
                    fn += 1
                    
            else: # (i,j) is a negative
                if clustList[i] == clustList[j]:
                    fp += 1
                else:
                    tn += 1 
                    
    return (tp, fn, tn, fp)


def accuracy(classes, clusters):
    """
    Computes accuracy of a clustering solution
    
    @param classes: list of true class labels
    @param clusters: list of cluster labels   
    @return: accuracy
    """
    (tp, fn, tn, fp) = computeErrors(classes, clusters)     
    if (tp + tn) != 0:
        return float(tp + tn) / (tp + fp + tn + fn)
    else:   
        return 0.0

def sensitivity(classes, clusters):
    """
    Computes sensitivity of a clustering solution
    
    @param classes: list of true class labels
    @param clusters: list of cluster labels   
    @return: sensitivity
    """
    (tp, fn, tn, fp) = computeErrors(classes, clusters)     
    if(tp + fn) != 0:
        return  float(tp)/(tp + fn)
    else:
        return 0.0

def specificity(classes, clusters):
    """
    Computes specificity of a clustering solution
    
    @param classes: list of true class labels
    @param clusters: list of cluster labels   
    @return: specificity
    """

    (tp, fn, tn, fp) = computeErrors(classes, clusters)     
    if (tp+fp) != 0.0:
        return float(tp)/(tp + fp)
    else:
        return 0.0


def random_vector(nr,normal=1.0):
    """
    Returns a random probability vector of length 'nr'.
    Can be used to generate random parametrizations of a multinomial distribution with
    M = 'nr'.
    
    @param nr: lenght of output vector
    @param normal: sum over output vector, default 1.0
    
    @return: list with random entries sampled from a uniform distribution on [0,1] and normalized to 'normal'
    """

    alpha = numpy.array([1.0]*nr)

    p = _C_mixextend.wrap_gsl_dirichlet_sample(alpha,nr)

    if float(normal) != 1.0:
        p = p*normal
    return p.tolist()


def variance(data):
    mean = data.mean()
    s = 0.0
    for i in range(len(data)):
        s += (data[i] - mean)**2
    return s/(len(data)-1)    

def entropy(p):
    """
    Returns the Shannon entropy for the probilistic vector 'p'.
    
    @param p: 'numpy' vector that sums up to 1.0
    """
    res = 0.0
    for i in range(len(p)):
        if p[i] != 0.0:
            res += p[i] * math.log(p[i],2)
    return -res

def get_loglikelihood(mix_model, data):  # old implementation XXX
    # matrix of posterior probs: components# * (sequence positions)#
    l = numpy.zeros((mix_model.G,len(data)),dtype='Float64')
    col_sum = numpy.zeros(len(data),dtype='Float64')
    for i in range(mix_model.G):            
        l[i] = numpy.log(mix_model.pi[i]) + mix_model.components[i].pdf(data)
    for j in range(len(data)):
        col_sum[j] = sumlogs(l[:,j]) # sum over jth column of l    
    log_p = numpy.sum(col_sum)    
    return log_p

def get_posterior(mix_model, data,logreturn=True): 
    # matrix of posterior probs: components# * (sequence positions)#
    log_l = numpy.zeros((mix_model.G,len(data)),dtype='Float64')
            
    # computing log posterior distribution 
    for i in range(mix_model.G):            
        log_l[i] = math.log(mix_model.pi[i]) + mix_model.components[i].pdf(data) 


    # computing data log likelihood as criteria of convergence
    # log_l is normalized in-place and likelihood is returned as log_p
    log_p = _C_mixextend.get_normalized_posterior_matrix(log_l)

    if logreturn == True:                
        return log_l
    else:
        return numpy.exp(log_l)
        
        
        
## function sumlogs is borrowed from GQLMixture.py
def sumlogs_purepy(a):
    """ Given a Numeric.array a of log p_i, return log(sum p_i)

        Uses (assuming p_1 is maximal):
        log(\Sum p_i) = log(p_1) + log( 1 + \Sum_{i=2} exp(log(p_i) - log(p_1)))

        NOTE: The sumlogs functions returns the sum for values != -Inf

    """
    m = max(a) # Maximal value must not be unique
    result = 0.0
    #minus_infinity = -float(1E300)
    for x in a:
        if x >= m: # For every maximal value
            result += 1.0
        else:
            if x == float('-inf'): # zero probability, hence
                # -Inf log prob. Doesnt contribute
                continue
            x = x - m
            # Special case to avoid numerical problems
            if x < -1.0e-16: # <=> |x| >  1.0e-16
                result += numpy.exp(x) 
            else: # |x| <  1.0e-16 => exp(x) = 1
                result += 1.0

    result = numpy.log(result)
    result += m
    return result


def sumlogs(a):
    """
    Call to C extension function sum_logs.
    """
    m = max(a) # Maximal value must not be unique
    result = _C_mixextend.sum_logs(a,m)
    return result 


def matrixSumlogs(mat):
    """ 
    Call to C extension function matrix_sum_logs 
    """
    return _C_mixextend.matrix_sum_logs(mat)


def dict_intersection(d1, d2):
    """
    Computes the intersections between the key sets of two Python dictionaries. 
    Returns another dictionary with the intersection as keys.
    
    @param d1: dictionary object
    @param d2: dictionary object

    @return: dictionary with keys equal to the intersection of keys between d1 and d2.
    
    """
    int_dict = {}
    for e in d2:
        if d1.has_key(e):
            int_dict[e] = 1

    return int_dict





#----------------------------------- File IO --------------------------------------------------
# XXX  While functional the flat file based file IO is somewhat crude. XXX
# XXX  The whole thing ought to be redone in XML at some point.        XXX

def writeMixture( model, fileName,silent=False):
    """
    Stores model parameters in file 'fileName'. 
    
    @param model: MixtureModel object
    @param fileName: file name the model is to be written to
    """
    f = open(fileName, 'w')
    if isinstance(model,labeledBayesMixtureModel):
        head = 'labelBayesMix'
    elif isinstance(model,BayesMixtureModel)    :
        head = "BayesMix"
    elif isinstance(model,MixtureModel):
        head = "Mix"
    else:
        raise TypeError    
    
    if not model.struct: 
        l = str(';'+head) +";"+str(model.G) + ";" + str(model.pi.tolist())+ ";"+ str(model.compFix) +"\n"
    else:
        l = str(';'+head) +";"+str(model.G) + ";" + str(model.pi.tolist())+ ";"+ str(model.compFix) + ";" + str(model.leaders)+";"+str(model.groups) +"\n"
    
    f.write(l)
    for i in range(model.G):
        l = model.components[i].flatStr(0)
        f.write(l)
    
    if head == "BayesMix" or head == 'labelBayesMix':
        l = model.prior.flatStr(0)
        f.write(l)
    
    if not silent:
        print "Model written to file " + str(fileName)+"."
    f.close()
        
def readMixture(fileName):
    """
    Reads model from file 'fileName'.
    
    @param fileName: file to be read
    
    @return: MixtureModel object
    """
    f = open(fileName, 'r')
    s = chomp(f.readline())
    struct = 0 
    if len(split(s,';')) == 5:  # MixtureModel object
        [offset,head,G,pi,compFix] = split(s,';')
        leaders = None
        groups = None
    elif len(split(s,';')) == 7:  # BayesMixtureModel object
        struct = 1
        [offset,head,G,pi,compFix,leaders,groups] = split(s,';')
    else:
        raise IOError, 'Flat file format not recognized.'
    if leaders and groups:
        mixModel = parseMix(f,head,int(G),simple_eval(pi),simple_eval(compFix),simple_eval(leaders),simple_eval(groups))
    else:
        mixModel = parseMix(f,head,int(G),simple_eval(pi),simple_eval(compFix))            
    f.close()
    return mixModel

def parseMix(fileHandle,mtype,G,pi,compFix,leaders = None,groups = None):
    """
    Parses a flat file for a mixture model. Internal function, is invoked from
    readMixture.
    
    """
    components = []
    while len(components) < G:
        components.append(parseFile(fileHandle))            

    if mtype == 'Mix':
        m = MixtureModel(G,pi,components,compFix=compFix)

    elif mtype == 'labelBayesMix':
        prior = parseFile(fileHandle)
        if sum(compFix) > 0: # XXX pass compFix if it is not trivial
            m = labeledBayesMixtureModel(G,pi,components,prior,compFix=compFix)
        else:
            m = labeledBayesMixtureModel(G,pi,components,prior)

    elif mtype == 'BayesMix':
        prior = parseFile(fileHandle)
        if sum(compFix) > 0: # XXX pass compFix if it is not trivial
            m = BayesMixtureModel(G,pi,components,prior,compFix=compFix)
        else:
            m = BayesMixtureModel(G,pi,components,prior)
    
    else:
        raise TypeError
    if leaders and groups:
        m.initStructure()            
        m.leaders = leaders
        m.groups = groups
        for i in range(m.dist_nr):
            for lead in m.leaders[i]:
                for g in m.groups[i][lead]:
                    if not m.components[lead][i] == m.components[g][i]:
                        raise IOError, 'Incompatible CSI structure and parameter values in parseMix.'
                    m.components[g][i] = m.components[lead][i]
    return m 

def parseProd(fileHandle,true_p):
    """
    Internal function. Parses product distribution.
    """
    distList = []
    p = 0
    while p < true_p:
        d = parseFile(fileHandle)
        distList.append(d)
        p += d.p
    return ProductDistribution(distList)

def parseMixPrior(fileHandle,nr_dist,structPrior,nrCompPrior):
    c = []
    piPrior = parseFile(fileHandle)
    for i in range(nr_dist):
        p = parseFile(fileHandle)
        c.append(p)
    return MixtureModelPrior(structPrior, nrCompPrior,piPrior, c)

def parseDirichletMixPrior(fileHandle,G,M,pi):
    dC = [ parseFile(fileHandle) for i in range(G)  ]
    return DirichletMixturePrior(G,M,pi,dC)    
    

def parseFile(fileHandle):
    """
    Internal function. Parses flat files.
    """
    s = chomp(fileHandle.readline())
    l = split(s,';')

    if l[1] == "Mix":
        [offset,head, G, pi,compFix] = l
        return parseMix(fileHandle,head,int(G),simple_eval(pi),simple_eval(compFix))
    elif l[1] ==  "Norm":
        [offset,head, mu, sigma] = l
        return NormalDistribution(float(mu),float(sigma))
    elif l[1] ==  "Exp":
        [offset,head, lambd] = l
        return ExponentialDistribution( float(lambd) )
    elif l[1] ==  "Mult":
        [offset,head, N,M,phi,alphabet,parFix] = l
        alph = Alphabet(simple_eval(alphabet))
        return MultinomialDistribution(int(N), int(M), simple_eval(phi),alph,simple_eval(parFix) )
    elif l[1] ==  "Discrete":
        [offset,head,M,phi,alphabet,parFix] = l
        alph = Alphabet(simple_eval(alphabet))
        return DiscreteDistribution( int(M), simple_eval(phi),alph,simple_eval(parFix) )
    elif l[1] ==  "MultiNormal":
        [offset,head,p,mu,sigma] = l
        # XXX the tokenize package used in simple_eval cannot deal with negative values in
        # mu or sigma. A hack solution to that would be to change simple_eval to a direct
        # call to eval in the line below. This carries all the usual implications for security.
        return MultiNormalDistribution( int(p), simple_eval(mu), simple_eval(sigma) )
    elif l[1] ==  "Dirichlet":
        [offset,head,M,alpha] = l
        return DirichletDistribution(int(M),simple_eval(alpha))
    elif l[1] ==  "DirichletPr":
        [offset,head,M,alpha] = l
        return DirichletPrior(int(M),simple_eval(alpha))
    elif l[1] ==  "NormalGamma":
        [offset,head,mu, kappa, dof, scale] = l
        return NormalGammaPrior( float(mu),float(kappa), float(dof),float(scale) )
    elif l[1] ==  "PriorForDirichlet":
        [offset,head,M, eta] = l
        return PriorForDirichletDistribution(int(M),simple_eval(eta))
    elif l[1] ==  "Prod":
        [offset,head, p] = l
        return parseProd(fileHandle,int(p))
    elif l[1] == "MixPrior":
	    #;MixPrior;4;0.7;0.7        
        [offset,head,nr_dist,structPrior,nrCompPrior] =  l  
        return parseMixPrior(fileHandle,int(nr_dist),float(structPrior),float(nrCompPrior) )
    elif l[1] == "DirichMixPrior":
    	#;DirichMixPrior;3;5;[ 0.3  0.3  0.4]
        [offset,head,G,M,pi] =  l          
        return parseDirichletMixPrior(fileHandle,int(G),int(M),simple_eval(pi))
    else:
        raise TypeError, "Unknown keyword: " + str(l[1])

def chomp(string):
    """
    Removes a newline character from the end of the string if present

    @param string: input string

    @return: the argument without tailing newline.
    
    """
    if string[-1] == "\n" or string[-1] == "\r":
        return string[0:-1]
#    elif string[len(string)-4:] == "<cr>":    
#        return string[0:-4]        
    else:
        return string

# The following functions come courtesy of the people at comp.lang.python
def sequence(next, token, end):
    out = []
    token = next()
    while token[1] != end:
        out.append(atom(next, token))
        token = next()
        if token[1] == "," or token[1] == ":":
            token = next()
    return out

def atom(next, token):
    if token[1] == "(":
        return tuple(sequence(next, token, ")"))
    elif token[1] == "[":
        return sequence(next, token, "]")
    elif token[1] == "{":
        seq = sequence(next, token, "}")
        res = {}
        for i in range(0, len(seq), 2):
            res[seq[i]] = seq[i+1]
        return res
    elif token[0] in (tokenize.STRING, tokenize.NUMBER):
        return eval(token[1]) # safe use of eval!
    raise SyntaxError("malformed expression (%s)" % token[1])

def simple_eval(source):
    src = cStringIO.StringIO(source).readline
    src = tokenize.generate_tokens(src)
    src = (token for token in src if token[0] is not tokenize.NL)
    res = atom(src.next, src.next())
    if src.next()[0] is not tokenize.ENDMARKER:
        raise SyntaxError("bogus data after expression")
    return res 
