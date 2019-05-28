#!/usr/bin/env python

import sys
import os
import math
import re
import cPickle
from xml.dom import minidom

import itertools
import tempfile

class PredictorError(Exception):
    """Exception raised for errors in the ARBMatrix and SMMMatrix classes."""
    def __init__(self, value):
            self.value = value
    def __str__(self):
            return self.value

class SturnioloMatrix:
    """Can load and save Sturniolo matrices in multiple formats and use them to score sequences """
    def __init__(self):
        self.offset= None
        self.length= None
        self.mat={}

    def loadTextfile(self, infile):
        lines=infile.readlines()
        self.mat.clear()
        self.length=int(lines[0].split()[1])
        for line in lines[1:21]:
            entries = line.split()
            numbers = []
            for e in entries[1:]:
                numbers.append(float(e))
            if len(numbers)!=self.length:
                raise PredictorError("Invalid number of columns in Sturniolo matrix: " + str(len(numbers)), " expected: " + str(self.length))
            self.mat[line[0]]=tuple(numbers)
        self.offset = float(0)

    def predict(self,sequence):
        scores=[]
        for nterm in range(len(sequence)-(self.length-1)):
            score=self.offset
            for pos in range(self.length):
                amino_acid=sequence[nterm+pos]
                try:
	            score+=self.mat[amino_acid][pos]
                except:
                    raise PredictorError("""Invalid character '%c' in sequence '%s'""" % (amino_acid, sequence))
	    if (score < -12.5):
		score = -12.5
            scores.append(score)
        return(tuple(scores))

    def consensus_predict(self,sequence):
        scores=[]
        for nterm in range(len(sequence)-14):
            score=-13
	    core=""
	    for i in range(7):
		temp_score = 0
	        for pos in range(self.length):
                    amino_acid=sequence[nterm+pos+i]
                    try:
	                temp_score+=self.mat[amino_acid][pos]
                    except:
                        raise PredictorError("""Invalid character '%c' in sequence '%s'""" % (amino_acid, sequence))
      	        if (temp_score < -12.5):
		    temp_score = -12.5
		if score < temp_score:
		    score = temp_score
		    core = sequence[nterm+i:nterm+i+9]
 	    consensus_score=(core, score)
            scores.append(consensus_score)
        return(tuple(scores))

    def log_predict(self,sequence):
        scores=[]
        for nterm in range(len(sequence)-(self.length-1)):
            score=self.offset
            for pos in range(self.length):
                amino_acid=sequence[nterm+pos]
                try:
                    score+=self.mat[amino_acid][pos]
                except:
                    raise PredictorError("""Invalid character '%c' in sequence '%s'""" % (amino_acid, sequence))
            scores.append(score)
        return(tuple(scores))

    def loadXML(self, xml):
        fin = minidom.parseString(xml)

        alphabet = fin.getElementsByTagName("Alphabet")[0].firstChild.data
        self.length = int(fin.getElementsByTagName("SequenceLength")[0].firstChild.data)
        smat = fin.getElementsByTagName("SeqMatrix")[0]
        self.offset = float(smat.getElementsByTagName("Offset")[0].firstChild.data)
        mat={}
        mcoefs=smat.getElementsByTagName("MatCoef")
        for mcoef in mcoefs:
            pos=mcoef.getElementsByTagName("Position")[0].firstChild.data
            let=mcoef.getElementsByTagName("Letter")[0].firstChild.data
            val=mcoef.getElementsByTagName("Value")[0].firstChild.data
            mat[(let,int(pos))]=float(val)
        self.mat.clear()

        for let in alphabet:
            vals=[]
            for pos in range(1,self.length+1):
                vals.append(mat[(let,pos)])
            self.mat[let]=tuple(vals)

    def saveTextFile(self, outfile):
        outfile.write("NumCols:\t" + str(self.length) +"\n")
        for letter in sorted(self.keys()):
            outfile.write(letter)
            for val in self.mat[letter]:
                outfile.write("\t" + str(val))
            outfile.write("\n")
        outfile.write(str(self.offset))

    def pickle_dump(self, file_name):
        fout = open(file_name,"wb")
        cPickle.dump(self.length, fout)
        cPickle.dump(self.mat,fout)
        cPickle.dump(self.offset,fout)
        fout.close()

    def pickle_load(self, file_name):
        fin = open(file_name,"rb")
        self.length = cPickle.load(fin)
        self.mat = cPickle.load(fin)
        self.offset = cPickle.load(fin)
        fin.close()

class COMBLIBMatrix:
    """Can load and save ARB matrices in multiple formats and use them to score sequences """
    def __init__(self):
        self.slope= None
        self.intercept= None
        self.length= None
        self.mat={}
	self.geomean = 1.0

    def predict(self,sequence):
        scores=[]
        for nterm in range(len(sequence)-(self.length-1)):
            score=0
            for pos in range(self.length):
                amino_acid=sequence[nterm+pos]
                try:
                    score += self.mat[amino_acid][pos]
                except:
                    raise PredictorError("""Invalid character '%c' in sequence '%s'""" % (amino_acid, sequence))
	    score = 0 - score
	    EIC50 = math.pow(10, score)*self.geomean
	    logEIC50 = math.log10(self.EIC50)	    
	    logPIC = (logEIC50 - self.intercept)/self.slope

	    PIC = math.pow(10, logPIC)
	    if PIC >= 1000000:
	    	PIC = 1000000
	    if PIC <= 0.01:
	    	PIC = 0.01
            scores.append(PIC)
        return(tuple(scores))

    def consensus_predict(self,sequence):
        scores=[]
        for nterm in range(len(sequence)-14):
            score=100000000
	    core=""
	    for i in range(7):
		temp_score = 0
		temp_tscore = 0
                for pos in range(self.length):
                	amino_acid=sequence[nterm+pos+i]
	                try:
        	            temp_score += self.mat[amino_acid][pos]
                	except:
	                    raise PredictorError("""Invalid character '%c' in sequence '%s'""" % (amino_acid, sequence))
	    	temp_tscore = 0 - temp_score
	        EIC50 = math.pow(10, temp_tscore)*float(self.geomean)
	        logEIC50 = math.log10(EIC50)	
         	logPIC = (logEIC50 - self.intercept)/self.slope
	        PIC = math.pow(10, logPIC)
	        if PIC >= 1000000:
	    		PIC = 1000000
	        if PIC <= 0.01:
	    		PIC = 0.01
		if PIC <= score:
		    score = PIC
		    core = sequence[nterm+i:9+i+nterm]
	    consensus_score = (core, score)
            scores.append(consensus_score)
        return(tuple(scores))

    def pickle_dump(self, file_name):
        fout=open(file_name,"wb")
        cPickle.dump(self.length, fout)
        cPickle.dump(self.mat,fout)
        cPickle.dump(self.slope,fout)
        cPickle.dump(self.intercept,fout)
	cPickle.dump(self.geomean, fout)
        fout.close()

    def pickle_load(self, file_name):
        fin = open(file_name,"rb")
        self.length = cPickle.load(fin)
        self.mat = cPickle.load(fin)
        self.slope = cPickle.load(fin)
        self.intercept = cPickle.load(fin)
        self.geomean = cPickle.load(fin)
        fin.close()

class SMMPredictor:
    def __init__(self, smm_dir):
        self.length = 9
        self.smm_dir = smm_dir
        self.path_executable= smm_dir + "/netMHCII"

    def peptide_predictions(self, sequence, allele):
        ic50s=[]
        infile = tempfile.NamedTemporaryFile(prefix='input_', dir='/tmp',)
        outfile = tempfile.NamedTemporaryFile(prefix='output_', dir='/tmp',)
        infile.write(">test\n"+sequence + "\n")
        
        if not re.compile("H2").match(allele): 
            allele = "HLA-%s" %allele.replace('*','').replace(':','')
        else:
            allele = allele.replace('H2','H-2')
        
        data = os.path.join(self.smm_dir, 'data', allele)
        
        infile.seek(0)
        if os.path.exists(data):
            os.spawnl (os.P_WAIT, self.path_executable, "netMHCII", "-a", allele, infile.name, outfile.name)
        else:
            sys.exit(0)
        infile.close()
        
        outfile = open(outfile.name, "r")
        lines = outfile.readlines()
        outfile.close()

        all_lists=[]
        newlist = [list(g) for k,g in itertools.groupby(lines,lambda x: re.match('-',x)) if not k]
        
        i = 2
        while i < len(newlist):
            all_lists.append(newlist[i])
            i = i + 3
    
        check_list=[]
        for i in xrange(len(all_lists)):
            ic50s=[]
            for l in all_lists[i]:
                mynewlist = l.split()
                allele = mynewlist[0]
                score = float(mynewlist[5])
                core = mynewlist[3]
                temp_list = (core, score)
                ic50s.append(temp_list)
            check_list.append(ic50s)
 
        outfile.close()
        return(check_list)

    def predict(self, sequence):
        peptides=[]
        for nterm in range(len(sequence)-(self.length-1)):
            peptides.append(sequence[nterm:nterm+self.length])
        return(self.peptide_predictions(peptides))

    def multi_predict(self, sequences):
        peptides=[]
        for sequence in sequences:
            for nterm in range(len(sequence)-(self.length-1)):
                peptides.append(sequence[nterm:nterm+self.length])
        scores=self.peptide_predictions(peptides)
        result=[]
        start=0
        for sequence in sequences:
            if len(sequence)>=self.length:
                end=start+len(sequence)-(self.length-1)
                result.append(scores[start:end])
                start=end
            else:
                result.append([])

        return(result)

class NNPredictor:
    def __init__(self, nn_dir):
        self.length = 9
        self.nn_dir = nn_dir
        self.path_executable= nn_dir + "/netMHCII-2.2"

    def peptide_predictions(self, sequence, allele):
        ic50s=[]
        infile = tempfile.NamedTemporaryFile(prefix='input_', dir='/tmp',)
        outfile = tempfile.NamedTemporaryFile(prefix='output_', dir='/tmp',)
        infile.write(">test\n"+sequence + "\n")
        
        if not re.compile("H2").match(allele): 
            allele = "HLA-%s" %allele.replace('*','').replace(':','')
        else:
            allele = allele.replace('H2','H-2')
        
        data = os.path.join(self.nn_dir, 'data', allele)
        
        infile.seek(0)
        if os.path.exists(data):
            os.spawnl (os.P_WAIT, self.path_executable, "netMHCII-2.2", "-a", allele, infile.name, outfile.name)
        else:
            sys.exit(0)
        infile.close()
        
        outfile = open(outfile.name, "r")
        lines = outfile.readlines()
        outfile.close()
        
        all_lists=[]
        newlist = [list(g) for k,g in itertools.groupby(lines,lambda x: re.match('-',x)) if not k]

        i = 2
        while i < len(newlist):
            all_lists.append(newlist[i])
            i = i + 3

        check_list=[]
        for i in xrange(len(all_lists)):
            ic50s=[]
            for l in all_lists[i]:
                mynewlist = l.split()
                allele = mynewlist[0]
                score = float(mynewlist[5])
                core = mynewlist[3]
                temp_list = (core, score)
                ic50s.append(temp_list)
            check_list.append(ic50s)

        outfile.close()
        return(check_list)
    def predict(self, sequence):
        peptides=[]
        for nterm in range(len(sequence)-(self.length-1)):
            peptides.append(sequence[nterm:nterm+self.length])
        return(self.peptide_predictions(peptides))

    def multi_predict(self, sequences):
        peptides=[]
        for sequence in sequences:
            for nterm in range(len(sequence)-(self.length-1)):
                peptides.append(sequence[nterm:nterm+self.length])
        scores=self.peptide_predictions(peptides)
        result=[]
        start=0
        for sequence in sequences:
            if len(sequence)>=self.length:
                end=start+len(sequence)-(self.length-1)
                result.append(scores[start:end])
                start=end
            else:
                result.append([])

        return(result)

class Net2Predictor:
    def __init__(self, net_dir):
        self.length = 9
        self.path_executable= net_dir + "/netMHCIIpan"

    def peptide_predictions(self, sequence, allele):
        ic50s=[]
        
        infile = tempfile.NamedTemporaryFile(prefix='input_', dir='/tmp',)
        outfile = tempfile.NamedTemporaryFile(prefix='output_', dir='/tmp',)
        infile.write(">test\n"+sequence + "\n")

        if not re.compile("H2").match(allele):
            allele = allele.replace(':','')
            if "DRB" in allele: allele = allele.replace('*','_')
            else: allele = "HLA-%s" %allele.replace('*','')
        else:
            allele = allele.replace('H2','H-2')
            
        infile.seek(0)
        
        import subprocess as sub
        proc = sub.Popen([self.path_executable, '-a', allele, '-f', infile.name], stdout=sub.PIPE, stderr=sub.PIPE)
        output, errors = proc.communicate()
        if errors:
            msg = "NetMHCIIpan Predictor did not execute. Path: %s, input :%s, allele: %s, error: %d." % (self.path_executable, infile.name, allele, errors)
            raise PredictorError(msg)
        
        lines = output.splitlines(True)
        
        all_lists=[]
        newlist = [list(g) for k,g in itertools.groupby(lines,lambda x: re.match('-',x)) if not k]

        i = 2
        while i < len(newlist):
            all_lists.append(newlist[i])
            i = i + 4    
        
        check_list=[]
        for i in xrange(len(all_lists)):
            ic50s=[]
            for l in all_lists[i]:
                mynewlist = l.split()
                allele = mynewlist[1]
                score = float(mynewlist[7])
                core = mynewlist[5]
                temp_list = (core, score)
                ic50s.append(temp_list)
            check_list.append(ic50s)

        outfile.close()
        return(check_list)

	def predict(self, sequence):
		peptides=[]
		for nterm in range(len(sequence)-(self.length-1)):
			peptides.append(sequence[nterm:nterm+self.length])
		return(self.peptide_predictions(peptides))

	def multi_predict(self, sequences):
		peptides=[]
		for sequence in sequences:
			for nterm in range(len(sequence)-(self.length-1)):
				peptides.append(sequence[nterm:nterm+self.length])
		scores=self.peptide_predictions(peptides)
		result=[]
		start=0
		for sequence in sequences:
			if len(sequence)>=self.length:
				end=start+len(sequence)-(self.length-1)
				result.append(scores[start:end])
				start=end
			else:
				result.append([])

		return(result)
     
    
