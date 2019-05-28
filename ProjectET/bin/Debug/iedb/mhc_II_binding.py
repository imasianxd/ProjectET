#!/usr/bin/env python

from bisect import bisect, bisect_left
import os
import traceback, time
import math
from seqpredictorII import SMMPredictor, SturnioloMatrix, COMBLIBMatrix, NNPredictor, Net2Predictor
import ConfigParser
import sys
import re
import csv

class InputError(Exception):
    """Exception raised for errors in the input."""

    def __init__(self, value):
            self.value = value

    def __str__(self):
            return self.value

class UnexpectedInputError(Exception):
    """Exception raised for errors in the input."""

    def __init__(self, value):
            self.value = value

    def __str__(self):
            return self.value

class Proteins:
    """Contains a list of protein sequences and names and several conversion functions."""
    def __init__(self, form=None):
        self.sequences=[]
        self.names=[]
        if not form==None:
            self.extractForm(form)

    def add_protein(self,sequence, name=""):
        sequence=sequence.strip().upper()
        for amino_acid in sequence:
            if not amino_acid in "ACDEFGHIKLMNPQRSTVWY":
                raise InputError("Sequence: '%s'<br /> Contains an invalid character: '%c'<br/> At position %d" %(sequence, amino_acid, sequence.find(amino_acid)))
        self.sequences.append(sequence)

        if name=="":
            name = "sequence %d" %(len(self.sequences))
        self.names.append(name)


    def extractFasta(self, fasta):
        input_sequences = fasta.split(">")
        if len(input_sequences)<2:
            raise InputError("Invalid Fasta format: </br>No '>' found.")
        for i in input_sequences[1:]:
            if(len(i) > 0):
                end_of_name = i.find("\n")
                if end_of_name == -1:
                        raise InputError("Invalid Fasta format: </br>No Protein sequence found between two >names ")
                name = i[:end_of_name]
                seq = i[end_of_name:].split()
                self.add_protein("".join(seq),name)

    def extractForm(self, form):
        sequence_format = form["sequence_format"]

        if form['sequence_file'] != "":
            FILE = open(form["sequence_file"], "r")
            input_sequences = FILE.read()
        else:
            input_sequences = form["sequence_text"]

        if sequence_format == "auto":
            if ">" in input_sequences:
                sequence_format="fasta"
            else:
                seq = input_sequences.split()
                sequence_format="space_separated"
                for s in seq:
                    if len(s)>25:
                        sequence_format="one_sequence"
                        break

        if sequence_format == "fasta":
            self.extractFasta(input_sequences)
        elif sequence_format == "one_sequence":
            seq = input_sequences.split()
            self.add_protein("".join(seq))
        elif sequence_format == "space_separated":
            seqs = input_sequences.split()
            for seq in seqs:
                self.add_protein(seq)

    def xmlEncodeMHCPredictInput(self):
        mhcPredictInput  = """<?xml version="1.0" encoding="UTF-8"?><MhcPredictInput>"""
        for (name, sequence) in zip(self.names, self.sequences):
            mhcPredictInput += '<Predict><Sequence name="' + name + '">' + sequence + '</Sequence></Predict>'
        mhcPredictInput += """</MhcPredictInput>"""
        return mhcPredictInput


class PredictionTable:
    """Generates a table from a set of sequences and predictions"""
    def __init__(self, form):
        self.processing_column={"Proteasome_score":4, "TAP_score":5, "MHC_score":6, "Processing_score":7, "Total_score":8}
        self.sort_output = form["sort_output"]
        self.row_data = []

    def extractCutoff(self, form):
        if form["cutoff_type"] == "none":
            self.use_cutoff=False
        else:
            self.use_cutoff=True
            try:
                self.cutoff_value = float(form["cutoff_value"])
            except:
                raise InputError("Invalid entry for output cutoff <br/>Please enter a numeric value.")
            if form["cutoff_type"] != "MHC_IC50":
                try:
                    self.cutoff_column = self.processing_column[form["cutoff_type"]]
                except:
                    raise UnexpectedInputError("Invalid cutoff option: " + form["cutoff_type"])

    def sortRows(self):
        if self.sort_output == "MHC_IC50":
            if form["pred_method"] == "Sturniolo":
                self.row_data.sort(lambda x,y: cmp(x[-1], y[-1]))
                self.row_data.reverse()
            else:
                self.row_data.sort(lambda x,y: cmp(x[-1], y[-1]))
        elif self.sort_output =="position_in_sequence":
            pass
        else:
            try:
                col = self.processing_column[self.sort_output]
            except:
                raise UnexpectedInputError("Invalid sorting option: " + self.sort_output)
            self.row_data.sort(lambda x,y: -cmp(x[col], y[col]))

    def add_rows_binding(self, allele, pep_length, proteins, scores):
        if form["pred_method"] == "consensus3" or form["pred_method"] == "IEDB_recommended":
            for (i,(sequence, predictions)) in enumerate(zip(proteins.sequences,scores)):
                dummy = []
                for (k, (dummy)) in enumerate(predictions):
                    peptide_sequence = sequence[k : k + 15]
                    peptide_position = "%s:%s-%s" %(i+1, k+1, k+15)
                    try:
                        temp_index = peptide_sequence.index(peptide_sequence)
                    except:
                        raise ValueError("Core sequence and peptide sequence doesn't match!!!")
                    dummy_list = []
                    dummy_list.append(allele)
                    dummy_list.append(peptide_position)
                    dummy_list.append(peptide_sequence)
                    for item in dummy:
                        dummy_list.append(item)
                    self.row_data.append(tuple(dummy_list))
        else:
            for (i,(sequence, predictions)) in enumerate(zip(proteins.sequences,scores)):
                for (k, (core, prediction)) in enumerate(predictions):
                    peptide_sequence = sequence[k : k + 15]
                    core_sequence = core
                    try:
                        temp_index = peptide_sequence.index(core_sequence)
                    except:
                        raise ValueError("Core sequence and peptide sequence doesn't match. you are here!!!")
                    core_position = "%s:%s-%s" %(i+1, k + 1 + temp_index, k + int(pep_length) + temp_index)
                    self.row_data.append((allele, core_position, core_sequence, peptide_sequence, prediction))

    def format_binding(self, proteins, results):
        for(pep_length, allele, scores) in results:
            self.add_rows_binding(allele,pep_length, proteins, scores)
        self.sortRows()
        return self.row_data

    def add_rows_processing(self, allele, pep_length, proteins, mhc_scores, tap_scores, prot_scores):
        prot_offset = pep_length-6
        for (i,(sequence, mhc, tap, prot)) in enumerate(zip(proteins.sequences,mhc_scores, tap_scores, prot_scores)):
            for nterm in range(len(sequence)-(pep_length-1)):
                mhc_ic50  = mhc[nterm]
                mhc_score = -math.log10(mhc_ic50)
                tap_score = tap[nterm]
                prot_score = prot[nterm+prot_offset]
                proc_score = prot_score + tap_score
                total_score = mhc_score+tap_score+ prot_score
                sequence_position = "%s:%s-%s" %(i+1, nterm + 1, nterm + int(pep_length))
                scanned_sequence = sequence[nterm : nterm + pep_length]
                row = (allele, sequence_position, pep_length, scanned_sequence, prot_score, tap_score, mhc_score, proc_score, total_score, mhc_ic50)
                self.row_data.append(row)

    def format_processing(self, proteins, mhc_results, tap_results, prot_scores):
        for (pep_length, allele, mhc_scores) in mhc_results:
            for (pep_length2, tap_scores) in tap_results:
                if pep_length==pep_length2:
                    self.add_rows_processing(allele, pep_length, proteins, mhc_scores, tap_scores, prot_scores)
                    break
        self.sortRows()
        return self.row_data

    def html_table_input_sequences(self, proteins):
        result="""
        <table id='input_table' class="result" border="2" >
            <thead><tr><th>#</th><th>Name</th><th>Sequence</th></thead>
            <tbody>
        """
        for i,(name, sequence) in enumerate(zip(proteins.names, proteins.sequences)):
            if i % 2==0:
                rowstyle = 'class="even_row"'
            else:
                rowstyle = ""
            blocklen = 50
            if len(sequence)>blocklen:
                seqblock=""
                for block in range(0,len(sequence)-blocklen,blocklen):
                    seqblock+="""%s<br/>""" %sequence[block:block+blocklen]
                seqblock+=sequence[block+blocklen:]
                result+='<tr %s><td>%d</td><td>%s</td><td class ="sequence">%s</td></tr>' % (rowstyle, i+1,name, seqblock)
            else:
                result+='<tr %s><td>%d</td><td>%s</td><td class = "sequence">%s</td></tr>' % (rowstyle, i+1,name, sequence)
        result+="</tbody></table>"
        return result


class MHCBindingPredictions:
    """Class to generate MHC binding predictions"""
    def __init__(self, form):
        self.method = form["pred_method"]
        form_allele = form["allele"]
        # remove 'HLA-' prefix if exists and change '/' to '-'
        form_allele = form_allele.replace("HLA-", "").replace("/","-")
        
        #Validate method selection
        method_list = {"consensus3" : 1, "comblib":2, "smm_align":3, "sturniolo":4, "nn_align":5, "NetMHCIIpan":6, "IEDB_recommended":7}
        if not method_list.has_key(self.method):
            raise UnexpectedInputError("Selected prediction method '%s' does not exist. To see all available methods, please type command:\n$ python mhc_II_binding.py method\n" % self.method)

        #Validate selections and / or determine values for "All lengths" / "All alleles" selection
        self.tool_selection = []
        self.tool_selection.append((form_allele, 9))

        if len(self.tool_selection)==0:
        	raise(UnexpectedInputError("No prediction method available for allele='%s', method='%s'. To see all available methods and corresponding allele list, please type command:\n$ python mhc_II_binding.py allele\n" % (form_allele, self.method)))

    def tool_location(self, allele, peplength):
        return("./tool_data/MHCII/%s/%s-%s" %(self.method, allele, peplength))
        
    def comblib_tool_location(self, allele, peplength):
        return("./tool_data/MHCII/%s/%s-%s" %(self.method, allele, peplength))
        
    def consensus_tool_location(self, allele, peplength, method):
        return("./tool_data/MHCII/%s/%s-%s" %(method, allele, peplength))

    def consensus_perc_location(self, allele,  method):
        return("./tool_data/MHCII/%s/perc_%s" %(method, allele))
        
    def getMedian(self, values):
        if len(values) % 2 == 1:
            return values[(len(values)+1)/2-1]
        else:
            lower = values[len(values)/2-1]
            upper = values[len(values)/2]
            return (float(lower + upper)) / 2

    def predict(self,proteins):
        results=[]
        if self.method=="consensus3":
            con_status = {}
            con_list = open("./tool_data/MHCII/consensus/consensus_allele_list.txt").read().strip().split("\n")
            for con_element in con_list:
                con_arr = con_element.split("\t")
                con_status[con_arr[0]] = con_arr[1:]
            FILE = open("/home/dorjee/projects/djangotools/helper_scripts/standalone_tools/for_download/mhc_ii/smm_dir.txt", "r")
            smm_dir=FILE.readline().strip()
            FILE = open("/home/dorjee/projects/djangotools/helper_scripts/standalone_tools/for_download/mhc_ii/nn_dir.txt", "r")
            nn_dir=FILE.readline().strip()
            FILE = open("/home/dorjee/projects/djangotools/helper_scripts/standalone_tools/for_download/mhc_ii/net2_dir.txt", "r")
            net2_dir=FILE.readline().strip()
            smm_predictor=SMMPredictor(smm_dir)
            nn_predictor=NNPredictor(nn_dir)
            comb_predictor=COMBLIBMatrix()
            sturniolo_predictor=SturnioloMatrix()
            net2_predictor = Net2Predictor(net2_dir)

            for (allele_list,length) in self.tool_selection:
                alleles = allele_list.split(',')
                for allele in alleles:
                    length = 15
                    
                    if not allele in con_status: 
                        raise UnexpectedInputError("No prediction method consensus3 available for " + allele + ". To see all available alleles, please type command:\n$ python mhc_II_binding.py allele\n")
                    
                    if con_status[allele][1] == "1":
                        try:
                            comb_predictor.pickle_load(self.consensus_tool_location(allele, 9, "comblib") + ".txt")
                        except:
                            raise UnexpectedInputError("No prediction method comblib available for " + allele + ". To see all available alleles, please type command:\n$ python mhc_II_binding.py allele\n")
                        try:
                            comb_perc_file = self.consensus_perc_location(allele, "comblib") + ".txt"
                        except:
                            raise UnexpectedInputError("Selected comb ranking file does not exist.")
                        comb_percs = []
                        comb_data = open(comb_perc_file, 'r').readlines()
                        for i in range(len(comb_data)):
                            comb_percs.append(float(comb_data[i].strip()))
    
                    if con_status[allele][2] == "1":
                        try:
                            smm_perc_file = self.consensus_perc_location(allele, "smm_align") + ".txt"
                        except:
                            raise UnexpectedInputError("Selected smm_align ranking file does not exist.")      
                        smm_percs = []
                        smm_data = open(smm_perc_file, 'r').readlines()
                        for i in range(len(smm_data)):
                            smm_percs.append(float(smm_data[i].strip()))
    
                        try:
                            nn_perc_file = self.consensus_perc_location(allele, "nn_align") + ".txt"
                        except:
                            raise UnexpectedInputError("Selected nn_align ranking file does not exist.")      
                        nn_percs = []
                        nn_data = open(nn_perc_file, 'r').readlines()
                        for i in range(len(nn_data)):
                            nn_percs.append(float(nn_data[i].strip()))
    
                    if con_status[allele][3] == "1":
                        try:
                            sturniolo_predictor.pickle_load(self.consensus_tool_location(allele, 9, "sturniolo") + ".txt")
                        except:
                            raise UnexpectedInputError("No prediction method sturniolo available for " + allele + ". To see all available methods, please type command:\n$ python mhc_II_binding.py allele\n")
                        try:
                            sturniolo_perc_file = self.consensus_perc_location(allele, "sturniolo") + ".txt"
                        except:
                            raise UnexpectedInputError("Selected sturniolo ranking file does not exist.")
                        sturniolo_percs = []
                        sturniolo_data = open(sturniolo_perc_file, 'r').readlines()
                        for i in range(len(sturniolo_data)):
                            sturniolo_percs.append(float(sturniolo_data[i].strip()))
    
                    scores=[]
                    for sequence in proteins.sequences:
                        seq_scores=[]
                        temp_scores = []
                        method_list = []
                        if con_status[allele][1] == "1":
                            comb_scores = comb_predictor.consensus_predict(sequence)
                            temp_scores.append(comb_scores)
                            method_list.append("comb")
                        else: method_list.append("-comb-")
    
                        if con_status[allele][2] == "1":
                            smm_scores = smm_predictor.peptide_predictions(sequence, allele)
                            smmscores = tuple([a for a in smm_scores[0]])
                            temp_scores.append(smmscores)
                            method_list.append("smm")
    
                            nn_scores = nn_predictor.peptide_predictions(sequence, allele)
                            if not nn_scores:
                                method_list.append("-nn-")
                            else: 
                                nnscores = tuple([a for a in nn_scores[0]])
                                method_list.append("nn")
                        else: 
                            blank = "-smm-","-nn-"
                            method_list.extend(blank)
    
                        if con_status[allele][1] == "0" and con_status[allele][3] == "1":
                            sturniolo_scores = sturniolo_predictor.consensus_predict(sequence)
                            temp_scores.append(sturniolo_scores)
                            method_list.append("sturniolo")
                        else: method_list.append("-sturniolo-")
    
                        
                        check_len = 0
                        for check_i in range(len(temp_scores)):
                            check_len += len(temp_scores[check_i]) - len(temp_scores[-1])
    
                        if(check_len != 0):
                            if len(sequence) > 15:
                                raise ValueError("Methods return different number of predictions - %f"  % len(temp_scores[1])) 
                        
                        temp_list = temp_scores[0]
                        for i in range(len(temp_list)):
                            temp_result = []
                            consensus_percs=[]
                            for m_name in method_list: 
                                if m_name == "comb":
                                    (temp_core, temp_score) = comb_scores[i]
                                    temp_result.append(temp_core)
                                    temp_result.append(round(temp_score,2))
                                    perc = 0.01*bisect_left(comb_percs, temp_score)
                                    if perc < 0.01:
                                        perc = 0.01
                                    temp_result.append(perc)
                                    consensus_percs.append(perc)
                                elif m_name == "smm":
                                    (temp_core, temp_score) = smmscores[i]
                                    temp_result.append(temp_core)
                                    temp_result.append(round(temp_score,2))
                                    perc = 0.01*bisect_left(smm_percs, temp_score)
                                    if perc < 0.01:
                                        perc = 0.01
                                    temp_result.append(perc)
                                    consensus_percs.append(perc)
                                elif m_name == "nn":
                                    (temp_core, temp_score) = nnscores[i]
                                    temp_result.append(temp_core)
                                    temp_result.append(round(temp_score,2))
                                    perc = 0.01*bisect_left(nn_percs, temp_score)
                                    if perc < 0.01:
                                        perc = 0.01
                                    temp_result.append(perc)
                                    consensus_percs.append(perc)
                                elif m_name == "sturniolo":
                                    (temp_core, temp_score) = sturniolo_scores[i]
                                    temp_result.append(temp_core)
                                    temp_result.append(round(temp_score,2))
                                    perc = 100 - 0.01*bisect(sturniolo_percs, temp_score)
                                    if perc < 0.01:
                                        perc = 0.01
                                    temp_result.append(perc)
                                    consensus_percs.append(perc)
                                else: 
                                    blank = "-","-","-"
                                    temp_result.extend(blank)
                            consensus_percs.sort()
                            consensus_perc=self.getMedian(consensus_percs)
                            temp_result.insert(0, consensus_perc)
                            seq_scores.append(tuple(temp_result))
                        scores.append(tuple(seq_scores))
                    results.append((length, allele, scores))
                    
        elif self.method=="IEDB_recommended":
            con_status = {}
            con_list = open("/home/dorjee/projects/djangotools/helper_scripts/standalone_tools/for_download/mhc_ii/tool_data/MHCII/IEDB_recommended/recommended_allele_list.txt").read().strip().split("\n")
            for con_element in con_list:
                con_arr = con_element.split("\t")
                con_status[con_arr[0]] = con_arr[1:]
            FILE = open("/home/dorjee/projects/djangotools/helper_scripts/standalone_tools/for_download/mhc_ii/smm_dir.txt", "r")
            smm_dir=FILE.readline().strip()
            FILE = open("/home/dorjee/projects/djangotools/helper_scripts/standalone_tools/for_download/mhc_ii/nn_dir.txt", "r")
            nn_dir=FILE.readline().strip()
            FILE = open("/home/dorjee/projects/djangotools/helper_scripts/standalone_tools/for_download/mhc_ii/net2_dir.txt", "r")
            net2_dir=FILE.readline().strip()
            smm_predictor=SMMPredictor(smm_dir)
            nn_predictor=NNPredictor(nn_dir)
            comb_predictor=COMBLIBMatrix()
            sturniolo_predictor=SturnioloMatrix()
            net2_predictor = Net2Predictor(net2_dir)

            for (allele_list,length) in self.tool_selection:
                alleles = allele_list.split(',')
                for allele in alleles:
                    length = 15
                    
                    if not allele in con_status:
                    	raise UnexpectedInputError("No prediction method IEDB_recommended available for " + allele + ". To see all available alleles, please type command:\n$ python mhc_II_binding.py allele\n")
                    
                    if con_status[allele][0] == "1":
                        try:
                            comb_predictor.pickle_load(self.consensus_tool_location(allele, 9, "comblib") + ".txt")
                        except:
                            raise UnexpectedInputError("No prediction method comblib available for " + allele + ". To see all available methods, please type command:\n$ python mhc_II_binding.py allele\n")
                            
                        try:
                            comb_perc_file = self.consensus_perc_location(allele, "comblib") + ".txt"
                        except:
                            raise UnexpectedInputError("Selected comblib ranking file does not exist.")
                        comb_percs = []
                        comb_data = open(comb_perc_file, 'r').readlines()
                        for i in range(len(comb_data)):
                            comb_percs.append(float(comb_data[i].strip()))
    
                    if con_status[allele][1] == "1":
                        try:
                            smm_perc_file = self.consensus_perc_location(allele, "smm_align") + ".txt"
                        except:
                            raise UnexpectedInputError("Selected smm_align ranking file does not exist.")      
                        smm_percs = []
                        smm_data = open(smm_perc_file, 'r').readlines()
                        for i in range(len(smm_data)):
                            smm_percs.append(float(smm_data[i].strip()))
    
                        try:
                            nn_perc_file = self.consensus_perc_location(allele, "nn_align") + ".txt"
                        except:
                            raise UnexpectedInputError("Selected nn_align ranking file does not exist.")      
                        nn_percs = []
                        nn_data = open(nn_perc_file, 'r').readlines()
                        for i in range(len(nn_data)):
                            nn_percs.append(float(nn_data[i].strip()))
    
                    if con_status[allele][2] == "1":
                        try:
                            sturniolo_predictor.pickle_load(self.consensus_tool_location(allele, 9, "sturniolo") + ".txt")
                        except:
                            raise UnexpectedInputError("No prediction method sturniolo available for " + allele + ". To see all available methods, please type command:\n$ python mhc_II_binding.py allele\n")
                            
                        try:
                            sturniolo_perc_file = self.consensus_perc_location(allele, "sturniolo") + ".txt"
                        except:
                            raise UnexpectedInputError("Selected sturniolo tanking file does not exist.")
                        sturniolo_percs = []
                        sturniolo_data = open(sturniolo_perc_file, 'r').readlines()
                        for i in range(len(sturniolo_data)):
                            sturniolo_percs.append(float(sturniolo_data[i].strip()))
    
                    if con_status[allele][3] == "1":
                        try:
                            net2_perc_file = self.consensus_perc_location(allele, "netmhc2pan") + ".txt"
                        except:
                            raise UnexpectedInputError("You are currently using nn_align ranking file. May be selected netmhc2pan ranking file does not exist.")
                        net2_percs = []
                        net2_data = open(net2_perc_file, 'r').readlines()
                        for i in range(len(net2_data)):
                            net2_percs.append(float(net2_data[i].strip()))
    
                    scores=[]
                    for sequence in proteins.sequences:
                        seq_scores=[]
                        temp_scores = []
                        method_list = []
                        if con_status[allele][0] == "1":
                            comb_scores = comb_predictor.consensus_predict(sequence)
                            temp_scores.append(comb_scores)
                            method_list.append("comb.lib.")
                        else: method_list.append("-")
    
                        if con_status[allele][1] == "1":
                            smm_scores = smm_predictor.peptide_predictions(sequence, allele)
                            smmscores = tuple([a for a in smm_scores[0]])
                            temp_scores.append(smmscores)
                            method_list.append("smm")
    
                            nn_scores = nn_predictor.peptide_predictions(sequence, allele)
                            if not nn_scores:
                                method_list.append("-")
                            else: 
                                nnscores = tuple([a for a in nn_scores[0]])
                                method_list.append("nn")
                        else: 
                            blank = "-","-"
                            method_list.extend(blank)
    
                        if con_status[allele][0] == "0" and con_status[allele][1] == "0" and con_status[allele][2] == "0" and con_status[allele][3] == "1": 
                            net2_scores = net2_predictor.peptide_predictions(sequence, allele)
                            net2scores = tuple([a for a in net2_scores[0]])
                            temp_scores.append(net2scores)
                            method_list.append("NetMHCIIpan")
                        else: method_list.append("-")
    
                        if con_status[allele][0] == "0" and con_status[allele][2] == "1":
                            sturniolo_scores = sturniolo_predictor.consensus_predict(sequence)
                            temp_scores.append(sturniolo_scores)
                            method_list.append("sturniolo")
                        else: method_list.append("-")
    
                        
                        check_len = 0
                        for check_i in range(len(temp_scores)):
                            check_len += len(temp_scores[check_i]) - len(temp_scores[-1])
    
                        if(check_len != 0):
                            if len(sequence) > 15:
                                raise ValueError("Methods return different number of predictions - %f"  % len(temp_scores[1])) 
                        
                        temp_list = temp_scores[0]
                        for i in range(len(temp_list)):
                            temp_result = []
                            consensus_percs=[]
                            method_used=[]
                            for m_name in method_list: 
                                if m_name == "comb.lib.":
                                    (temp_core, temp_score) = comb_scores[i]
                                    temp_result.append(temp_core)
                                    temp_result.append(round(temp_score,2))
                                    perc = 0.01*bisect_left(comb_percs, temp_score)
                                    if perc < 0.01:
                                        perc = 0.01
                                    temp_result.append(perc)
                                    consensus_percs.append(perc)
                                elif m_name == "smm":
                                    (temp_core, temp_score) = smmscores[i]
                                    temp_result.append(temp_core)
                                    temp_result.append(round(temp_score,2))
                                    perc = 0.01*bisect_left(smm_percs, temp_score)
                                    if perc < 0.01:
                                        perc = 0.01
                                    temp_result.append(perc)
                                    consensus_percs.append(perc)
                                elif m_name == "nn":
                                    (temp_core, temp_score) = nnscores[i]
                                    temp_result.append(temp_core)
                                    temp_result.append(round(temp_score,2))
                                    perc = 0.01*bisect_left(nn_percs, temp_score)
                                    if perc < 0.01:
                                        perc = 0.01
                                    temp_result.append(perc)
                                    consensus_percs.append(perc)
                                elif m_name == "NetMHCIIpan":
                                    (temp_core, temp_score) = net2scores[i]
                                    temp_result.append(temp_core)
                                    temp_result.append(round(temp_score,2))
                                    perc = 0.01*bisect_left(net2_percs, temp_score)
                                    if perc < 0.01:
                                        perc = 0.01
                                    temp_result.append(perc)
                                    consensus_percs.append(perc) 
                                elif m_name == "sturniolo":
                                    (temp_core, temp_score) = sturniolo_scores[i]
                                    temp_result.append(temp_core)
                                    temp_result.append(round(temp_score,2))
                                    perc = 100 - 0.01*bisect(sturniolo_percs, temp_score)
                                    if perc < 0.01:
                                        perc = 0.01
                                    temp_result.append(perc)
                                    consensus_percs.append(perc)
                                else: 
                                    blank = "-","-","-"
                                    temp_result.extend(blank)
                                found = re.search(r'\-',m_name)
                                
                                if not found:
                                    method_used.append(m_name)
                            consensus_percs.sort()
                            consensus_perc=self.getMedian(consensus_percs)
                            temp_result.insert(0, consensus_perc)
                            method_used = ",".join([i for i in method_used])
                            temp_result.append(method_used)
                            seq_scores.append(tuple(temp_result))
                        scores.append(tuple(seq_scores))
                    results.append((length, allele, scores))
                    
        elif self.method=="smm_align":
            FILE = open("/home/dorjee/projects/djangotools/helper_scripts/standalone_tools/for_download/mhc_ii/smm_dir.txt", "r")
            smm_dir=FILE.readline().strip()
            predictor=SMMPredictor(smm_dir)
            for (allele,length) in self.tool_selection:
                length = int(length)
            try:
                scores=[]
                for sequence in proteins.sequences:
                    list_of_lists = predictor.peptide_predictions(sequence, allele)
                    for lis in list_of_lists:
                        scores.append(tuple(lis))
            except:
                raise UnexpectedInputError("No prediction method smm_align available for " + allele + ". To see all available methods, please type command:\n$ python mhc_II_binding.py allele\n")
                   
            alleles=allele.split(',')
            for i, a in enumerate(alleles):
                results.append((length, a, scores))
                
        elif self.method=="nn_align":
            FILE = open("/home/dorjee/projects/djangotools/helper_scripts/standalone_tools/for_download/mhc_ii/nn_dir.txt", "r")
            nn_dir=FILE.readline().strip()
            predictor=NNPredictor(nn_dir)
            for (allele,length) in self.tool_selection:
                length = int(length)
            try:
                scores=[]
                for sequence in proteins.sequences:
                    list_of_lists = predictor.peptide_predictions(sequence, allele)
                    for lis in list_of_lists:
                        scores.append(tuple(lis))
            except:
                raise UnexpectedInputError("No prediction method nn_align available for " + allele + ". To see all available methods, please type command:\n$ python mhc_II_binding.py allele\n")
            
            alleles=allele.split(',')
            for i, a in enumerate(alleles):
                results.append((length, a, scores))

        elif self.method == "comblib":
            predictor=COMBLIBMatrix()
            for (alleles, length) in self.tool_selection:
                allele_list = alleles.split(",")
                for allele in allele_list: 
                    try:
                        predictor.pickle_load(self.comblib_tool_location(allele, length) + ".txt")
                    except:
                        raise UnexpectedInputError("No prediction method comblib available for " + allele + ". To see all available methods, please type command:\n$ python mhc_II_binding.py allele\n")
                    scores=[]
                    for sequence in proteins.sequences:
                        scores.append(predictor.consensus_predict(sequence))
                    results.append((length, allele, scores))
                    
        elif self.method=="NetMHCIIpan":
            FILE = open("/home/dorjee/projects/djangotools/helper_scripts/standalone_tools/for_download/mhc_ii/net2_dir.txt", "r")
            net_dir=FILE.readline().strip()
            predictor=Net2Predictor(net_dir)
            for (allele,length) in self.tool_selection:
                length = int(length)
                try:
                    scores=[]
                    for sequence in proteins.sequences:
                        list_of_lists = predictor.peptide_predictions(sequence, allele)
                        for lis in list_of_lists:
                            scores.append(tuple(lis))
                except:
                    raise UnexpectedInputError("No prediction method NetMHCIIpan available for " + allele + ". To see all available methods, please type command:\n$ python mhc_II_binding.py allele\n")
                    
            alleles=allele.split(',')
            for i, a in enumerate(alleles):
                results.append((length, a, scores))
                
        else:
            if self.method=="sturniolo":
                predictor=SturnioloMatrix()
            else:
                raise UnexpectedInputError("Selected prediction method '%s' does not exist.\nTo see all available methods, please type command:\n$ python mhc_II_binding.py method\n" % self.method)
                
            for (allele, length) in self.tool_selection:
                try:
                    predictor.pickle_load(self.tool_location(allele, length) + ".txt")
                except:
		    print "hello"
                    raise UnexpectedInputError("No prediction method sturniolo available for " + allele + ". To see all available methods, please type command:\n$ python mhc_II_binding.py allele\n")
                scores=[]
                for sequence in proteins.sequences:
                    scores.append(predictor.consensus_predict(sequence))
                results.append((length, allele, scores))
                        
        return results

def print_input_page(sequence):
        config_parser = ConfigParser.ConfigParser()
        config_parser.read("../setup.cfg")
        html_path=config_parser.get("path", "html")
        template = open(html_path + "/html/mhc_II_binding.seq","r").read()
        print "Content-Type: text/html"
        print ""
        print template % sequence


def print_result_page(title, body="", head=""):
    print "Content-Type: text/html"
    print ""
    print """<html xmlns="http://www.w3.org/1999/xhtml" xml:lang="en" lang="en">
    <head>
        <title>MHC-II binding predictions</title>
        <link rel="stylesheet" type="text/css" href="/files/css/iedb.css">
        %s
    </head>
    <body id="tab1">
        <a href = "/" id = "header1">IEDB Analysis Resource</a>
        <ul id="tabnav">
            <li class="tab1"><a href = "../html/mhc_II_binding.html">Home</a></li>
            <li class="tab2"><a href = "../html/tutorial_mhc_II_binding.html">Help</a></li>
            <li class="tab3"><a href = "../html/example_mhc_II_binding.html">Example</a></li>
	    <li class="tab4"><a href = "../html/reference_mhc_II_binding.html">Reference</a></li>


        </ul>
        <div id="center">
            <h2>%s<h2>
            %s
        </div>
    </body></html>""" % (head, title, body)

def get_old_name(new_allele):
    ''' gets a list of old HLA alleles, given the new ones '''
    reader = csv.reader(open("/home/dorjee/projects/djangotools/helper_scripts/standalone_tools/for_download/mhc_ii/tool_data/MHCII/HLA_nomenclature_table"), delimiter='\t')
    row_dict = dict((row[0], row[1]) for row in reader)
    if row_dict.has_key(new_allele):
        return row_dict[new_allele]
    else:
        return new_allele

def get_new_name(old_allele):
    ''' gets the new HLA allele, given the old one '''
    if re.search("HLA.*", old_allele):
        reader = csv.reader(open("/home/dorjee/projects/djangotools/helper_scripts/standalone_tools/for_download/mhc_ii/tool_data/MHCII/HLA_nomenclature_table"), delimiter='\t')
        row_dict = dict((row[0], row[1]) for row in reader)
        for key, value in dict.items(row_dict):
            if old_allele == value:
                return key
    else:
        return old_allele
    

def main(form):
    if form.has_key("sequence"):
        print_input_page(form["sequence"])
    elif not form.has_key("sequence_format"):
        print_input_page("")
    else:
        try:
            proteins=Proteins(form)
            mhc_predictor=MHCBindingPredictions(form)
            table = PredictionTable(form)
            mhc_scores = mhc_predictor.predict(proteins)
            table_rows = table.format_binding(proteins, mhc_scores)
            output_format= form["output_format"]
            if output_format!="xhtml" and output_format!="ascii":
                raise UnexpectedInputError("Invalid output format '" + output_format + "'")
            con_status = {}
            con_list = open("./tool_data/MHCII/consensus/consensus_allele_list.txt").read().strip().split("\n")
            for con_element in con_list:
                con_arr = con_element.split("\t")
                con_status[con_arr[0]] = con_arr[1:]
        except InputError, inst:
            print_result_page("Input Error", """<p>%s</p>""" % inst)

        except Exception, inst:
            sys.exit(inst)
        else:
            too_long_flag=False
            if len(table_rows)>10000:
                output_format="ascii"
                too_long_flag=True

            if output_format=="ascii":
                if form["pred_method"] == "consensus3":
                    print "allele\tseq_num\tstart\tend\tpeptide\tconsensus_percentile_rank\tcomblib_core\tcomblib_score\tcomblib_rank\tsmm_align_core\tsmm_align_ic50\tsmm_align_rank\tnn_align_core\tnn_align_ic50\tnn_align_rank\tsturniolo_score\tsturniolo_rank"
                elif form["pred_method"] == "IEDB_recommended":
                    print "allele\tseq_num\tstart\tend\tpeptide\tmethod\tpercentile_rank\tcomblib_core\tcomblib_score\tcomblib_rank\tsmm_align_core\tsmm_align_ic50\tsmm_align_rank\tnn_align_core\tnn_align_ic50\tnn_align_rank\tnetmhciipan_core\tnetmhciipan_ic50\tnetmhciipan_rank\tsturniolo_core\tsturniolo_score\tsturniolo_rank"
                elif form["pred_method"] == "sturniolo":
                    print "allele\tseq_num\tstart\tend\tcore_peptide\tpeptide\tscore"
                else:
                    print "allele\tseq_num\tstart\tend\tcore_peptide\tpeptide\tic50"
                
                for table_row in table_rows:
                    table_row = list(table_row)

                    seq_pos = table_row[1].split(':')
                    pos = seq_pos[1].split('-')
                    table_row[1] = seq_pos[0]
                    table_row.insert(2,pos[0])
                    table_row.insert(3,pos[1])
                    
                    if 'H2' not in table_row[0]:
                        table_row[0] = "HLA-%s" % table_row[0].replace("-", "/")
                        
                    if form["pred_method"] == "consensus3":
                        print """%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s""" % tuple(table_row)
                    elif form["pred_method"] == "IEDB_recommended":
                        method_used = table_row[-1].upper().replace("STURNIOLO","Sturniolo").replace('NETMHCIIPAN','NetMHCIIpan')
                        method_list = method_used.split(',')
                        if len(method_list) <= 1:  table_row.insert(5,method_used)
                        else:
                            method_used = "Consensus("+method_used+")"
                            table_row.insert(5,method_used)
                        print """%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s""" % tuple(table_row[:-1])
                    elif form["pred_method"] == "sturniolo":
                        print """%s\t%s\t%s\t%s\t%s\t%s\t%f""" % tuple(table_row)
                    else:
                        print """%s\t%s\t%s\t%s\t%s\t%s\t%f""" % tuple(table_row)

def commandline_help():
    print " ________________________________________________________________________________________"
    print "|****************************************************************************************|"
    print "| * List all available commands.                                                         |"
    print "| python mhc_II_binding.py                                                               |"
    print "|________________________________________________________________________________________|"
    print "| * List all available mhc_II prediction methods.                                        |"
    print "| python mhc_II_binding.py method                                                        |"
    print "|________________________________________________________________________________________|"
    print "| * List all alleles.                                                                    |"
    print "| python mhc_II_binding.py allele                                                        |"
    print "|________________________________________________________________________________________|"
    print "| * Make predictions given a file containing a list of sequences.                        |"
    print "| python mhc_II_binding.py prediction_method_name allele_name input_sequence_file_name   |"
    print "| Example: python mhc_II_binding.py consensus3 HLA-DRB1*03:01 test.fasta                 |"
    print "|________________________________________________________________________________________|"
    print "| * You may also redirect (pipe) the input file into the script.                         |"
    print "| Example: echo -e test.fasta | python mhc_II_binding.py consensus3 HLA-DRB1*03:01       |"
    print "|________________________________________________________________________________________|"

def commandline_method():
    '''Return all available prediction methods.'''
    print 
    print "Available methods are:"
    print "----------------------"
    print "comblib"
    print "consensus3"
    print "IEDB_recommended"
    print "NetMHCIIpan"
    print "nn_align"
    print "smm_align"
    print "sturniolo"
    print
    
def commandline_allele():
    '''Return all available alleles.'''
    print """
  ================================================================================================================ 
 | Alleles available for consensus3 method:                                                                       | 
 |----------------------------------------------------------------------------------------------------------------| 
 | HLA-DPA1*01/DPB1*04:01    | HLA-DRB1*01:01 | HLA-DRB1*04:04 | HLA-DRB1*08:04 | HLA-DRB1*11:21 | HLA-DRB1*13:27 | 
 | HLA-DPA1*01:03/DPB1*02:01 | HLA-DRB1*01:02 | HLA-DRB1*04:05 | HLA-DRB1*08:06 | HLA-DRB1*11:28 | HLA-DRB1*13:28 | 
 | HLA-DPA1*02:01/DPB1*01:01 | HLA-DRB1*03:01 | HLA-DRB1*04:08 | HLA-DRB1*08:13 | HLA-DRB1*13:01 | HLA-DRB1*15:01 | 
 | HLA-DPA1*02:01/DPB1*05:01 | HLA-DRB1*03:05 | HLA-DRB1*04:10 | HLA-DRB1*08:17 | HLA-DRB1*13:02 | HLA-DRB1*15:02 | 
 | HLA-DPA1*03:01/DPB1*04:02 | HLA-DRB1*03:06 | HLA-DRB1*04:21 | HLA-DRB1*11:01 | HLA-DRB1*13:04 | HLA-DRB1*15:06 | 
 | HLA-DQA1*01:01/DQB1*05:01 | HLA-DRB1*03:07 | HLA-DRB1*04:23 | HLA-DRB1*11:02 | HLA-DRB1*13:05 | HLA-DRB3*01:01 | 
 | HLA-DQA1*01:02/DQB1*06:02 | HLA-DRB1*03:08 | HLA-DRB1*04:26 | HLA-DRB1*11:04 | HLA-DRB1*13:07 | HLA-DRB4*01:01 | 
 | HLA-DQA1*03:01/DQB1*03:02 | HLA-DRB1*03:09 | HLA-DRB1*07:01 | HLA-DRB1*11:06 | HLA-DRB1*13:11 | HLA-DRB5*01:01 | 
 | HLA-DQA1*04:01/DQB1*04:02 | HLA-DRB1*03:11 | HLA-DRB1*07:03 | HLA-DRB1*11:07 | HLA-DRB1*13:21 | HLA-DRB5*01:05 | 
 | HLA-DQA1*05:01/DQB1*02:01 | HLA-DRB1*04:01 | HLA-DRB1*08:01 | HLA-DRB1*11:14 | HLA-DRB1*13:22 | H2-IAb         | 
 | HLA-DQA1*05:01/DQB1*03:01 | HLA-DRB1*04:02 | HLA-DRB1*08:02 | HLA-DRB1*11:20 | HLA-DRB1*13:23 | H2-IAd         |
 |----------------------------------------------------------------------------------------------------------------| 
 | Alleles available for comblib method:                                                                          |  
 |----------------------------------------------------------------------------------------------------------------| 
 | HLA-DPA1*01/DPB1*04:01    | HLA-DRB1*01:01 |                                                                   |                                                             | 
 | HLA-DPA1*01:03/DPB1*02:01 | HLA-DRB1*07:01 |                                                                   |               
 | HLA-DPA1*02:01/DPB1*01:01 | HLA-DRB1*09:01 |                                                                   | 
 | HLA-DPA1*02:01/DPB1*05:01 | HLA-DRB3*01:01 |                                                                   | 
 | HLA-DPA1*03:01/DPB1*04:02 | HLA-DRB4*01:01 |                                                                   |               
 | HLA-DQA1*01:01/DQB1*05:01 |                |                                                                   | 
 | HLA-DQA1*01:02/DQB1*06:02 |                |                                                                   | 
 | HLA-DQA1*03:01/DQB1*03:02 |                |                                                                   |  
 | HLA-DQA1*04:01/DQB1*04:02 |                |                                                                   | 
 | HLA-DQA1*05:01/DQB1*02:01 |                |                                                                   | 
 | HLA-DQA1*05:01/DQB1*03:01 |                |                                                                   | 
 |----------------------------------------------------------------------------------------------------------------| 
 | Alleles available for smm_align method:                                                                        | 
 |----------------------------------------------------------------------------------------------------------------| 
 | HLA-DPA1*01/DPB1*04:01    | HLA-DRB1*01:01 | HLA-DRB1*15:01 |                                                  | 
 | HLA-DPA1*01:03/DPB1*02:01 | HLA-DRB1*03:01 | HLA-DRB3*01:01 |                                                  | 
 | HLA-DPA1*02:01/DPB1*01:01 | HLA-DRB1*04:01 | HLA-DRB4*01:01 |                                                  | 
 | HLA-DPA1*02:01/DPB1*05:01 | HLA-DRB1*04:04 | HLA-DRB5*01:01 |                                                  | 
 | HLA-DPA1*03:01/DPB1*04:02 | HLA-DRB1*04:05 | H2-IAb         |                                                  | 
 | HLA-DQA1*01:01/DQB1*05:01 | HLA-DRB1*07:01 | H2-IAd         |                                                  | 
 | HLA-DQA1*01:02/DQB1*06:02 | HLA-DRB1*08:02 | H2-IEd         |                                                  | 
 | HLA-DQA1*03:01/DQB1*03:02 | HLA-DRB1*09:01 |                |                                                  | 
 | HLA-DQA1*04:01/DQB1*04:02 | HLA-DRB1*11:01 |                |                                                  | 
 | HLA-DQA1*05:01/DQB1*02:01 | HLA-DRB1*12:01 |                |                                                  | 
 | HLA-DQA1*05:01/DQB1*03:01 | HLA-DRB1*13:02 |                |                                                  | 
 |----------------------------------------------------------------------------------------------------------------| 
 | Alleles available for nn_align method:                                                                         | 
 |----------------------------------------------------------------------------------------------------------------| 
 | HLA-DPA1*01/DPB1*04:01    | HLA-DRB1*01:01 | HLA-DRB1*15:01 |                                                  | 
 | HLA-DPA1*01:03/DPB1*02:01 | HLA-DRB1*03:01 | HLA-DRB3*01:01 |                                                  | 
 | HLA-DPA1*02:01/DPB1*01:01 | HLA-DRB1*04:01 | HLA-DRB4*01:01 |                                                  | 
 | HLA-DPA1*02:01/DPB1*05:01 | HLA-DRB1*04:04 | HLA-DRB5*01:01 |                                                  | 
 | HLA-DPA1*03:01/DPB1*04:02 | HLA-DRB1*04:05 | H2-IAb         |                                                  | 
 | HLA-DQA1*01:01/DQB1*05:01 | HLA-DRB1*07:01 | H2-IAd         |                                                  | 
 | HLA-DQA1*01:02/DQB1*06:02 | HLA-DRB1*08:02 |                |                                                  | 
 | HLA-DQA1*03:01/DQB1*03:02 | HLA-DRB1*09:01 |                |                                                  | 
 | HLA-DQA1*04:01/DQB1*04:02 | HLA-DRB1*11:01 |                |                                                  | 
 | HLA-DQA1*05:01/DQB1*02:01 | HLA-DRB1*12:01 |                |                                                  | 
 | HLA-DQA1*05:01/DQB1*03:01 | HLA-DRB1*13:02 |                |                                                  | 
 |----------------------------------------------------------------------------------------------------------------| 
 | Alleles available for sturniolo method:                                                                         | 
 |----------------------------------------------------------------------------------------------------------------| 
 | HLA-DRB1*01:01 | HLA-DRB1*04:21 | HLA-DRB1*11:07 | HLA-DRB1*13:28 |                                            | 
 | HLA-DRB1*01:02 | HLA-DRB1*04:23 | HLA-DRB1*11:14 | HLA-DRB1*15:01 |                                            | 
 | HLA-DRB1*03:01 | HLA-DRB1*04:26 | HLA-DRB1*11:20 | HLA-DRB1*15:02 |                                            | 
 | HLA-DRB1*03:05 | HLA-DRB1*07:01 | HLA-DRB1*11:21 | HLA-DRB1*15:06 |                                            | 
 | HLA-DRB1*03:06 | HLA-DRB1*07:03 | HLA-DRB1*11:28 | HLA-DRB5*01:01 |                                            | 
 | HLA-DRB1*03:07 | HLA-DRB1*08:01 | HLA-DRB1*13:01 | HLA-DRB5*01:05 |                                            | 
 | HLA-DRB1*03:08 | HLA-DRB1*08:02 | HLA-DRB1*13:02 |                |                                            | 
 | HLA-DRB1*03:09 | HLA-DRB1*08:04 | HLA-DRB1*13:04 |                |                                            | 
 | HLA-DRB1*03:11 | HLA-DRB1*08:06 | HLA-DRB1*13:05 |                |                                            | 
 | HLA-DRB1*04:01 | HLA-DRB1*08:13 | HLA-DRB1*13:07 |                |                                            | 
 | HLA-DRB1*04:02 | HLA-DRB1*08:17 | HLA-DRB1*13:11 |                |                                            | 
 | HLA-DRB1*04:04 | HLA-DRB1*11:01 | HLA-DRB1*13:21 |                |                                            | 
 | HLA-DRB1*04:05 | HLA-DRB1*11:02 | HLA-DRB1*13:22 |                |                                            | 
 | HLA-DRB1*04:08 | HLA-DRB1*11:04 | HLA-DRB1*13:23 |                |                                            | 
 | HLA-DRB1*04:10 | HLA-DRB1*11:06 | HLA-DRB1*13:27 |                |                                            | 
 |----------------------------------------------------------------------------------------------------------------| 
 | Alleles available for NetMHCIIpan method:                                                                      | 
 |----------------------------------------------------------------------------------------------------------------| 
 |  Please see the NetMHCIIpan webserver for information on the                                                   | 
 |  more than 500 alleles available:                                                                              | 
 |  http://www.cbs.dtu.dk/services/NetMHCIIpan/                                                                   | 
  ----------------------------------------------------------------------------------------------------------------
"""

if __name__ == '__main__':
    import select
    debug = False
    
    method = allele = infile = None
    
    if (len(sys.argv) == 1):
        commandline_help()
    elif ((len(sys.argv) == 2) and (sys.argv[1] == "method")):
        commandline_method()
    elif ((len(sys.argv) == 2) and (sys.argv[1] == "allele")):
        commandline_allele()
    else:
        method = sys.argv[1]
        if len(sys.argv) > 2:
            allele = sys.argv[2]
    
    # If there's input ready, do something, else do something
    # else. Note timeout is zero so select won't block at all.
    if sys.stdin in select.select([sys.stdin], [], [], 0)[0]:
        infile = sys.stdin.readline().strip()
        if infile: pass
        else: # an empty line means stdin has been closed
            exit(0)
    else:
        if len(sys.argv) > 3:
            infile = sys.argv[3]
            
    if not method or not allele or not infile:
        print "* To run the predction, you must specify method, allele and input file name."
        exit(0)
    
    allele = allele.replace("_","-").replace("H-2","H2")
    seq = [('sequence_format', 'auto'), ('sort_output', 'position_in_sequence'), ('cutoff_type', 'none'), ('output_format', 'ascii'), ('allele', allele), ('sequence_file', infile), ('pred_method', method)]
    form = dict(seq)
    if not debug:
        main(form)
    else:
        calltime = str(time.ctime(time.time()))
        cgilog = open("tool_data/tool_log.txt","a")
        cgilog.write('\n--- '+ calltime +' ---\n')
        cgilog.write('mhc_II_binding.py\n')
        cgilog.flush()
        cgilog.write('Form keys: %d\n' % len(form.keys()))
        for key in form.keys():
            cgilog.write('%s - "%s" - "%s"\n' % (calltime,key, form[key]))
        cgilog.close()

        try:
            main(form)
        except Exception, inst:
            cgilog = open("tool_data/tool_log.txt","a")
            cgilog.write('XXX '+ calltime +' XXX\n')
            cgilog.write("EXCEPTION: '%s'\n" % str(inst))
            traceback.print_exc(None, cgilog)
            cgilog.close()
        else:
            cgilog = open("tool_data/tool_log.txt","a")
            cgilog.write('... '+ calltime +' ...\n')
            cgilog.close()
            
