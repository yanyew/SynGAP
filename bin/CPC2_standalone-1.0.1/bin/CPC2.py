#!/usr/bin/env python

# This version is edited at 20190124 for compatibility on OS system

import sys
import os
import re
if (sys.version_info[0] == 2):
	import commands
else:
	import subprocess
import time
from optparse import OptionParser,OptionGroup
import six
import numpy as np
from Bio.Seq import Seq
from Bio.SeqUtils import ProtParam
import seqio

def __main():
	start_time = time.time()
	usage = "usage: %prog [options] -i input.fasta -o output_file"
	description = "Contact: Kang Yujian <kangyj@mail.cbi.pku.edu.cn>"
	parser = OptionParser(usage,version="%prog 0.1",description = description)
	Common_group = OptionGroup(parser,"Common Options")
	Common_group.add_option("-i",dest="fasta",help="input sequence in fasta format [Required]",metavar="FILE",type="string",default=None)
	Common_group.add_option("-o",dest="outfile",help="output file [Default: cpc2output.txt]",metavar="FILE",type="string",default="cpc2output.txt")
	Common_group.add_option("-r",dest="reverse",help="also check the reverse strand [Default: FALSE]",action="store_true")
	Common_group.add_option("--ORF",dest="ORF",help="output the start position of longest ORF [Default: FALSE]",action="store_true")
	parser.add_option_group(Common_group)
	(options, args) = parser.parse_args()
	if options.fasta == None:
		parser.print_help()
		return -1
	else:
		if not os.path.isfile(options.fasta):
			sys.stderr.write("[ERROR] %s is not a file\n"%options.fasta)
			return -1
	if options.reverse:
		strand = "-"
	else:
		strand = "+"
	if options.ORF:
		output_orf = 1
	else:
		output_orf = 0
	if calculate_potential(options.fasta,strand,output_orf,options.outfile):
		return 1
	sys.stderr.write("[INFO] cost time: %ds\n"%(time.time()-start_time))
	return 0

class FindCDS:
	'''
	Find the most like CDS in a given sequence 
	The most like CDS is the longest ORF found in the sequence
	When having same length, the upstream ORF is printed
	modified from source code of CPAT 1.2.1 downloaded from https://sourceforge.net/projects/rna-cpat/files/?source=navbar
	'''
	def __init__(self,seq):
		self.seq = seq
		self.result = (0,0,0,0,0)
		self.longest = 0
		self.basepair = {"A":"T","T":"A","U":"A","C":"G","G":"C","N":"N","X":"X"}

	def _reversecompliment(self):
		return "".join(self.basepair[base] for base in self.seq)[::-1]

	def get_codons(self,frame_number):
		'''
		Record every nucleotide triplet and its coordinate position for input sequence in one frame
		'''
		coordinate = frame_number
		while coordinate + 3 <= len(self.seq):
			yield (self.seq[coordinate:coordinate+3], coordinate)
			coordinate += 3 
	
	def find_longest_in_one(self,myframe,direction,start_codon,stop_codon):
		'''
		find the longest ORF in one reading myframe
		'''
		triplet_got = self.get_codons(myframe)	
		starts = start_codon
		stops = stop_codon
		'''
		Extend sequence by triplet after start codon encountered
		End ORF extension when stop codon encountered
		'''
		while True:
			try: 
				codon,index = next(triplet_got)
			except StopIteration:
				break 
			if codon in starts and codon not in stops:
				'''
				find the ORF start
				'''
				orf_start = index
				end_extension = False
				while True:
					try: 
						codon,index = next(triplet_got)
					except StopIteration:
						end_extension = True
						integrity = -1
					if codon in stops:
						integrity = 1
						end_extension = True
					if end_extension:
						orf_end = index + 3
						Length = (orf_end - orf_start)
						if Length > self.longest:
							self.longest = Length
							self.result = [direction,orf_start,orf_end,Length,integrity]
						if Length == self.longest and orf_start < self.result[1]:
							'''
							if ORFs have same length, return the one that if upstream
							'''
							self.result = [direction,orf_start,orf_end,Length,integrity]
						break

	def longest_orf(self,direction,start_codon={"ATG":None}, stop_codon={"TAG":None,"TAA":None,"TGA":None}):
		return_orf = ""
		for frame in range(3):
			self.find_longest_in_one(frame,"+",start_codon,stop_codon)
		return_orf = self.seq[self.result[1]:self.result[2]][:]
		start_coordinate = self.result[1]
		strand_direction = "+"
		orf_integrity = self.result[4]
		'''
		Also check reverse chain if -r is chosen
		'''
		if direction == "-":
			self.seq = self._reversecompliment()
			for frame in range(3):
				self.find_longest_in_one(frame,"-",start_codon,stop_codon)
			if self.result[0] == "-":
				return_orf = self.seq[self.result[1]:self.result[2]][:]
				start_coordinate = self.result[1]
				strand_direction = "-"
				orf_integrity = self.result[4]
		return return_orf,start_coordinate,strand_direction,orf_integrity


class Fickett:
	'''
	calculate Fickett TESTCODE for full sequence
	NAR 10(17) 5303-531
	modified from source code of CPAT 1.2.1 downloaded from https://sourceforge.net/projects/rna-cpat/files/?source=navbar 
	'''
	def __init__(self):
		'''new compiled Fickett look-up table'''
		self.position_parameter  = [1.9,1.8,1.7,1.6,1.5,1.4,1.3,1.2,1.1,0.0]
		self.content_parameter  = [0.33,0.31,0.29,0.27,0.25,0.23,0.21,0.19,0.17,0]
		self.position_probability = {
			"A":[0.51,0.55,0.57,0.52,0.48,0.58,0.57,0.54,0.50,0.36],
			"C":[0.29,0.44,0.55,0.49,0.52,0.60,0.60,0.56,0.51,0.38],
			"G":[0.62,0.67,0.74,0.65,0.61,0.62,0.52,0.41,0.31,0.17],
	        "T":[0.51,0.60,0.69,0.64,0.62,0.67,0.58,0.48,0.39,0.24],
			}
		self.position_weight = {"A":0.062,"C":0.093,"G":0.205,"T":0.154}
		self.content_probability = {
			"A":[0.40,0.55,0.58,0.58,0.52,0.48,0.45,0.45,0.38,0.19],
			"C":[0.50,0.63,0.59,0.50,0.46,0.45,0.47,0.56,0.59,0.33],
			"G":[0.21,0.40,0.47,0.50,0.52,0.56,0.57,0.52,0.44,0.23],
			"T":[0.30,0.49,0.56,0.53,0.48,0.48,0.52,0.57,0.60,0.51]
			}
		self.content_weight = {"A":0.084,"C":0.076,"G":0.081,"T":0.055}


	def look_up_position_probability(self,value, base):
		'''
		look up positional probability by base and value
		'''
		if float(value) < 0:
			return None
		for idx,val in enumerate (self.position_parameter):
			if (float(value) >= val):
				return float(self.position_probability[base][idx]) * float(self.position_weight[base])

	def look_up_content_probability(self,value, base):
		'''
		look up content probability by base and value
		'''
		if float(value) < 0:
			return None
		for idx,val in enumerate (self.content_parameter):
			if (float(value) >= val):
				return float(self.content_probability[base][idx]) * float(self.content_weight[base])

	def fickett_value(self,dna):
		'''
		calculate Fickett value from full RNA transcript sequence
		'''
		if len(dna) < 2:
			return 0
		fickett_score=0
		dna=dna
		total_base = len(dna)
		A_content = float(dna.count("A"))/total_base
		C_content = float(dna.count("C"))/total_base
		G_content = float(dna.count("G"))/total_base
		T_content = float(dna.count("T"))/total_base

		phase_0 = dna[::3]
		phase_1 = dna[1::3]
		phase_2 = dna[2::3]
		
		phase_0_A = phase_0.count("A")
		phase_1_A = phase_1.count("A")
		phase_2_A = phase_2.count("A")
		phase_0_C = phase_0.count("C")
		phase_1_C = phase_1.count("C")
		phase_2_C = phase_2.count("C")
		phase_0_G = phase_0.count("G")
		phase_1_G = phase_1.count("G")
		phase_2_G = phase_2.count("G")
		phase_0_T = phase_0.count("T")
		phase_1_T = phase_1.count("T")
		phase_2_T = phase_2.count("T")

		A_content = float(phase_0_A + phase_1_A + phase_2_A)/total_base
		C_content = float(phase_0_C + phase_1_C + phase_2_C)/total_base
		G_content = float(phase_0_G + phase_1_G + phase_2_G)/total_base
		T_content = float(phase_0_T + phase_1_T + phase_2_T)/total_base
		A_position= np.max([phase_0_A,phase_1_A,phase_2_A])/(np.min([phase_0_A,phase_1_A,phase_2_A]) +1.0)
		C_position= np.max([phase_0_C,phase_1_C,phase_2_C])/(np.min([phase_0_C,phase_1_C,phase_2_C]) +1.0)
		G_position= np.max([phase_0_G,phase_1_G,phase_2_G])/(np.min([phase_0_G,phase_1_G,phase_2_G]) +1.0)
		T_position= np.max([phase_0_T,phase_1_T,phase_2_T])/(np.min([phase_0_T,phase_1_T,phase_2_T]) +1.0)

		fickett_score += self.look_up_content_probability(A_content,"A")
		fickett_score += self.look_up_content_probability(C_content,"C")
		fickett_score += self.look_up_content_probability(G_content,"G")
		fickett_score += self.look_up_content_probability(T_content,"T")
		
		fickett_score += self.look_up_position_probability(A_position,"A")
		fickett_score += self.look_up_position_probability(C_position,"C")
		fickett_score += self.look_up_position_probability(G_position,"G")
		fickett_score += self.look_up_position_probability(T_position,"T")
			
		return fickett_score

#===================

def mRNA_translate(mRNA):
	return Seq(mRNA).translate()

def protein_param(putative_seqprot):
	return putative_seqprot.isoelectric_point()

def calculate_potential(fasta,strand,output_orf,outfile):
	'''
	Calculate three features: putative peptide length,pI and Fickett
	And assess coding potential based on SVM model
	'''
	strinfoAmbiguous = re.compile("X|B|Z|J|U",re.I)
	ptU = re.compile("U",re.I)
## merged by Yang Ding on 2019-11-23
## 1. all python 2-"file"'s are replaced with python 3-"open"'s
## 2. keep kangyj's check on output_orf
	ftmp_feat = open(outfile + ".feat","w")
	ftmp_svm = open(outfile + ".tmp.1","w")
	ftmp_result = open(outfile,"w")
	if output_orf == 1:
		my_header = ["#ID","transcript_length","peptide_length","Fickett_score","pI","ORF_integrity","ORF_Start","coding_probability","label"]
	else:
		my_header = ["#ID","transcript_length","peptide_length","Fickett_score","pI","ORF_integrity","coding_probability","label"]
	ftmp_result.write("\t".join(map(str,my_header))+"\n")
	fickett_obj = Fickett()
	for seq in seqio.fasta_read(fasta):
		seqid = seq.id
		seqRNA = ptU.sub("T",str(seq.seq).strip())
		'''seqRNA:transcript full sequence'''
		seqRNA = seqRNA.upper()
		seqCDS,start_pos,orf_strand,orf_fullness = FindCDS(seqRNA).longest_orf(strand)
		'''seqCDS:longest ORF'''
		seqprot = mRNA_translate(seqCDS)
		pep_len = len(seqprot) #pep_len = len(seqprot.strip("*"))
		newseqprot = strinfoAmbiguous.sub("",str(seqprot))
		'''exclude ambiguous amio acid X, B, Z, J, Y in peptide sequence'''
		fickett_score = fickett_obj.fickett_value(seqRNA)
		protparam_obj = ProtParam.ProteinAnalysis(str(newseqprot.strip("*")))
		if pep_len > 0:
			#fickett_score = fickett_obj.fickett_value(seqCDS)
			start_pos = start_pos + 1
			isoelectric_point = protein_param(protparam_obj)
		else:
			#fickett_score = 0.0
			start_pos = 0
			orf_fullness = -1
			isoelectric_point = 0.0
		if output_orf == 1:
			output_line = [seqid,len(seqRNA),pep_len,fickett_score,isoelectric_point,orf_fullness,start_pos]
		else:
			output_line = [seqid,len(seqRNA),pep_len,fickett_score,isoelectric_point,orf_fullness]
		ftmp_result.write("\t".join(map(str,output_line))+"\n")
		ftmp_feat.write("\t".join(map(str,[seqid,len(seqRNA),pep_len,fickett_score,isoelectric_point,orf_fullness]))+"\n")
		ftmp_svm.write("".join(map(str,["999"," 1:",pep_len," 2:",fickett_score," 3:",isoelectric_point," 4:",orf_fullness]))+"\n")
	ftmp_result.close()
	ftmp_feat.close()
	ftmp_svm.close()
	#return 0

	'''
	calculate the coding probability using LIBSVM
	'''
	sys.stderr.write("[INFO] Predicting coding potential, please wait ...\n")
	
	'''
		set directories and check depending tools existance
	'''
	script_dir,filename = os.path.split(os.path.abspath(sys.argv[0]))
	data_dir = script_dir + "/../data/"
	lib_dir = script_dir + "/../libs/"
	app_svm_scale = lib_dir + "libsvm/libsvm-3.18/svm-scale"
	app_svm_predict = lib_dir + "libsvm/libsvm-3.18/svm-predict"
	os.system('test -x '+ app_svm_scale + ' || echo \"[ERROR] No excutable svm-scale on CPC2 path!\" > /dev/stderr')
	os.system('test -x '+ app_svm_predict + ' || echo \"[ERROR] No excutable svm-predict on CPC2 path!\" > /dev/stderr')
	
	cmd = app_svm_scale + ' -r ' + data_dir + 'cpc2.range ' + outfile + '.tmp.1 > ' + outfile + '.tmp.2 &&'
	cmd = cmd + app_svm_predict + ' -b 1 -q ' + outfile + '.tmp.2 ' + data_dir + 'cpc2.model ' + outfile + '.tmp.out'
	#cmd = cmd + 'awk -vOFS="\\t" \'{if ($1 == 1){print $2,"coding"} else if ($1 == 0){print $2,"noncoding"}}\' ' + outfile + '.tmp.1 > ' + outfile + '.tmp.2 &&'
	#cmd = cmd + 'paste ' + outfile + '.feat ' + outfile + '.tmp.2 >>' + outfile
	if (sys.version_info[0] == 2):
		(exitstatus, outtext) = commands.getstatusoutput(cmd)
	else:
		(exitstatus, outtext) = subprocess.getstatusoutput(cmd)
	
	'''deal with the output'''
	#print outfile + '.tmp.out'
	tmp_file = open(outfile + '.tmp.out','r')
	prob = {}
	i = 0
	# get libsvm output
	for line in tmp_file:
		i = i + 1
		array = line.split(' ')
		if array[0] == "1":
			label = "coding"
		else:
			label = "noncoding"
		prob[i] = str(array[1]) + '\t' + label
	tmp_file.close()
	# paste to features
	tmp_file = open(outfile,'r')
	out_file = open(outfile + '.txt','w')
	i = 0
	for line in tmp_file:
		i = i + 1
		if i == 1:
			out_file.write(line)
		else:
			line = line.rstrip('\n')
			out_file.write(line)
			out_file.write('\t' + prob[i] + '\n')
	tmp_file.close()
	out_file.close()
	#	subprocess.call("Rscript " + outfile + '.r', shell=True)
	#except:
	#	pass
	if exitstatus == 0:
		os.system('rm -f ' + outfile + '.tmp.1 ' + outfile + '.tmp.2 ' + outfile + '.tmp.out ' + outfile)
		rm_cmd = "rm -f " + outfile + '.feat'
		if (sys.version_info[0] == 2):
			commands.getstatusoutput(rm_cmd)
		else:
			subprocess.getstatusoutput(rm_cmd)
		sys.stderr.write("[INFO] Running Done!\n")
		return 0
	else:
		sys.stderr.write("[ERROR] Prediction error!\n")
		return -1

if __name__ == "__main__":
	sys.exit(__main())
