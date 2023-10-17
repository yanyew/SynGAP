'''
This module deals with input or output in commom formats
'''

import sys
import compress
import os
from Bio import SeqIO

def variant_snpindel_pop(total_fn):
	f = compress.gz_file(total_fn,"r")
	for line in f:
		if line.startswith("#"):continue
		chrom,position1,position2,ref,alt,qual,group_test_pvalue,depth_ref,depth_alt,depth_ref_samples,depth_alt_samples,genotype,other = line.rstrip("\n").split("\t")
		yield [chrom,position1,position2,ref,alt,qual,group_test_pvalue,depth_ref,depth_alt,depth_ref_samples,depth_alt_samples,genotype]
	f.close()


def merge_region(regions):
	# regions must be sorted
	mergedregion = []
	if len(regions) > 0:
		initstart,initend = regions[0]	
	for start,end in regions:
		if start <= initend:
			initend = end
		else:
			mergedregion.append([initstart,initend])
			initstart = start
			initend = end
	mergedregion.append([initstart,initend])
	return mergedregion

def arf_read(arffn):
	f = compress.gz_file(arffn,"r")
	for line in f:
		if line.startswith("#"):continue
		rname,rleng,rstart,rend,rseq,gname,gleng,gstart,gend,gseq,gstrand,nmismatch,mathclabel = line.rstrip("\n").split("\t")
		yield [rname,rleng,rstart,rend,rseq,gname,gleng,gstart,gend,gseq,gstrand,nmismatch,mathclabel]
	f.close()

def fileread(fn):
	if not os.path.isfile(fn):
		sys.stderr.write("[Error] '%s' is not a file\n"%fn)
		exit(1)
	if fn.endswith(".gz"):
		f = compress.gz_file(fn,"r")
	elif fn.endswith(".bz2"):
		f = compress.bz2file(fn)
	return f

def bed6_parse(fn):
	f = compress.gz_file(fn,"r")
	for line in f:
		if line.startswith("#"):continue
		chrom,start,end,name,score,strandother = line.rstrip().split("\t",5)
		yield [chrom,start,end,name,score,strandother]
	f.close()

def gff3_parse(fn):
	f = compress.gz_file(fn,"r")
	for line in f:
		if line.startswith("#"):continue
		chrom,source,seqtype,start,end,score,strand,phase,attributes = line.rstrip("\n").split("\t")
		yield [chrom,source,seqtype,start,end,score,strand,phase,attributes]
	f.close()

def refgene_parse(fn):
	f = compress.gz_file(fn,"r")
	for line in f:
		if line.startswith("#"):continue
		num,nm_name,chrom,strand,exon_s,exon_e,cds_s,cds_e,exon_num,exonstarts,exonends,uniq_id,symbol, kown1, kown2, exon_status = line.rstrip().split("\t")
		yield[num,nm_name,chrom,strand,exon_s,exon_e,cds_s,cds_e,exon_num,exonstarts,exonends,uniq_id,symbol, kown1, kown2, exon_status]
	f.close()

def miRNA_target_parse(fn):
	f = compress.gz_file(fn,"r")
	for line in f:
		if line.startswith("#"):continue
		arr = line.rstrip("\n").split("\t")
		microRNAid,detalmciroRNA,target_Genes = arr[0:3]
		UTR = arr[-3]
		pairing = arr[-2]
		miseq = arr[-1]
		yield[microRNAid,detalmciroRNA,target_Genes,UTR,pairing,miseq]
	f.close()


def sigfile_parse(fn):
	f = compress.gz_file(fn,"r")
	for line in f:
		if line.startswith("#"):continue
		anno1,anno2,fc,rawp,fdr = line.rstrip("\n").split("\t")
		yield [anno1,anno2,fc,rawp,fdr]
	f.close()

def soap_aln_parse(fn):
	f = compress.gz_file(fn,"r")
	for line in f:
		if line.startswith("#"):continue
		seqid,seqread,qual,mcounts,PEtag,length,strand,chrom,sitestart1,mismatch,cigar,match = line.rstrip("\n").split("\t")
		yield [seqid,seqread,qual,mcounts,PEtag,length,strand,chrom,sitestart1,mismatch,cigar,match]
	f.close()

def fasta_read(fn):
	"""
	ID: gi|2765658|emb|Z78533.1|CIZ78533
	Name: gi|2765658|emb|Z78533.1|CIZ78533
	Description: gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA
	Number of features: 0
	Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...GGG', SingleLetterAlphabet())
	"""
	f = compress.gz_file(fn,"r")
	for seq in SeqIO.parse(f,"fasta"):
		yield seq
	f.close()

def fastq_read(fn,qual):
	if qual not in ['fastq-sanger','fastq-solexa','fastq-illumina']:
		sys.stderr.write("[ERROR] Unknown quality\n")
		exit(1)
	
	f = fileread(fn)
	"""
	rec.seq
	rec.letter_annotations['phred_quality']
	"""
	for seq in SeqIO.parse(f,qual):
		yield seq
	f.close()

def bwt_parse(fn):
	f = compress.gz_file(fn,"r")
	for line in f:
		query_id,strand,subject_id,pos,seq,qual,score,mismatch = line.rstrip("\n").split("\t")
		yield [query_id,strand,subject_id,pos,seq,qual,score,mismatch]
	f.close()

def blast6_parse(fn):
	f = compress.gz_file(fn,"r")
	for line in f:
		if line.startswith("#"):continue
		try:
			query_id, subject_id, identity, alignment_length, mismatches, gap_opens, qstart, qend, sstart, send, evalue, bitscore = line.rstrip("\n").split("\t")
		except:
			sys.stderr.write("[WARN] blast can not parse '%s'"%line)
			continue
		yield [query_id, subject_id, identity, alignment_length, mismatches, gap_opens, qstart, qend, sstart, send, evalue, bitscore]
	f.close()

def gtf_parse(fn,add="chr"):
	f = compress.gz_file(fn,"r")
	for line in f:
		if line.startswith("#"):continue
		chrom,rnatype,region_type,start,end,score,strand,codon,commnet = line.rstrip("\n").split("\t")
		yield[add+chrom.lstrip("chr"),rnatype,region_type,start,end,score,strand,codon,commnet]
	f.close()

## def blast_parse
## def sam or bam parse ...
##

if __name__ == "__main__":
	a = [[1,10],[17,22],[40,44],[42,47],[46,100],[101,408]]
	print (a)
	print (merge_region(a))



