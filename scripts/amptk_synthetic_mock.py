#!/usr/bin/env python

from __future__ import print_function
from __future__ import division
from builtins import str
from builtins import range
from past.utils import old_div
import sys
import random
import math
import re
from natsort import natsorted
from Bio.SeqUtils import GC

def merge_dicts(*dict_args):
    '''
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    '''
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result

def randomSeq(length, GC):
    global RandomSeq
    pctGC = length * GC / 2
    pctGC = math.trunc(pctGC)
    pctAT = old_div(length,2) - pctGC
    Seq = "A"*pctAT + "T"*pctAT + "G"*pctGC + "C"*pctGC
    SeqList = list(Seq)
    random.shuffle(SeqList)
    RandomSeq = "".join(SeqList)
    return RandomSeq

def pashuffle(string, perc=10):
    data = list(string)
    for index, letter in enumerate(data):
        if random.randrange(0, 100) < old_div(perc,2):
            new_index = random.randrange(0, len(data))
            data[index], data[new_index] = data[new_index], data[index]
    return "".join(data)

def main():
	SSU = "CTTGGTCATTTAGAGGAAGTAAAAGTCGTAACAAGGTTTCCGTAGGTGAACCTGCGG"
	FES = "CAACAACGGATCTCTTGGTTCTCGCATCGATGAAGAACGCAGCNNNNNNNNNNNNNNNNNNGTGAATCATCGAATCTTTGAACGCACATTGCGCTCCTTGGTAT"
	LSU = "TGACCTCAAATCAGGTAGGAGTACCCGCTGAACTTAAGCATATCAATAAGCGGAGGA"

	'''
	ITS3 = GCATCGATGAAGAACGCAGC
	ITS3_KYO1 = 
	fITS9 = GAACGCAGCRAAIIGYGA
	probe   GAACGCAGCGAAATGCGA
	58A1F	GCATCGATGAAGAACGC
	58A2F	ATCGATGAAGAACGCAG
	58A2R	CTGCGTTCTTCATCGAT
	'''
	#ITS1: 57 bp from SSU + 43 bp for 5.8S = 100
	#5.8S: 158 bp
	#ITS2: 57 bp from LSU + 43 bp for 5.8S = 100

	#qPCR probe 18 bp random assignment, store in dictionary for re.sub later
	probes = {}
	for i in range(1,20):
		randomSeq(18, 0.65)
		probes[i] = RandomSeq

	#create the ~ 50% GC size fragments
	count = 1
	seqs = []
	for i in 100,150,200,250,300,350:
		randomSeq(i, 0.55)
		Tup = [count,i,RandomSeq]
		seqs.append(Tup)
		count += 1    

	Seqs = {}
	forward = 0
	reverse = 5
	for i in range(1,7):
		its1GC = math.trunc(GC(seqs[forward][2]))
		its2GC = math.trunc(GC(seqs[reverse][2]))
		its1Len = len(seqs[forward][2]) + 100
		its2Len = len(seqs[reverse][2]) + 100
		sequence = SSU + seqs[forward][2] + FES + seqs[reverse][2] + LSU
		qPCR = probes.get(i)
		sequence = re.sub('NNNNNNNNNNNNNNNNNN', qPCR, sequence)
		name = "SynMock_" + str(i) + ";qPCR_probe=" + qPCR + ";Len_ITS1=" + str(its1Len) + ";GC_ITS1=" + str(its1GC) + ";Len_ITS2=" + str(its1Len) + ";GC_ITS2=" + str(its2GC)
		Seqs[name] = sequence
		forward += 1
		reverse -= 1
 
	#now vary GC content
	count = 1
	GeeCee = []
	for p in 0.40, 0.47, 0.63, 0.70:
		randomSeq(190, p)
		Tup = [count,i,RandomSeq]
		GeeCee.append(Tup)
		count += 1

	GCSeqs = {}
	forward = 0
	reverse = 3
	for i in 7,8,9,10:
		its1GC = math.trunc(GC(GeeCee[forward][2]))
		its2GC = math.trunc(GC(GeeCee[reverse][2]))
		its1Len = len(GeeCee[forward][2]) + 100
		its2Len = len(GeeCee[reverse][2]) + 100
		sequence = SSU + GeeCee[forward][2] + FES + GeeCee[reverse][2] + LSU
		qPCR = probes.get(i)
		sequence = re.sub('NNNNNNNNNNNNNNNNNN', qPCR, sequence)
		name = "SynMock_" + str(i) + ";qPCR_probe=" + qPCR + ";Len_ITS1=" + str(its1Len) + ";GC_ITS1=" + str(its1GC) + ";Len_ITS2=" + str(its1Len) + ";GC_ITS2=" + str(its2GC)
		GCSeqs[name] = sequence
		forward += 1
		reverse -= 1

	#lets get some closely related paralogs, partial shuffle
	Paralog = {}
	paralog1 = pashuffle(seqs[2][2])
	paralog2 = pashuffle(seqs[3][2])
	sequence1 = SSU + paralog1 + FES + paralog2 + LSU
	sequence2 = SSU + paralog2 + FES + paralog1 + LSU
	para1GC = math.trunc(GC(sequence1))
	para2GC = math.trunc(GC(sequence2))
	para1Len = len(paralog1) + 100
	para2Len = len(paralog2) + 100
	qPCR1 = probes.get(11)
	qPCR2 = probes.get(12)
	sequence1 = re.sub('NNNNNNNNNNNNNNNNNN', qPCR1, sequence1)
	sequence2 = re.sub('NNNNNNNNNNNNNNNNNN', qPCR2, sequence2)
	name1 = "SynMock_11" + ";qPCR_probe=" + qPCR1 + ";Len_ITS1=" + str(para1Len) + ";GC_ITS1=" + str(para1GC) + ";Len_ITS2=" + str(para2Len) + ";GC_ITS2=" + str(para2GC)
	name2 = "SynMock_12" + ";qPCR_probe=" + qPCR2 + ";Len_ITS1=" + str(para2Len) + ";GC_ITS1=" + str(para2GC) + ";Len_ITS2=" + str(para1Len) + ";GC_ITS2=" + str(para1GC)
	Paralog[name1] = sequence1
	Paralog[name2] = sequence2
	

	#merge dictionaries
	final = merge_dicts(Seqs,GCSeqs,Paralog)

	for key,value in natsorted(iter(final.items())):
		print(">%s\n%s" % (key, value))

if __name__ == "__main__":
	main()