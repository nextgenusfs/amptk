#!/usr/bin/env python
#script written by Robert Edgar http://drive5.com
import sys
import os
import uc
import die
import fasta

if len(sys.argv)<2:
    print "Usage: " + sys.argv[0] + " table.uc > table.otu.txt"
    os._exit(1)

FileName = sys.argv[1]
OutName = sys.argv[2]

def GetSampleId(Label):
	Fields = Label.split(";")
	for Field in Fields:
		if Field.startswith("barcodelabel="):
			return Field[13:]
	die.Die("barcodelabel= not found in read label '%s'" % Label)

def OnRec():
	global OTUs, Samples, OTUTable
	if uc.Type != 'H':
		return

	OTUId = uc.TargetLabel
	if OTUId not in OTUIds:
		OTUIds.append(OTUId)
		OTUTable[OTUId] = {}

	SampleId = GetSampleId(uc.QueryLabel)
	if SampleId not in SampleIds:
		SampleIds.append(SampleId)

	N = fasta.GetSizeFromLabel(uc.QueryLabel, 1)
	try:
		OTUTable[OTUId][SampleId] += N
	except:
		OTUTable[OTUId][SampleId] = N

OTUIds = []
SampleIds = []
OTUTable = {}

uc.ReadRecs(FileName, OnRec)

handle = open(OutName, 'w')

s = "OTUId"
for SampleId in SampleIds:
	s += "\t" + SampleId
s = s + '\n'
handle.write(s)

for OTUId in OTUIds:
	s = OTUId
	for SampleId in SampleIds:
		try:
			n = OTUTable[OTUId][SampleId]
		except:
			n = 0
		s += "\t" + str(n)
	s = s + '\n'
	handle.write(s)
handle.close()