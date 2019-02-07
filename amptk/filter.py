#!/usr/bin/env python

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import sys
import os
import argparse
import inspect
import subprocess
import sys
import math
from Bio import SeqIO
from natsort import natsorted
import pandas as pd
import numpy as np
from amptk import amptklib

class colr(object):
    GRN = '\033[92m'
    END = '\033[0m'
    WARN = '\033[93m'

class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=50)      

def main(args):
	parser=argparse.ArgumentParser(prog='amptk-filter.py',
		description='''Script inspects output of amptk-OTU_cluster.py and 
		determines useful threshold for OTU output based on a spike-in 
		mock community.''',
		epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
		formatter_class=MyFormatter)

	parser.add_argument('-i','--otu_table', required=True, help='Input OTU table')
	parser.add_argument('-f','--fasta', required=True, help='Input OTUs (multi-fasta)')
	parser.add_argument('-b','--mock_barcode', help='Barocde of Mock community')
	parser.add_argument('-p','--index_bleed',  help='Index Bleed filter. Default: auto')
	parser.add_argument('-t','--threshold', default='max', choices=['sum','max','top25','top10','top5'],help='Threshold to use when calculating index-bleed')
	parser.add_argument('-c','--calculate', default='all', choices=['all', 'in'], help='Calculate index-bleed, if synthetic mock use all otherwise use in')
	parser.add_argument('-s','--subtract', default=0, help='Threshold to subtract')
	parser.add_argument('-n','--normalize', default='y', choices=['y','n'], help='Normalize OTU table prior to filtering')
	parser.add_argument('-m','--mc', help='Multi-FASTA mock community')
	parser.add_argument('-d','--drop', nargs='+', help='samples to drop from table after index-bleed filtering')
	parser.add_argument('--ignore', nargs='+', help='Ignore OTUs during index-bleed')
	parser.add_argument('--delimiter', default='tsv', choices=['csv','tsv'], help='Delimiter')
	parser.add_argument('--col_order', nargs='+', dest="col_order", help='Provide space separated list')
	parser.add_argument('--keep_mock', action='store_true', help='Keep mock sample in OTU table (Default: False)')
	parser.add_argument('--show_stats', action='store_true', help='Show stats datatable STDOUT')
	parser.add_argument('--negatives', nargs='+', help='Negative Control Sample names')
	parser.add_argument('-o','--out', help='Base output name')
	parser.add_argument('--min_reads_otu', default=2, type=int, help='Minimum number of reads per OTU for experiment')
	parser.add_argument('--min_samples_otu', default=1, type=int, help='Minimum number of samples per OTU for experiment')
	parser.add_argument('-u','--usearch', dest="usearch", default='usearch9', help='USEARCH9 EXE')
	parser.add_argument('--debug', action='store_true', help='Remove Intermediate Files')
	args=parser.parse_args(args)
	
	parentdir = os.path.join(os.path.dirname(amptklib.__file__))

	if not args.out:
		#get base name of files
		base = args.otu_table.split(".otu_table")[0]
	else:
		base = args.out

	#remove logfile if exists
	log_name = base + '.amptk-filter.log'
	amptklib.removefile(log_name)

	amptklib.setupLogging(log_name)
	FNULL = open(os.devnull, 'w')
	cmd_args = " ".join(sys.argv)+'\n'
	amptklib.log.debug(cmd_args)
	print("-------------------------------------------------------")

	#initialize script, log system info and usearch version
	amptklib.SystemInfo()
	#Do a version check
	usearch = args.usearch
	amptklib.versionDependencyChecks(usearch)

	#check if otu_table is empty
	amptklib.log.info("Loading OTU table: %s" % args.otu_table)
	check = os.stat(args.otu_table).st_size
	if check == 0:
		amptklib.log.error("Input OTU table is empty")
		sys.exit(1)
	#get the OTU header info (depending on how OTU table was constructed, this might be different, so find it as you need for indexing)
	with open(args.otu_table, 'r') as f:
		first_line = f.readline()
		OTUhead = first_line.split('\t')[0]

	if args.delimiter == 'csv':
		delim = str(',')
		ending = '.csv'
	elif args.delimiter == 'tsv':
		delim = str('\t')
		ending = '.txt'

	#setup outputs
	sorted_table = base+'.sorted'+ending
	normal_table_pct = base+'.normalized.pct'+ending
	normal_table_nums = base+'.normalized.num'+ending
	subtract_table = base+'.normalized.subtract'+ending
	filtered_table = base+'.normalized'+ending
	final_table = base+'.final'+ending
	final_binary_table = base+'.final.binary'+ending
	stats_table = base+'.stats'+ending

	#load OTU table into pandas DataFrame
	df = pd.read_csv(args.otu_table, sep='\t')
	df.set_index(OTUhead, inplace=True)
	headers = df.columns.values.tolist()
	if headers[-1] == 'taxonomy' or headers[-1] == 'Taxonomy':
		otuDict = df[headers[-1]].to_dict()
		del df[headers[-1]]
	else:
		otuDict = False

	#parse OTU table to get count data for each OTU
	AddCounts = {}
	OTUcounts = df.sum(1)
	for x in OTUcounts.index:
		AddCounts[x] = int(OTUcounts[x])

	#now add counts to fasta header
	FastaCounts = base+'.otus.counts.fa'
	OTU_tax = {}
	with open(FastaCounts, 'w') as outfile:
		with open(args.fasta, 'r') as infile:
			for rec in SeqIO.parse(infile, 'fasta'):
				if ';' in rec.id: #this should mean there is taxonomy, so split it
					ID = rec.id.split(';',1)[0]
					tax = rec.id.split(';',1)[-1]
					OTU_tax[ID] = tax
					if ID in AddCounts:
						count = AddCounts.get(ID)
					else:
						count = 0
					outfile.write('>%s;size=%i\n%s\n' % (ID, count, rec.seq))
				else: #no tax, just process
					if rec.id in AddCounts:
						count = AddCounts.get(rec.id)
					else:
						count = 0
					outfile.write('>%s;size=%i\n%s\n' % (rec.id, count, rec.seq))

	amptklib.log.info('OTU table contains {:,} samples, {:,} OTUs, and {:,} reads counts'.format(len(df.columns.values.tolist()), len(df.index), int(df.values.sum())))

	#setup output files/variables
	mock_out = base + '.mockmap.txt'

	if args.mock_barcode: #if user passes a column name for mock
		#check if mock barcode is valid
		validBCs = df.columns.values.tolist()
		if not args.mock_barcode in validBCs:
			amptklib.log.error("%s not a valid barcode." % args.mock_barcode)
			amptklib.log.error("Valid barcodes: %s" % (' '.join(validBCs)))
			sys.exit(1)
		if args.col_order and not args.mock_barcode in args.col_order:
			amptklib.log.error("Error: %s not listed in --col_order." % args.mock_barcode)
			sys.exit(1)
		#make sure there is a --mc passed here otherwise throw error
		if not args.mc:
			amptklib.log.error("If using the -b,--barcode option you must specify a fasta file of mock community via the --mc option")
			sys.exit(1)
		#get default mock community value
		if args.mc == "mock3":
			mock = os.path.join(parentdir, 'DB', 'amptk_mock3.fa')
		elif args.mc == "mock2":
			mock = os.path.join(parentdir, 'DB', 'amptk_mock2.fa')
		elif args.mc == "mock1":
			mock = os.path.join(parentdir, 'DB', 'amptk_mock1.fa')
		elif args.mc == "synmock":
			mock = os.path.join(parentdir, 'DB', 'amptk_synmock.fa')
		else:
			mock = os.path.abspath(args.mc)

		#open mock community fasta and count records
		mock_ref_count = amptklib.countfasta(mock)
	
		#load OTU lengths into dictionary
		SeqLength = amptklib.fastalen2dict(args.fasta)
	
		#map OTUs to mock community, this is fast enough, but running twice, first to get only top hit, then
		amptklib.log.info("Mapping OTUs to Mock Community (USEARCH)")
		cmd = [usearch, '-usearch_global', mock, '-strand', 'plus', '-id', '0.65', '-db', FastaCounts, '-userout', mock_out, '-userfields', 'query+target+id+ql+tl+alnlen+caln+mism+diffs', '-maxaccepts', '0', '-maxrejects', '0']
		amptklib.runSubprocess(cmd, amptklib.log)

		#generate dictionary for name change
		'''
		If args.calculate is set to all, that means the script is trying to measure a synthetic
		mock of some kind.  if that is the case, then chimeras are < 95% identical to mock members
		and variants would be hits in between, i.e 95% > but not the best hit.
		'''
		Results = {}
		errorrate = {}
		with open(mock_out, 'r') as map:
			for line in map:
				line = line.replace('\n', '')
				cols = line.split('\t')
				MockID = cols[0]
				hit = cols[1].split(';size=')
				otuID = hit[0]
				abundance = int(hit[1])
				pident = float(cols[2])
				length = int(cols[4])
				mism = int(cols[7])
				diffs = int(cols[8])
				score = abundance * pident * length
				if not otuID in errorrate:
					errorrate[otuID] = [MockID,diffs]
				else:
					olderror = errorrate.get(otuID)
					if diffs < olderror[1]:
						errorrate[otuID] = [MockID,diffs]
				if not MockID in Results:
					Results[MockID] = [(otuID,abundance,pident,length,mism,diffs,score)]
				else:
					Results[MockID].append((otuID,abundance,pident,length,mism,diffs,score))
	
		found_dict = {}
		chimeras = []
		variants = []
		missing = []
		for k,v in natsorted(list(Results.items())):
			besthit = []
			#v is a list of tuples of results, parse through to get best hit
			for y in v:
				if y[2] >= 97.0:
					besthit.append(y)
				elif y[2] >= 95.0 and y[2] < 97.0:
					if not y[0] in variants:
						variants.append(y[0])
				else:
					if not y[0] in chimeras:
						chimeras.append(y[0])
			if len(besthit) > 0:
				besthit.sort(key=lambda x: x[1], reverse=True)
				best = sorted(besthit[:3], key=lambda x: x[6], reverse=True)
				found_dict[k] = best[0]
			else:
				missing.append(k)
			
		#make name change dict
		annotate_dict = {}
		seen = []
		for k,v in natsorted(list(found_dict.items())):
			ID = v[0].replace('_chimera', '')
			newID = k+'_pident='+str(v[2])+'_'+v[0]
			annotate_dict[ID] = newID
			if not v[0] in seen:
				seen.append(v[0])
		if args.calculate == 'all':
			chimeras = [x for x in chimeras if x not in seen]
			variants = [x for x in variants if x not in seen]
			for i in chimeras:
				annotate_dict[i] = i+'_suspect_mock_chimera'
			for x in variants:
				annotate_dict[x] = x+'_suspect_mock_variant'
		if len(missing) > 0:
			amptklib.log.info("%i mock missing: %s" % (len(missing), ', '.join(missing)))
	else:
		otu_new = args.fasta

	#rename OTUs
	if args.mock_barcode:
		df.rename(index=annotate_dict, inplace=True)

	#sort the table
	df2 = df.reindex(index=natsorted(df.index))
	if not args.col_order:
		amptklib.log.info("Sorting OTU table naturally")
		df = df2.reindex(columns=natsorted(df2.columns))
	else:
		amptklib.log.info("Sorting OTU table by user defined order (--col_order)")
		col_headers = args.col_order
		#check if all names in headers or not
		for i in col_headers:
			if not i in df2.columns.values:
				col_headers.remove(i)
		df = df2.reindex(columns=col_headers)
	SortedTable = df
	if otuDict:
		df['Taxonomy'] = pd.Series(otuDict)
		df.to_csv(sorted_table, sep=delim)
		del df['Taxonomy']
	else:
		df.to_csv(sorted_table, sep=delim)

	#get sums of columns
	fs = df.sum(axis=0)
	#fs.to_csv('reads.per.sample.csv')
	otus_per_sample_original = df[df > 0].count(axis=0, numeric_only=True)
	filtered = pd.DataFrame(df, columns=fs.index)
	filt2 = filtered.loc[(filtered != 0).any(1)]
	tos = filt2.sum(axis=1)
	fotus = tos[tos >= args.min_reads_otu] #valid allele must be found atleast from than 2 times, i.e. no singletons
	if len(fotus.index) < len(tos.index):
		diff = len(tos.index) - len(fotus.index)
		amptklib.log.info("Removing {:,} OTUs according to --min_reads_otu {:,}".format(diff, args.min_reads_otu))
	filt3 = pd.DataFrame(filt2, index=fotus.index)

	if args.normalize == 'y':
		#normalize the OTU table
		normal = filt3.truediv(fs)
		if otuDict:
			normal['Taxonomy'] = pd.Series(otuDict)
			normal.to_csv(normal_table_pct, sep=delim)
			del normal['Taxonomy']
		else:
			normal.to_csv(normal_table_pct, sep=delim)
		#normalize back to read counts, pretend 100,000 reads in each
		norm_round = np.round(normal.multiply(100000), decimals=0)
		if otuDict:
			norm_round['Taxonomy'] = pd.Series(otuDict)
			norm_round.to_csv(normal_table_nums, sep=delim)
			del norm_round['Taxonomy']
		else:
			norm_round.to_csv(normal_table_nums, sep=delim)
		amptklib.log.info("Normalizing OTU table to number of reads per sample")
	else:
		norm_round = filt3

	if args.mock_barcode:
		#now calculate the index-bleed in both directions (into the mock and mock into the other samples)
		mock = []
		sample = []
		#get names from mapping
		for k,v in list(annotate_dict.items()):
			if not '_suspect_mock_' in v:
				mock.append(v)
		for i in norm_round.index:
			if not i in mock:
				sample.append(i)
		if args.ignore:
			mock = [x for x in mock if x not in args.ignore]
			sample = [x for x in sample if x not in args.ignore]
		#first calculate bleed out of mock community
		#slice normalized dataframe to get only mock OTUs from table
		mock_df = pd.DataFrame(norm_round, index=mock)
		#if there are samples to drop, make sure they aren't being used in this calculation
		if args.drop:
			mock_df.drop(args.drop, axis=1, inplace=True)
		#get total number of reads from mock OTUs from entire table
		total = np.sum(np.sum(mock_df,axis=None))
		#now drop the mock barcode sample
		mock_df.drop(args.mock_barcode, axis=1, inplace=True)
		#get number of reads that are result of bleed over
		bleed1 = np.sum(np.sum(mock_df,axis=None))
		#calculate rate of bleed by taking num reads bleed divided by the total
		bleed1max = bleed1 / float(total)
	
		#second, calculate bleed into mock community
		#get list of mock OTUs not found in any other sample -> these are likely chimeras
		mock_only = pd.DataFrame(norm_round, index=list(norm_round.index), columns=[args.mock_barcode])
		mock_OTUs_zeros = mock_only.loc[(mock_only==0).all(axis=1)]
		theRest = [x for x in list(norm_round.columns.values) if x not in [args.mock_barcode]]
		non_mocks = pd.DataFrame(norm_round, index=sample, columns=theRest)
		non_mock_zeros = non_mocks.loc[(non_mocks==0).all(axis=1)]
		zeros = [x for x in list(non_mock_zeros.index) if x not in list(mock_OTUs_zeros.index)]
		if len(zeros) > 0:
			amptklib.log.info("Found {:,} mock chimeras (only in mock sample and not mapped to mock sequences) excluding from index-bleed calculation".format(len(zeros)))
			amptklib.log.debug('{:}'.format(', '.join(zeros)))
		#now get updated list of samples, dropping chimeras
		samples_trimmed = [x for x in sample if x not in zeros]
		#slice the OTU table to get all OTUs that are not in mock community from the mock sample
		sample_df = pd.DataFrame(norm_round, index=samples_trimmed, columns=[args.mock_barcode])
		#get total number of reads that don't belong in mock
		bleed2 = np.sum(np.sum(sample_df,axis=None))
		#now pull the entire mock sample
		mock_sample = pd.DataFrame(norm_round, columns=[args.mock_barcode])
		#calcuate bleed into mock by taking num reads that don't belong divided by the total, so this is x% of bad reads in the mock
		bleed2max = bleed2 / float(np.sum(mock_sample.sum(axis=1)))
		#autocalculate the subtraction filter by taking the maximum value that doesn't belong
		subtract_num = max(sample_df.max())

		#get max values for bleed
		#can only use into samples measurement if not using synmock
		if args.calculate == 'all':
			if bleed1max > bleed2max:
				bleedfilter = math.ceil(bleed1max*1000)/1000
			else:
				bleedfilter = math.ceil(bleed2max*1000)/1000
			amptklib.log.info("Index bleed, mock into samples: %f%%.  Index bleed, samples into mock: %f%%." % (bleed1max*100, bleed2max*100))
		else:
			bleedfilter = math.ceil(bleed2max*1000)/1000
			amptklib.log.info("Index bleed, samples into mock: %f%%." % (bleed2max*100))
		
	else:
		bleedfilter = args.index_bleed #this is value needed to filter MiSeq, Ion is likely less, but shouldn't effect the data very much either way.

	if args.index_bleed:
		args.index_bleed = float(args.index_bleed)
		amptklib.log.info("Overwriting auto detect index-bleed, setting to %f%%" % (args.index_bleed*100))
		bleedfilter = args.index_bleed
	else:
		if bleedfilter:
			amptklib.log.info("Will use value of %f%% for index-bleed OTU filtering." % (bleedfilter*100))
		else:
			bleedfilter = 0 #no filtering if you don't pass -p or -b 
			amptklib.log.info("No spike-in mock (-b) or index-bleed (-p) specified, thus not running index-bleed filtering") 

	if bleedfilter > 0.05:
		amptklib.log.info("Index bleed into samples is abnormally high (%f%%), if you have biological mock you should use `--calculate in`" % (bleedfilter*100))


	#to combat barcode switching, loop through each OTU filtering out if less than bleedfilter threshold
	cleaned = []
	for row in norm_round.itertuples():
		result = [row[0]]
		if args.threshold == 'max':
			total = max(row[1:]) #get max OTU count from table to calculate index bleed from.
		elif args.threshold == 'sum':
			total = sum(row[1:])
		elif args.threshold == 'top25':
			top = sorted(row[1:], key=int, reverse=True)
			topn = int(round(len(row[1:])*0.25))
			total = sum(top[:topn])
		elif args.threshold == 'top10':
			top = sorted(row[1:], key=int, reverse=True)
			topn = int(round(len(row[1:])*0.10))
			total = sum(top[:topn])
		elif args.threshold == 'top5':
			top = sorted(row[1:], key=int, reverse=True)
			topn = int(round(len(row[1:])*0.05))
			total = sum(top[:topn])
		sub = total * bleedfilter
		for i in row[1:]:
			if i < sub:
				i = 0
			result.append(i)
		cleaned.append(result)

	header = [OTUhead]
	for i in norm_round.columns:
		header.append(i)

	#create data frame of index bleed filtered results       
	final = pd.DataFrame(cleaned, columns=header)
	final.set_index(OTUhead, inplace=True)
	
	if args.drop: #if user has passed samples to drop, do it here, subtract drop list from Header
		amptklib.log.info("Dropping %i samples from table: %s" % (len(args.drop), ', '.join(args.drop)))
	
		colsdrop = []
		for x in args.drop:
			if x in header:
				colsdrop.append(x)
		#now drop those columns
		final.drop(colsdrop, axis=1, inplace=True)

	if args.subtract != 'auto':
		subtract_num = int(args.subtract)
	else:
		try:
			subtract_num = int(subtract_num)
			amptklib.log.info("Auto subtract filter set to %i" % subtract_num)
		except NameError:
			subtract_num = 0
			amptklib.log.info("Error: to use 'auto' subtract feature, provide a sample name to -b,--mock_barcode.")
	if subtract_num != 0:
		amptklib.log.info("Subtracting %i from OTU table" % subtract_num)
		sub = final.subtract(subtract_num)
		sub[sub < 0] = 0 #if negative, change to zero
		sub = sub.loc[~(sub==0).all(axis=1)]
		sub = sub.astype(int)
		if otuDict:
			sub['Taxonomy'] = pd.Series(otuDict)
			sub.to_csv(subtract_table, sep=delim)
			del sub['Taxonomy']
		else:
			sub.to_csv(subtract_table, sep=delim)
		otus_if_sub = sub[sub > 0].count(axis=0, numeric_only=True)
		final = sub.astype(int)
	otus_per_sample = final[final > 0].count(axis=0, numeric_only=True)
	stats = pd.concat([fs, otus_per_sample_original, otus_per_sample], axis=1)
	stats.columns = ['reads per sample', 'original OTUs', 'final OTUs']
	stats.fillna(0, inplace=True)
	stats = stats.astype(int)
	if args.show_stats:
		print(stats.to_string())
	stats.to_csv(stats_table, sep=delim)
	#after all filtering, get list of OTUs in mock barcode
	if args.mock_barcode:
		mocks = final[args.mock_barcode]
		mocks = mocks.loc[~(mocks==0)].astype(int)
		totalmismatches = 0
		totallength = 0
		chimera_count = 0
		variant_count = 0
		for otu in mocks.index:
			count = mocks[otu]
			if 'suspect_mock' in otu:
				if 'chimera' in otu:
					chimera_count += 1
				if 'variant' in otu:
					variant_count += 1
				otu = otu.split('_',1)[0]
			else:
				otu = otu.split('_',-1)[-1]
			otu_length = SeqLength.get(otu)
			countlen = otu_length * count
			totallength += countlen
			if otu in errorrate:
				otu_diffs = errorrate.get(otu)[1]
				totaldiffs = otu_diffs * count
				totalmismatches += totaldiffs
			else:
				totalmismatches += countlen
		e_rate = totalmismatches / float(totallength) * 100
		amptklib.log.info(args.mock_barcode + ' sample has '+'{0:,}'.format(len(mocks))+' OTUS out of '+'{0:,}'.format(mock_ref_count)+ ' expected; '+'{0:,}'.format(variant_count)+ ' mock variants; '+ '{0:,}'.format(chimera_count)+ ' mock chimeras; Error rate: '+'{0:.3f}%'.format(e_rate))

	if not args.keep_mock:
		try:
			final.drop(args.mock_barcode, axis=1, inplace=True)
		except:
			pass
		
	#drop OTUs that are now zeros through whole table
	final = final.loc[~(final==0).all(axis=1)]
	final = final.astype(int)

	#output filtered normalized table
	if otuDict:
		final['Taxonomy'] = pd.Series(otuDict)
		final.to_csv(filtered_table, sep=delim)
		del final['Taxonomy']
	else:
		final.to_csv(filtered_table, sep=delim)

	#convert to binary
	final[final > 0] = 1

	#apply min_sample_otu here (most stringent filter, not sure I would use this unless you know what you are doing)
	los = final.sum(axis=1)
	fotus = los[los >= args.min_samples_otu]
	keep = fotus.index
	final2 = pd.DataFrame(final, index=keep)
	diff = len(final.index) - len(keep)
	if diff > 0:
		amptklib.log.info('Dropped {:,} OTUs found in fewer than {:,} samples'.format(diff, args.min_samples_otu))

	#drop samples that don't have any OTUs after filtering
	final3 = final2.loc[:, (final2 != 0).any(axis=0)]
	final3 = final3.astype(int)

	#get the actual read counts from binary table
	merge = {}
	for index, row in final3.items():
		merge[index] = []
		for i in range(0, len(row)):
			if row[i] == 0:
				merge[index].append(row[i])
			else:
				merge[index].append(SortedTable[index][row.index[i]])

	FiltTable = pd.DataFrame(merge, index=list(final3.index))
	FiltTable.index.name = '#OTU ID'

	#order the filtered table
	#sort the table
	FiltTable2 = FiltTable.reindex(index=natsorted(FiltTable.index))
	if not args.col_order:
		FiltTable = FiltTable2.reindex(columns=natsorted(FiltTable2.columns))
	else:
		col_headers = args.col_order
		#check if all names in headers or not
		for i in col_headers:
			if not i in FiltTable2.columns.values:
				col_headers.remove(i)
		FiltTable = FiltTable2.reindex(columns=col_headers)

	#check for negative samples and how many OTUs are in these samples
	#if found, filter the OTUs and alert user to rebuild OTU table, I could do this automatically, but would then require
	#there to be reads passed to this script which seems stupid.  Just deleting the OTUs is probably not okay....
	if args.negatives:
		if len(args.negatives) > 1: #if greater than 1 then assuming list of sample names
			Neg = args.negatives
		else:
			if os.path.isfile(args.negatives[0]): #check if it is a file or not
				Neg = []
				with open(args.negatives[0], 'r') as negfile:
					for line in negfile:
						line = line.replace('\n', '')
						Neg.append(line)
			else:
				Neg = args.negatives
		#Now slice the final OTU table, check if values are valid
		NotFound = []
		for i in Neg:
			if not i in FiltTable.columns.values:
				Neg.remove(i)
				NotFound.append(i)
		if len(NotFound) > 0:
			amptklib.log.info('Samples not found: %s' % ' '.join(NotFound))
		#slice table
		NegTable = FiltTable.reindex(columns=Neg)
		#drop those that are zeros through all samples, just pull out OTUs found in the negative samples
		NegTable = NegTable.loc[~(NegTable==0).all(axis=1)]
		NegOTUs = list(NegTable.index)
		#now make sure you aren't dropping mock OTUs as you want to keep those for filtering new OTU table
		NegOTUs = [item for item in NegOTUs if item not in mock]
	else:
		NegOTUs = []

	#check if negative OTUs exist, if so, then output updated OTUs and instructions on creating new OTU table
	if len(NegOTUs) > 0:
		amptklib.log.info("%i OTUs are potentially contamination" % len(NegOTUs))
		otu_clean = base + '.cleaned.otus.fa'
		with open(otu_clean, 'w') as otu_update:
			with open(args.fasta, "rU") as myfasta:
				for rec in SeqIO.parse(myfasta, 'fasta'):
					if not rec.id in NegOTUs:
						SeqIO.write(rec, otu_update, 'fasta')
		amptklib.log.info("Cleaned OTUs saved to: %s" % otu_clean)
		amptklib.log.info("Generate a new OTU table like so:\namptk remove -i %s --format fasta -l %s -o %s\nvsearch --usearch_global %s --db %s --strand plus --id 0.97 --otutabout newOTU.table.txt\n" % (base+'.demux.fq', ' '.join(Neg), base+'.cleaned.fa', base+'.cleaned.fa', otu_clean))

	else: #proceed with rest of script    
		#output final table
		if otuDict:
			FiltTable['Taxonomy'] = pd.Series(otuDict)
			FiltTable.to_csv(final_table, sep=delim)
			del FiltTable['Taxonomy']
		else:
			FiltTable.to_csv(final_table, sep=delim)
		finalSamples = FiltTable.columns.values.tolist()
		if 'Taxonomy' in finalSamples:
			numFinalSamples = len(finalSamples) - 1
		else:
			numFinalSamples = len(finalSamples)
		amptklib.log.info('Filtered OTU table contains {:,} samples, {:,} OTUs, and {:,} read counts'.format(numFinalSamples, len(FiltTable.index), FiltTable.values.sum()))
		if numFinalSamples < len(df.columns.values.tolist()):
			diffSamples = [item for item in headers if item not in FiltTable.columns.values.tolist()]
			amptklib.log.info('Samples dropped: %s' % (','.join(diffSamples)))
		#output binary table
		if otuDict:
			final3['Taxonomy'] = pd.Series(otuDict)
			final3.to_csv(final_binary_table, sep=delim)
		else:
			final3.to_csv(final_binary_table, sep=delim)

		#generate final OTU list for taxonomy
		amptklib.log.info("Finding valid OTUs")
		otu_new = base + '.filtered.otus.fa'
		with open(otu_new, 'w') as otu_update:
			with open(args.fasta, "rU") as myfasta:
				for rec in SeqIO.parse(myfasta, 'fasta'):
					if ';' in rec.id:
						rec.id = rec.id.split(';',1)[0]
					if args.mock_barcode:
						#map new names of mock
						if rec.id in annotate_dict:
							newname = annotate_dict.get(rec.id)
							rec.id = newname
							rec.description = ''
					if rec.id in final3.index:
						if rec.id in OTU_tax:
							otu_update.write('>%s;%s\n%s\n' % (rec.id, OTU_tax.get(rec.id), rec.seq))
						else:
							otu_update.write('>%s\n%s\n' % (rec.id, rec.seq))

				
		#tell user what output files are
		print("-------------------------------------------------------")
		print("OTU Table filtering finished")
		print("-------------------------------------------------------")
		print("OTU Table Stats:      %s" % stats_table)
		print("Sorted OTU table:     %s" % sorted_table)
		if not args.debug:
			for i in [normal_table_pct, normal_table_nums, subtract_table, mock_out, FastaCounts]:
				amptklib.removefile(i)
		else:   
			print("Normalized (pct):     %s" % normal_table_pct)
			print("Normalized (10k):     %s" % normal_table_nums)
			if args.subtract != 0:
				print("Subtracted table:     %s" % subtract_table)
		print("Normalized/filter:    %s" % filtered_table)
		print("Final Binary table:   %s" % final_binary_table)
		print("Final OTU table:      %s" % final_table)
		print("Filtered OTUs:        %s" % otu_new)
		print("-------------------------------------------------------")

		if 'darwin' in sys.platform:
			print(colr.WARN + "\nExample of next cmd:" + colr.END + " amptk taxonomy -f %s -i %s -m mapping_file.txt -d ITS2\n" % (otu_new, final_table))
		else:
			print("\nExample of next cmd: amptk taxonomy -f %s -i %s -m mapping_file.txt -d ITS2\n" % (otu_new, final_table))
        
if __name__ == "__main__":
	main(args)