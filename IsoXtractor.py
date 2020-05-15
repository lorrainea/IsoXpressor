#IsoXtractor
#Copyright (C) 2020 Lorraine A. K. Ayad
#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>

import os
import sys
import random
import string
import argparse
import math
import subprocess

import reads
import analysis

def parseArgs():
	parser = argparse.ArgumentParser(description='Assessing transcriptional activity within isochores', epilog="", formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-g', '--genome', type=str, required=True, help="Genome file name")
	parser.add_argument('-r', '--reads_dir', type=str, required=True, help="The directory name where the read files exist.")
	parser.add_argument('-o', '--output_dir', type=str, required=True, help="Output directory.")
	parser.add_argument('-c', '--conditions', type=str, required=True, help="The conditions file identifying which read belongs to which condition.")
	parser.add_argument('-w', '--window_size', type=int, required=False, help="Set window size of isochores in bp for IsoSegmenter (default: 100,000).")
	parser.add_argument('-s', '--seed_errors', type=int, required=False, help="Maximum number of seed errors for REAL (default: 2).")
	parser.add_argument('-e', '--total_errors', type=int, required=False, help="Total number of errors for REAL (default: 5).")
	parser.add_argument('-l', '--seed_length', type=int, required=False, help="Length of seed in bp for REAL (default: 32).")
	parser.add_argument('-a', '--statistical_analysis', type=str, required=False, help="'TPM' (Transcripts Per Kilobase Million) or 'RPKM' (Reads Per Kilobase Million) (default: 'TPM').")
	parser.add_argument('-t', '--threads', type=int, required=False, help="Number of threads to use (default: 1).")

	args = parser.parse_args()

	if args.window_size == None:
		args.window_size = 100000
	if args.seed_errors == None:
		args.seed_errors = 2
	if args.total_errors == None:
		args.total_errors = 5
	if args.seed_length == None:
		args.seed_length = 32
	if args.threads == None:
		args.threads = 1
	if args.statistical_analysis != "TPM" or args.statistical_analysis != "RPKM":
		args.statistical_analysis = "TPM"


	if not os.path.exists(args.genome):
		print("-Error: Genome file ".join(args.genome).join(" does not exist"))
		exit()
	else: 
		path = args.genome  # this is a path given by user of where genome folder is

	if not os.path.exists(args.reads_dir):
		print("-Error: Reads folder ".join(args.reads_dir).join(" does not exist"))
		exit()
	else:
		path_reads = args.reads_dir # path where reads are contained
		
	if not os.path.isfile(args.conditions):
		print("-Error: Conditions file ".join(args.conditions).join(" does not exist"))	
		exit()

	if os.path.exists(args.output_dir):
		print("-Error: Output folder ".join(args.output_dir).join(" exists. Please choose a different name."))
		exit()

	#check genome ends in .fa
	if( path.endswith(".fa") == 0 ):
		print( "-Error: Genome fasta filename must end in .fa")
		exit()

	return args

def splitGenome(path, output):
	#split genome into chromosomes and store in folder chromosomes
	os.makedirs(os.path.join(output, "chromosomes"))
	print("-Splitting genome into chromosomes")
	text_file = open(path, "r").read() #file path of genome provided by user, open filename inside directory
	seq = text_file.splitlines()
	for k in range(0 , len(seq)):
		if(seq[k].startswith(">")):
			file_in = str(seq[k][1:])
			chr_dir = os.path.join(output, "chromosomes")
			text_file2 = open(os.path.join(chr_dir, str(file_in)+".fasta"), "w")
			text_file2.close()
		text_file3 = open(os.path.join(chr_dir, str(file_in)+".fasta"), "a")
		text_file3.write(seq[k]+"\n")
		text_file3.close()


def runReal(path, path_reads, output, args):
	#run REAL on all reads
	aligned = os.makedirs(os.path.join(output, "aligned"))
	for filename_read in os.listdir(path_reads): 
		print("-Running REAL on ".join(str(filename_read)))
		real = './real'
		command = [real]
		command.append('-t')
		command.append(path)
		command.append('-s')
		command.append(str(args.seed_errors))
		command.append('-e')
		command.append(str(args.total_errors))
		command.append('-l')
		command.append(str(args.seed_length))
		command.append('-T')
		command.append(str(args.threads))
		command.append('-p')
		command.append(os.path.join(path_reads,filename_read))
		command.append('-o')
		command.append(os.path.join(os.path.join(output, "aligned"), str(filename_read)+".OUT"))
		subprocess.call(command)

def runIsosegmenter(output, args):
	chr_dir = os.path.join(output, "chromosomes")
	graphs = os.makedirs(os.path.join(output, "graphs"))
	isochores = os.makedirs(os.path.join(output, "isochores"))
	#run isosegmenter on all chromosomes 
	for filename in os.listdir(os.path.join(output, "chromosomes")):
		print("-Running isoSegmenter on "+str(filename))
		isosegmenter = './isoSegmenter/scripts/isoSegmenter.py'
		command = [isosegmenter]
		command.append('-i')
		command.append(os.path.join(chr_dir, str(filename)))
		command.append('--y_min')
		command.append('1')
		command.append('--y_max')	
		command.append('100')
		command.append('-g')
		command.append(os.path.join(os.path.join(output, "graphs"), str(filename)+".jpg"))
		command.append('--window_size')
		command.append(str(args.window_size))
		command.append('-o')
		command.append(os.path.join(os.path.join(output, "isochores"), str(filename)+".csv"))
		subprocess.call(command)	



def main():
	args = parseArgs()
	path = args.genome
	path_reads = args.reads_dir
	output = args.output_dir

	os.mkdir(output)

	splitGenome(path, output)
	runReal(path, path_reads, output, args)
	runIsosegmenter(output, args)
	reads.computeReadHead(output)
	reads.computeReads(args, output)
	reads.joinReads(output)
	lines = []
	cond_count = analysis.numConditions(args, lines)

	conditions =[[] for i in xrange(cond_count)]
	isochores = []
	isochore_family = []
	analysis.initialise(output, args, isochores, isochore_family, conditions, cond_count, lines)

	isochore_family = set( isochore_family )
	iso_family = list( isochore_family )
	conditions_val=[[] for i in xrange(cond_count)]
	total_reads = []
	analysis.computeConditions(conditions, conditions_val, isochores, total_reads, cond_count, lines)

	expression = [0] * cond_count
	final_table = [0] * ( 1+len(iso_family) )
	mean = [0] * len(iso_family)
	analysis.computeExpression(cond_count, total_reads, conditions_val, isochores, args, final_table, expression, iso_family, mean)
	analysis.avgExpression(output, isochores, expression, cond_count )
	analysis.chromAvgExpression(output, cond_count, isochores, expression)
	analysis.isoClassAvgExpression(final_table, iso_family, isochores, output, cond_count, expression, mean)

if __name__ == "__main__":
    main()
	
