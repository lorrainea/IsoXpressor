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

def parseArgs():
	parser = argparse.ArgumentParser(description='Assessing transcriptional activity within isochores', epilog="", formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-o', '--output_dir', type=str, required=True, help="Output directory.")
	parser.add_argument('-c', '--conditions', type=str, required=True, help="The conditions file identifying which read belongs to which condition.")
	parser.add_argument('-a', '--statistical_analysis', type=str, required=False, help="'TPM' (Transcripts Per Kilobase Million) or 'RPKM' (Reads Per Kilobase Million) (default: 'TPM').")
	args = parser.parse_args()

	if args.statistical_analysis != "TPM" and args.statistical_analysis != "RPKM":
		args.statistical_analysis = "TPM"


	if not os.path.isfile(args.conditions):
		print("-Error: Conditions file ".join(args.conditions).join(" does not exist"))	
		exit()


	if not os.path.exists(args.output_dir):
		print("-Error: Output folder ".join(args.output_dir).join(" does not exist."))
		exit()
	else: 
		output = args.output_dir  # this is a path given by user of where genome folder is

	return args


def numConditions(args, lines):
	conditions_file = open(args.conditions)

	cond_count = 0
	for line in conditions_file:
		lines.append(line)
		no_newline = line.strip('\n')
		cond = line.split(',')
		if( int(cond[1]) > cond_count ):
			cond_count = int(cond[1])
	conditions_file.close()

	return cond_count

def initialise(output, args, isochores, isochore_family, conditions, cond_count, lines):

	for i in range( 0, len(lines) ):
		no_newline = lines[i].strip('\n')
		cond = no_newline.split(',')
		conditions[ int(cond[1])-1 ].append( cond[0] )

	input = open(os.path.join(output,"reads.csv"))
	inputList = input.read()

	seq = inputList.splitlines()


	for i in range(0, len(seq)):
		isochores.append ( seq[i].split(',') )


	#compute number of isochore families
	for i in range(1, len(isochores) ):
		isochore_family.append ( isochores[i][1] ) #isochore family must be in column 1 of input

	isochore_family = set( isochore_family )

	iso_family = list( isochore_family )

################################################################################################

#identify which read header is in which column
def computeConditions(conditions, conditions_val, isochores, total_reads, cond_count, lines):
	#conditions_val=[[] for i in xrange(cond_count)]

	for i in range(0, len(conditions) ):
		for j in range( 0, len(conditions[i] ) ):
			conditions_val[i].append( 0 )

	#head = seq[0]
	header = isochores[0][6:]

	for j in range(0, cond_count):
		for k in range( 0, len(conditions[j]) ):
			for l in range( 0, len(header) ):
				if( str( conditions[j][k] ) == str( header[l] ) ):
					conditions_val[j][k] = l + 6


	#compute total reads for each column
	for j in range ( 6, 6+len(lines) ):
		total = 0
		for k in range( 1, len( isochores) ):
			total = total + float(isochores[k][j])
		total_reads.append(total)

###################################################################################


def computeExpression(cond_count, total_reads, conditions_val, isochores, args, final_table, expression, iso_family, mean, lines):
	#expression = [0] * cond_count #columns

	for j in range(0, cond_count ):
		expression[j] = [0] * len(isochores) #rows


	#compute average for each of the conditions
	if args.statistical_analysis == "NORM":
		for j in range(0, cond_count):
			for l in range (1, len(isochores) ):
				avg = 0
				#total_reads_counter = j*cond_count
				for k in range( 0, len(conditions_val[j]) ):
					formula = float(isochores[l][conditions_val[j][k]]) / ( total_reads[int(conditions_val[j][k])-6] * float(isochores[l][5]) ) 
					#lengths must be in column 5 of input

					if formula == 0:
						avg = avg + 0
					else: avg = avg + ( math.log( formula, 2.0 ) ) # compute log of formula above

				expression[j][l-1] = avg/ len(conditions_val[j])

	elif args.statistical_analysis == "TPM":
		for j in range(0, cond_count):
			rpk = 0
			for l in range (1, len(isochores) ):
				#total_reads_counter = j*cond_count
				for k in range( 0, len(conditions_val[j]) ):
					rpk = rpk +  float(isochores[l][conditions_val[j][k]]) / ( float(isochores[l][5]) / 1000.0 ) 

			rpk = float(rpk) / 1000000.0
			for l in range (1, len(isochores) ):
				avg = 0
				#total_reads_counter = j*cond_count
				for k in range( 0, len(conditions_val[j]) ):
					tpm = ( float(isochores[l][conditions_val[j][k]]) / ( float(isochores[l][5]) / 1000.0 ) ) / float(rpk)
					avg = avg + tpm

				expression[j][l-1] = avg/ len(conditions_val[j])

	elif args.statistical_analysis == "RPKM":
		for j in range ( 0, len(lines) ):
			total_reads[j] = float(total_reads[j]) / 1000000.0

		
		for j in range(0, cond_count):
			rpkm = 0
			for l in range (1, len(isochores) ):
				avg = 0
				#total_reads_counter = j*cond_count
				for k in range( 0, len(conditions_val[j]) ):
					rpkm = ( float(isochores[l][conditions_val[j][k]]) / float( total_reads[int(conditions_val[j][k])-6] ) ) / float(isochores[l][5])
					avg = avg + rpkm

				expression[j][l-1] = avg/ len(conditions_val[j])


	for j in range(0, 1+len(iso_family) ):
		final_table[j] = [0] * ( 3 * cond_count + 2 )

	#add titles to final table
	final_table[0][0] = "isochore family"
	c = 0
	inp = 1
	for j in range(0, 3*cond_count + 1 ):
		cond = "Condition "
		cond+= str(inp)
		if c == 0: cond+= " avg"
		if c == 1: cond+= " std dev"
		if c == 2: cond+= " var"
		if c == 2: 
			inp = inp + 1
			c = 0
		else: c = c + 1
		final_table[0][j+1] = cond

	final_table[0][3 * cond_count + 1 ] = " Count"

	for j in range(0, len(iso_family) ):
		mean[j] = [0] * cond_count

	#compute mean
	for i in range( 0, len(iso_family) ):
		for j in range(0, cond_count ):
			counter = 0
			total = 0
			for k in range( 1, len( isochores ) ):
				if isochores[k][1] == iso_family[i]: #isochore family in column 1
					counter = counter + 1
					total = total + expression[j][k-1]

			mean[i][j] = total / counter


#########################table 2 ####################################################

def avgExpression(output, isochores, expression, cond_count):
	print("-Computing average expression for isochores.")


	input_read = open(os.path.join(output,"reads.csv"))
	input_readlist = input_read.read()

	all_reads = input_readlist.splitlines()
	avg_exp = open(os.path.join(output,"avg_exp.csv"), 'w')

	split = all_reads[0].split(',')

	for j in range(0, 6 ):
		avg_exp.write( split[j] + ',' )
	for k in range(0, cond_count):
		avg_exp.write( "condition " + str(k+1) + ',' )
	avg_exp.write('\n')	


	for i in range(1, len(isochores) ) :
		split = all_reads[i].split(',')
		for j in range(0, 6 ):
			avg_exp.write( split[j] + ',' )
		for k in range(0, cond_count):
			avg_exp.write( str(expression[k][i-1]) + ',' )
			
		avg_exp.write('\n')


#########################table 3 ####################################################


def chromAvgExpression(output, cond_count, isochores, expression):
	print("-Computing average expression for chromosomes.")


	input_read = open(os.path.join(output,"avg_exp.csv"))
	input_readlist = input_read.read()

	all_reads = input_readlist.splitlines()


	chromosomes = []
	for i in range(1, len(all_reads) ) :
		if  all_reads[i].split(',')[0] not in chromosomes:
			chromosomes.append(all_reads[i].split(',')[0])

	final_table_chr = [0] * ( 1+len(chromosomes) )

	for j in range(0, 1+len(chromosomes) ):
		final_table_chr[j] = [0] * ( 3 * cond_count + 2 )

	#add titles to final table
	final_table_chr[0][0] = "chromosome"
	c = 0
	inp = 1
	for j in range(0, 3*cond_count + 1 ):
		cond = "Condition "
		cond+= str(inp)
		if c == 0: cond+= " avg"
		if c == 1: cond+= " std dev"
		if c == 2: cond+= " var"
		if c == 2: 
			inp = inp + 1
			c = 0
		else: c = c + 1
		final_table_chr[0][j+1] = cond

	final_table_chr[0][3 * cond_count + 1 ] = " Count"

	mean_chr = [0] * len(chromosomes)
	for j in range(0, len(chromosomes) ):
		mean_chr[j] = [0] * cond_count

	#compute mean
	for i in range( 0, len(chromosomes) ):
		for j in range(0, cond_count ):
			counter = 0
			total = 0
			for k in range( 1, len( isochores ) ):
				if isochores[k][0] == chromosomes[i]: #isochore family in column 1
					counter = counter + 1
					total = total + expression[j][k-1]

			if counter == 0:
				mean_chr[i][j] = total
			else:
				mean_chr[i][j] = total / counter


	for i in range(0, len(chromosomes)):
		final_table_chr[i+1][0] = chromosomes[i]
		cond = 0
		for j in range(0, cond_count ):
			counter = 0
			standard_dev_chr = 0
			variance = 0
			total = 0
			for k in range( 1, len( isochores ) ):
				if isochores[k][0] == chromosomes[i]: #isochore family in column 1
					counter = counter + 1
					total = total + expression[j][k-1]
					standard_dev_chr = standard_dev_chr + ( ( expression[j][k-1] - mean_chr[i][j] )** 2 )


			final_table_chr[i+1][cond + 1] = mean_chr[i][j] #average

			if counter == 0:
				final_table_chr[i+1][cond + 2] = math.sqrt( standard_dev_chr )
				final_table_chr[i+1][cond + 3] = standard_dev_chr
			else:
				final_table_chr[i+1][cond + 2] = math.sqrt( standard_dev_chr / counter )
				final_table_chr[i+1][cond + 3] = standard_dev_chr / counter
			final_table_chr[i+1][3*cond_count + 1] = counter #count
			cond = cond + 3

	avg_exp = open(os.path.join(output,"avg_exp_chromosome.csv"), 'w')

	for i in range (1+len(chromosomes) ):
		for j in range(0, 3 * cond_count + 2):
			avg_exp.write( str(final_table_chr[i][j])+',' )	
		avg_exp.write("\n")
			

###############################################################################

def isoClassAvgExpression(final_table, iso_family, isochores, output, cond_count, expression, mean):
	print("-Computing average expression for isochore classes.")	

	#compute final table
	for i in range( 0, len(iso_family) ):
		final_table[i+1][0] = iso_family[i]
		cond = 0
		for j in range(0, cond_count ):
			counter = 0
			standard_dev = 0
			variance = 0
			total = 0
			for k in range( 1, len( isochores ) ):
				if isochores[k][1] == iso_family[i]: #isochore family in column 1
					counter = counter + 1
					total = total + expression[j][k-1]
					standard_dev = standard_dev + ( ( expression[j][k-1] - mean[i][j] )** 2 )


			final_table[i+1][cond + 1] = mean[i][j] #average
			final_table[i+1][cond + 2] = math.sqrt( standard_dev / counter )
			final_table[i+1][cond + 3] = standard_dev / counter
			final_table[i+1][3*cond_count + 1] = counter #count
			cond = cond + 3


	exp = open(os.path.join(output,"expression.csv"), 'w')
	for i in range ( 1+len(iso_family) ):
		for j in range(0, 3 * cond_count + 2):
			exp.write( str(final_table[i][j])+',' )	
		exp.write("\n")
			

def main():
	args = parseArgs()
	output = args.output_dir
	conditions = [[]]
	isochores = []
	lines = []
	isochore_family = []
	cond_count = numConditions(args, lines)

	conditions =[[] for i in xrange(cond_count)]

	initialise(output, args, isochores, isochore_family, conditions, cond_count, lines)
	
	isochore_family = set( isochore_family )
	iso_family = list(isochore_family)
	conditions_val=[[] for i in xrange(cond_count)]	
	total_reads = []
	computeConditions(conditions, conditions_val, isochores, total_reads, cond_count, lines)
	mean = [0] * len(iso_family)
	expression = [0] * cond_count
	final_table = [0] * ( 1+len(iso_family) )
	computeExpression(cond_count, total_reads, conditions_val, isochores, args, final_table, expression, iso_family, mean, lines)
	avgExpression(output, isochores, expression, cond_count )
	chromAvgExpression(output, cond_count, isochores, expression)
	isoClassAvgExpression(final_table, iso_family, isochores, output, cond_count, expression, mean)
    
if __name__ == '__main__':
    main()
