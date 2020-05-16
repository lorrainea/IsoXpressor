#IsoXpressor
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
import multiprocessing
from functools import partial

def parseArgs():
	parser = argparse.ArgumentParser(description='Assessing transcriptional activity within isochores', epilog="", formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-o', '--output_dir', type=str, required=True, help="Output directory.")
	parser.add_argument('-t', '--threads', type=int, required=False, help="Number of threads to use (default: 1).")

	args = parser.parse_args()

	if args.threads == None:
		args.threads = 1

	if not os.path.exists(args.output_dir):
		print("-Error: Output folder ".join(args.output_dir).join(" does not exist."))
		exit()
	else: 
		output = args.output_dir  # this is a path given by user of where genome folder is

	return args

#os.path.join(os.path.join(output, "graphs"), str(filename)+".jpg"))
def align(aligned, output):
	aligned_file = open(os.path.join(os.path.join(output, 'aligned'), aligned), "r") #output from REAL
	print("-Computing reads for " + str(aligned) )
	temp_aligned = open(os.path.join(output, 'temp'+aligned+'.csv'), 'a')

	temp_aligned.write(aligned.split('.fasta')[0]+'\n'),

	for isochore in os.listdir(os.path.join(output, 'isochores')):
		aligned_iso = []
		print("-Computing reads for " + str(isochore) )	
		isochore_file = open(os.path.join(os.path.join(output, 'isochores'), isochore ), "r").read() #output from isoSegmenter
		isochores = isochore_file.splitlines()

		while True:
			line = aligned_file.readline().strip()

			if not line:
        			break

			re = line.split('\t')

			if( str(isochore).split('.fasta')[0] == re[7] ):
				aligned_iso.append(line)

			
		for i in range(1, len(isochores)):
			iso = isochores[i].split(',')
			counter = 0
			if ( iso[3] != "gap" ): 
				for j in range(0, len(aligned_iso)):
					re = aligned_iso[j].split('\t')

					if float(re[8]) <  float(iso[0]) and (float(re[8]) + float(re[5])  > float(iso[0])) and   ( ( float(re[8]) + float(re[5]) ) -  float(iso[0]) ) > ( float(re[8]) + float(re[5]) )/2:
						counter = counter +1 
					elif float(re[8]) > float(iso[0]) and ( ( float(re[8]) +  float(re[5]) ) < float(iso[1]) ):
						counter = counter +1
					elif float(re[8]) < float(iso[1]) and (float(re[8]) + float(re[5]) ) > iso[1] and float(iso[1]) - float(re[8]) > ( float(re[8]) + float(re[5]) )/2:
						counter=counter+1
				
				temp_aligned.write( str(counter)+'\n')	
				temp_aligned.flush()

		aligned_file.seek(0)

	temp_aligned.flush()

def computeReadHead(output):
	#print titles for number of reads file Chromosome
	temp1 = open(os.path.join(output, "temp1.csv"), 'a')
	temp1.write("Chromosome,Isochore Class,GC Level,Isochore Start,Isochore End,Isochore Size\n")

	#filling in column by column

	for isochore in os.listdir(os.path.join(output, 'isochores')):
		isochore_file = open(os.path.join(os.path.join(output, 'isochores'), isochore ), "r").read() #output from isoSegmenter
		isochores = isochore_file.splitlines()
			
		for i in range(1, len(isochores)):
			iso = isochores[i].split(',') 

			if ( iso[3] != "gap" ): 
				counter = 0
				if(',' in str(isochore) ):
					new_chrm = str(isochore).replace(',', ' ')
				else:
					new_chrm = str(isochore)
				temp1.write(new_chrm.split('.fasta')[0] +','+ iso[3]  +','+ iso[4]  +','+ iso[0] +','+ iso[1]  +','+ iso[2] +'\n' )


	temp1.flush()
	temp1.close()


def computeReads(args, output):
	pool = multiprocessing.Pool(args.threads) 
	align_partial = partial(align, output=output)
	pool.map(align_partial, os.listdir(os.path.join(output, 'aligned')))

	pool.close()
	pool.join()


def joinReads(output):
	for aligned in os.listdir(os.path.join(output,'aligned')):
		#combine temp files	
		#tempaligned.flush()
		tempOne = open(os.path.join(output,"temp1.csv"), 'r').read()	
		tempTwo = open(os.path.join(output,'temp'+aligned+'.csv'), 'r').read()
		temp3 = open(os.path.join(output,"temp3.csv"), 'w')

		lines = tempOne.splitlines()
		lines2 = tempTwo.splitlines()
		
		for i in range ( 0, len(lines) ):
			temp3.write( lines[i] + ',' + lines2[i] + '\n')
			temp3.flush()

		del(lines)
		del(lines2)

		os.remove(os.path.join(output,"temp1.csv"))
		os.remove(os.path.join(output,'temp'+aligned+'.csv'))
		os.rename(os.path.join(output,"temp3.csv"), os.path.join(output,"temp1.csv")) 
		

	os.rename(os.path.join(output,"temp1.csv"), os.path.join(output,"reads.csv")) 


def main():
	args = parseArgs()
	output = args.output_dir
	computeReadHead(output)
	computeReads(args, output)
	joinReads(output)
	
if __name__ == '__main__':
    main()

