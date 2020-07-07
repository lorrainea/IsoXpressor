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

import pandas as pd
import matplotlib.pyplot as plt
import argparse
import os
import sys

def parseArgs():
	parser = argparse.ArgumentParser(description='Assessing transcriptional activity within isochores', epilog="", formatter_class=argparse.RawDescriptionHelpFormatter)
	parser.add_argument('-o', '--output_dir', type=str, required=True, help="Output directory.")
	parser.add_argument('-c', '--conditions', type=str, required=True, help="The conditions file identifying which read belongs to which condition.")
	parser.add_argument('-a', '--statistical_analysis', type=str, required=True, help="'TPM' (Transcripts Per Kilobase Million) or 'RPKM' (Reads Per Kilobase Million).")
	args = parser.parse_args()

	if args.statistical_analysis != "TPM" and args.statistical_analysis != "RPKM":
		print("-Error: Statistical analysis should be TPM or RPKM ")	
		exit()

	if not os.path.isfile(args.conditions):
		print("-Error: Conditions file ".join(args.conditions).join(" does not exist"))	
		exit()

	if not os.path.exists(args.output_dir):
		print("-Error: Output folder ".join(args.output_dir).join(" does not exist."))
		exit()

	return args


def chromosome_profiles(args):

	print("-Creating chromosome profile graphs")
	with open(os.path.join(args.output_dir,"avg_exp.csv")) as f:
	    lines = f.readlines()

	graphs = os.path.join(args.output_dir, "chromosome_profiles")
	os.mkdir(graphs)

	cond_count = 0
	with open(args.conditions) as f2:
		conditions = f2.readlines()

	for line in conditions:
		cond = line.split(',')
		if( int(cond[1]) > cond_count ):
			cond_count = int(cond[1])

	labels = [0 for i in range(cond_count)]
	for line in conditions:
		cond = line.split(',')
		labels[int(cond[1])-1] = cond[2].strip('\n')

	for chrom in os.listdir(os.path.join(args.output_dir, 'chromosomes')):
		h1 = pd.DataFrame(columns=['start','gc','size'])
		h2 = pd.DataFrame(columns=['start','gc','size'])
		h3 = pd.DataFrame(columns=['start','gc','size'])
		l1 = pd.DataFrame(columns=['start','gc','size'])
		l2 = pd.DataFrame(columns=['start','gc','size'])

		exp = pd.DataFrame()
		for k in range(0,cond_count):
			exp = exp.append(pd.DataFrame(columns=['exp'+str(k)] ), sort=False)	
		start = []


		maxval = 0;
		chromosome = str(chrom).split('.fast')[0].replace(',',' ')
		for j in range(1, len(lines)):
			if chromosome == lines[j].split(',')[0]:
				if lines[j].split(',')[1] == "H1":
					h1 = h1.append( {'gc': float(lines[j].split(',')[2]) , 'start': float(lines[j].split(',')[3])/1000000.0 , 'size': float(lines[j].split(',')[5])/1000000.0 }, ignore_index=True)
				elif lines[j].split(',')[1] == "H2":
					h2 = h2.append( {'gc': float(lines[j].split(',')[2]) ,'start': float(lines[j].split(',')[3])/1000000.0, 'size': float(lines[j].split(',')[5])/1000000.0 }, ignore_index=True)
				elif lines[j].split(',')[1] == "H3":
					h3 = h3.append( {'gc': float(lines[j].split(',')[2]) ,'start': float(lines[j].split(',')[3])/1000000.0 ,'size': float(lines[j].split(',')[5])/1000000.0 }, ignore_index=True)
				elif lines[j].split(',')[1] == "L1":
					l1 = l1.append( {'gc': float(lines[j].split(',')[2]) ,'start': float(lines[j].split(',')[3])/1000000.0 ,'size': float(lines[j].split(',')[5])/1000000.0 }, ignore_index=True)
				elif lines[j].split(',')[1] == "L2":
					l2 = l2.append( {'gc': float(lines[j].split(',')[2]) ,'start': float(lines[j].split(',')[3])/1000000.0 ,'size': float(lines[j].split(',')[5])/1000000.0 }, ignore_index=True)
				start.append(float(lines[j].split(',')[3])/1000000.0)

				vals = []
				for k in range(cond_count):
					value = float(lines[j].split(',')[6+k])
					vals.append( value )
					if( value > maxval ):
						maxval = value
				
				a_series = pd.Series(vals, index = exp.columns)
				exp = exp.append(a_series, ignore_index=True)
		

		fig, ax_arr  = plt.subplots(cond_count+1)

		if( not h3['gc'].empty ):
			ax_arr[cond_count].bar(h3['start'], h3['gc'], color='red', width=h3['size'], linewidth=0, align='edge', label="H3")
		if( not h2['gc'].empty ):
			ax_arr[cond_count].bar(h2['start'], h2['gc'], color='orange', width=h2['size'], linewidth=0, align='edge', label="H2")
		if( not h1['gc'].empty ):
			ax_arr[cond_count].bar(h1['start'], h1['gc'], color='yellow', width=h1['size'], linewidth=0, align='edge', label="H1")
		if( not l2['gc'].empty ):
			ax_arr[cond_count].bar(l2['start'], l2['gc'], color='dodgerblue', width=l2['size'], linewidth=0, align='edge', label="L2")	
		if( not l1['gc'].empty ):
			ax_arr[cond_count].bar(l1['start'], l1['gc'], color='lightskyblue', width=l1['size'], linewidth=0, align='edge', label="L1")
		

		colours = ['purple','green','brown','black','red'] #append more colours depending on conditions

		for k in range(0,cond_count):
			ax_arr[k].plot(start, exp['exp'+str(k)], color=colours[k], label=labels[k]) 
			ax_arr[k].legend(prop={'size': 16})
			ax_arr[k].xaxis.set_ticks([])
			ax_arr[k].tick_params(axis='y', which='major', labelsize=15)

		ax_arr[cond_count].tick_params(axis='y', which='major', labelsize=15)

		center = (cond_count + 1)/2
		ax_arr[center].set_ylabel("Average " + args.statistical_analysis+" Expression" , fontsize=20)
		custom_ylim = (0, maxval+ (5.0/100.0)*maxval+1)
		for k in range(0, cond_count):
			plt.setp(ax_arr[k], ylim=custom_ylim)

		ax_arr[cond_count].tick_params(axis='x', which='major', labelsize=15)
		ax_arr[cond_count].set_xlabel("Position (MB)", fontsize=23)
		ax_arr[cond_count].set_ylabel("GC\nLevel\n(%)", fontsize=18)
		ax_arr[cond_count].legend(prop={'size': 9.5}, loc='right')

		fig.set_size_inches(21,11)
		plt.subplots_adjust(hspace=0.19)
		plt.savefig(os.path.join( graphs, chrom+'_'+args.statistical_analysis+'.jpg'))
		plt.close()


def main():
	args = parseArgs()
	chromosome_profiles(args)
    
if __name__ == '__main__':
    main()
