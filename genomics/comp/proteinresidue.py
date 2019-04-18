import os
import sys
import time
import pprint
import subprocess
import numpy as np
import pprint
import pickle
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

from math import log

amino_acids_list = ["A", "C", "D", "E", "F", "G", "H", "I", "K", "L",\
		    "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"]

def parse_alignment_file(aligned_fasta_file):
	"""
	Go through an aligned FASTA file, extract all the sequences,
	determine which of the sequence positions have residues (i.e.,
	not gaps) in 50% or more of the sequences ("well-represented positions")
	Return a tuple (seqs_byposition_dict, human_sequence), where
	(1) seqs_byposition_dict is a dictionary
	that contains each of the well-represented positions and the
	residue at that position in each of the sequences:
	seqs_byposition_dict[pos][seq_name], and
	(2) human_sequence is a string of the aligned human sequence
	for later analysis
	"""
	sequences_dict = {}

	with open(aligned_fasta_file, "r") as f:
		# since sequences span multiple lines, we'll keep track of
		# one as we're going through in prev_line_seq, and then
		# deposit it in the dictionary when we hit the next line
		# starting with ">"
		prev_line_seq = ""
		for line in f:
			# if line begins with ">", it's the name of a sequence,
			# so we want to store the previous sequence in the
			# dictionary and start a new record for this sequence
			if (line.strip()[0] == ">"):
				if (len(prev_line_seq) > 0):
					# add the previous sequence to the dict
					sequences_dict[seq_name] = prev_line_seq
					prev_line_seq = ""
				seq_name = line.strip()[1:].split("/")[0]

			# if line doesn't begin with ">", it's sequence that
			# we want to keep track of
			else:
				prev_line_seq += line.strip()

	# add the final sequence to the dictionary
	sequences_dict[seq_name] = prev_line_seq

	# Now go through each of the positions in the alignment and determine
	# which positions have residues (i.e., non-gaps) in at least 50% of
	# sequences - henceforth referred to as "well-represented positions" -
	# we'll only consider these positions going forward
	seq_names = sequences_dict.keys()
	length_of_seqs = len(sequences_dict[seq_names[0]])
	seqs_byposition_dict = {}

	plotted_position = 0
	for pos in range(length_of_seqs):
		tot_positions_with_res = 0
		for seq_name in seq_names:
			# see if that sequence has an amino acid or gap at pos
			if sequences_dict[seq_name][pos] in amino_acids_list:
				tot_positions_with_res += 1
		if (float(tot_positions_with_res)/len(seq_names) > 0.5):
			seqs_byposition_dict[pos] = {}
			# go through and include each of the sequences for that
			# position
			for seq_name in seq_names:
				seqs_byposition_dict[pos][seq_name] =\
				    sequences_dict[seq_name][pos]

	# get the human_sequence ("CGL_HUMAN") for later analysis
	human_sequence = sequences_dict["CGL_HUMAN"]
	
	return (seqs_byposition_dict, human_sequence)

def get_mono_freqs_and_di_counts(seqs_byposition_dict):
	"""
	With an input seqs_byposition_dict, a dictionary returned by the parse_alignment_file() function.
	Return a tuple of (mono_freqs_dict, di_counts_dict), two dictionaries that have
        (1) the (mono-)amino acid frequencies at each of the well-
	represented positions (2) the di-amino acid counts at each pair 
        of well-represented positions. The format is di_counts_dict[pos][sec_pos][aa_1 + aa_2],
        where pos is the left column, sec_pos is the right column, aa_1+aa_2 is a string of the di-a.a.
	"""
	# a list of the "well-represented positions"
	positions_list = seqs_byposition_dict.keys()
	positions_list.sort()
	seq_names_list = seqs_byposition_dict[positions_list[0]].keys()
	
	# the dictionaries that contain the frequences of mono- and di-
	# amino acids, which will be returned (the _counts_ dictionaries are
	# temporary dictionaries to store counts, which will then be normalized
	# in into _freqs_)
	mono_freqs_dict = {}
	di_counts_dict = {}
	for pos_num, pos in enumerate(positions_list):
		print("Calculating amino acid frequencies for original MSA position ",\
		    pos," of", positions_list[-1],\
		    " (well-represented position #", pos_num, ")")
		mono_counts_dict = {}
		mono_freqs_dict[pos] = {}
		di_counts_dict[pos] = {}
		
		# initialize the count and frequency to be 0 for each amino acid
		for aa in amino_acids_list:
			mono_counts_dict[aa] = 0
			mono_freqs_dict[pos][aa] = 0
		# go through each of the sequences and add the counts
		non_gaps_at_pos = 0
		for seq_name in seq_names_list:
			aa = seqs_byposition_dict[pos][seq_name]
			# make sure that it's an amino acid (not a gap)
			if aa in amino_acids_list:
				non_gaps_at_pos += 1
				mono_counts_dict[aa] += 1
		# normalize the counts into frequencies
		for aa in amino_acids_list:
			mono_freqs_dict[pos][aa] =\
			    mono_counts_dict[aa]/float(non_gaps_at_pos)
		
		# Get the di-a.a. counts by iterating through
		# all positions after pos to the end as the second position of
		# the di-a.a.
		for sec_pos in positions_list[(pos_num+1):]:
			di_counts_dict[pos][sec_pos] = {}
			for aa_1 in amino_acids_list:
				for aa_2 in amino_acids_list:
					di_counts_dict[pos][sec_pos][aa_1+aa_2]=0
			
			# go through each of the sequences and get the aa pair
			# in positions (pos, sec_pos)
			for seq_name in seq_names_list:
				aa_1 = seqs_byposition_dict[pos][seq_name]
				aa_2 = seqs_byposition_dict[sec_pos][seq_name]
				# make sure that both are amino acids (not gaps)
				if (aa_1 in amino_acids_list) and\
					    (aa_2 in amino_acids_list):
					di_counts_dict[pos][sec_pos][aa_1+aa_2]+=1
	
	return (mono_freqs_dict, di_counts_dict)

def get_information_content(mono_freqs_dict):
	"""
	Calculate the information content at each position and return
	a list, info_content_list, that contains tuples
	[(orig_pos_in_seq, info), (orig_pos_in_seq, info), ...],
	where orig_pos_in_seq is the index in the original MSA
	(i.e., the keys of mono_nts_freqs_dict) and info is a float of the
	information content at that position
	"""
	info_content_list = []

	# get a sorted list of all of the "well-represented positions"
	positions_list = mono_freqs_dict.keys()
	positions_list.sort()

	# go through each of the positions and get the information content
	largest_info = -1
	for pos in positions_list:
		info_at_position = log(20, 2)
		for aa in amino_acids_list:
			p_aa = mono_freqs_dict[pos][aa]
			if (p_aa > 0):
				info_at_position += (p_aa * log(p_aa, 2))
		info_content_list.append((pos, info_at_position))
		if (info_at_position > largest_info):
			largest_info = info_at_position
	# go through and get all indices that have largest_info
	print("Highest information content at position is: ", largest_info)
	for well_represented_pos, tupl in enumerate(info_content_list):
		orig_pos = tupl[0]
		info = tupl[1]
		if (info == largest_info):
			print("\nHighest information at original MSA position ",\
			    orig_pos, " (well-represented position # ",\
			    well_represented_pos, ")\n")

	return info_content_list
				
def plot_info_content_list(info_content_list):
	"""
	Create 2 plots of the information content at each position, one indexed by the 
        original MSA position and one indexed by the well-represented position 
	matplotlib must be installed; otherwise, use the output from info_content_list to plot
        the information at each position in Excel or another
        package of your choice
	"""

	print("Making 2 plots of the information content at each position...")

	max_info = -1
	min_info = 1000
	for tupl in info_content_list:
		info = tupl[1]
		if (info > max_info):
			max_info = info
		if (info < min_info):
			min_info = info

	# makes a figure, axis, and plots the points on it relative to the
	# well-positioned index
	fig = plt.figure()
	ax = fig.add_subplot(111)
	#plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	OrRd = plt.cm.OrRd
	for plot_pos, info_tupl in enumerate(info_content_list):
		info = info_tupl[1]
		plt.scatter(plot_pos, info,\
				color=OrRd((info-min_info)/(max_info-min_info)),\
				s = 10, edgecolor = 'none')

	sm = plt.cm.ScalarMappable(cmap=plt.cm.OrRd,\
		norm=plt.normalize(vmin=min_info, vmax=max_info))
	# include a bar that labels the colors with their information content
	sm._A = []
	ticks_list = []
	ticks_labels = []
	for i in range(11):
		ticks_list.append( min_info + i*(max_info - min_info)/10. )
		ticks_labels.append( "{0:.2f}".format(min_info+\
			i*(max_info - min_info)/10.))
	pp = plt.colorbar(sm, ticks = ticks_list)
	pp.ax.get_yaxis().set_ticklabels(ticks_labels)

	# set the title and labels
	ax.set_title("Information scatterplot")	
	ax.set_xlim((0, len(info_content_list) + 1))
	ax.set_ylim((0, 1.1*max_info))
	ax.set_xlabel("Well-represented Position",fontsize=14)
	ax.set_ylabel("Information",fontsize=14)
	# this will save the figure in the directory from which this script
	# is run
	fig.savefig("information_plot_well_represented_positions.pdf")
	print("\n\tPlot saved as information_plot_well_represented_positions.pdf")

	# makes a figure, axis, and plots the points on it relative to the
	# original MSA positions
	fig = plt.figure()
	ax = fig.add_subplot(111)
	#plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	OrRd = plt.cm.OrRd
	for plot_pos, info_tupl in enumerate(info_content_list):
		pos = info_tupl[0]
		info = info_tupl[1]
		plt.scatter(pos, info,\
				color=OrRd((info-min_info)/(max_info-min_info)),\
				s = 10, edgecolor = 'none')

	sm = plt.cm.ScalarMappable(cmap=plt.cm.OrRd,\
		norm=plt.normalize(vmin=min_info, vmax=max_info))
	# include a bar that labels the colors with their information content
	sm._A = []
	ticks_list = []
	ticks_labels = []
	for i in range(11):
		ticks_list.append( min_info + i*(max_info - min_info)/10. )
		ticks_labels.append( "{0:.2f}".format(min_info +\
			i*(max_info - min_info)/10.) )
	pp = plt.colorbar(sm, ticks = ticks_list)
	pp.ax.get_yaxis().set_ticklabels(ticks_labels)

	# set the title and labels
	ax.set_title("Information scatterplot")	
	ax.set_xlim((0, pos + 1))
	ax.set_ylim((0, 1.1*max_info))
	ax.set_xlabel("Original MSA Position",fontsize=14)
	ax.set_ylabel("Information",fontsize=14)
	# this will save the figure in the directory from which this script
	# is run
	fig.savefig("information_plot_original_MSA_positions.pdf")
	print("\n\tPlot saved as information_plot_original_MSA_positions.pdf")
	plt.clf()

def get_MI_at_pairs_of_positions(di_counts_dict):
	"""
	Calculate the MI at each pair of positions (i, j) in the multiple
	sequence alignment, where i < j are indexes corresponding to positions
	in the original MSA (i.e., are keys of mono_freqs_dict)
	Return a dictionary, mutual_information_dict, where
	mutual_information_dict[i][j] = MI_at_ij_positions
	When calculating the mutual information for a pair of positions,
	you can scale the di-a.a. counts directly into frequencies to
	get the joint distribution f^(i,j)_x,y, and you should sum over
	this joint distribution to get the marginal distributions
	f^(i)_x and f^(j)_y
	"""

	mutual_information_dict = {}
	
	
	print("\nMaking mutual_information_dict...")
	max_MI = -1
	positions_list = di_counts_dict.keys()
	positions_list.sort()
	length_of_seqs = len(positions_list)
	for pos_num, pos in enumerate(positions_list):
		# don't include the last a.a. since there's no a.a. following it
		# to be the 2nd a.a. of the di-a.a.
		if (pos_num != (length_of_seqs-1)):
			mutual_information_dict[pos] = {}
		for sec_pos in positions_list[(pos_num+1):]:
			mi = 0
			# for a given pair of positions, we calculate 3 things:
			#	1. the di-a.a. frequencies,
			#	2. the mono-a.a. frequencies at the first pos.,
			#	1. the mono-a.a. frequencies at the second pos.,
			di_freq = {}
			mono_freq_1st = {}
			mono_freq_2nd = {}

			# get the total number of dinucleotides
			total_di = 0					
			for aa_1 in amino_acids_list:
				for aa_2 in amino_acids_list:
					total_di += di_counts_dict[pos][sec_pos][aa_1+aa_2]

			# get the frequency of mono-a.a. in the first position
			for aa_1 in amino_acids_list:
				aa_1_count = 0
				for aa_2 in amino_acids_list:
					aa_1_count+=di_counts_dict[pos][sec_pos][aa_1+aa_2]
				try:
					mono_freq_1st[aa_1] = aa_1_count/float(total_di)
				except ZeroDivisionError:
					mono_freq_1st[aa_1] = 0.			

			# get the frequency of mono-a.a. in the second position
			for aa_2 in amino_acids_list:
				aa_2_count = 0
				for aa_1 in amino_acids_list:
					aa_2_count+=di_counts_dict[pos][sec_pos][aa_1+aa_2]
				try:
					mono_freq_2nd[aa_2] = aa_2_count/float(total_di)
				except ZeroDivisionError:
					mono_freq_2nd[aa_2] = 0.			

			# now that we have all of the frequencies, calculate the MI
			for aa_1 in amino_acids_list:
				for aa_2 in amino_acids_list:
					try:
						di_freq = di_counts_dict[pos][sec_pos][aa_1+aa_2]/\
							float(total_di)
						mono_1_freq = mono_freq_1st[aa_1]
						mono_2_freq = mono_freq_2nd[aa_2]
						# this term only contributes to the MI if
						# all terms are nonzero
						if ((mono_1_freq > 0) and\
							(mono_2_freq > 0) and\
							(di_freq > 0)):
							mi += (di_freq *\
							log(di_freq/(mono_1_freq*mono_2_freq),2))
					except ZeroDivisionError:
						pass
			mutual_information_dict[pos][sec_pos] = mi
			if (mi >= max_MI):
				max_MI = mi
	# go through and get the positions pairs where the maximum MI occurs
	list_of_indices_max_MI = []
	for pos_num, pos in enumerate(positions_list):
		for sec_pos in positions_list[(pos_num+1):]:
			if (mutual_information_dict[pos][sec_pos] == max_MI):
				list_of_indices_max_MI.append( (pos, sec_pos) )
	print("Maximum MI is ", max_MI)
	print("\t-The maximum MI occurs at positions:",
	pprint.pprint(list_of_indices_max_MI))
	
	return mutual_information_dict

def get_highest_MI_block_of_10(mutual_information_dict):
	"""
	Find the block of 10 consecutive position pairs (i.e., keys of
	mutual_information_dict) with the highest average MI:
	between (pos_1, pos_2)  (pos_2, pos_3) ... ... (pos_9, pos_10)		
	Return the list [pos_1, pos_2, ..., pos_10], where each of the pos_i
        are consecutive keys of mutual_information_dict
	"""

	highest_average_positions_list = []

	positions_list = mutual_information_dict.keys()
	positions_list.sort()
	
	# get all except the last 9 as candidate start positions
	# for the highest average block
	possible_blockof10_starts = positions_list[:-9]
	
	highest_total = -1.
	for pos_in_list, start in enumerate(possible_blockof10_starts):
		total = 0
		for i in range(9):
			first_pos = positions_list[pos_in_list + i]
			sec_pos = positions_list[pos_in_list + i + 1]
			total += mutual_information_dict[first_pos][sec_pos]
		# update the highest_average_positions_list if the block
		# starting at start is higher than any previous block
		if (total > highest_total):
			highest_total = total
			highest_average_positions_list = []
			for i in range(10):
				highest_average_positions_list.append(\
					positions_list[pos_in_list+i])
	
	return highest_average_positions_list

def plot_mutual_information(mutual_information_dict):
	"""
	Create a plot of the mutual information
	at each pair of positions.
	"""
	import matplotlib
	matplotlib.use('agg')
	import matplotlib.pyplot as plt
	from matplotlib.patches import Rectangle
	positions_list = mutual_information_dict.keys()
	positions_list.sort()
	min_pos = min(positions_list)
	max_pos = max(positions_list)
	num_positions = len(positions_list)

	array_to_plot = []
	for i in range(num_positions):
		temp_list = [0] * num_positions
		array_to_plot.append(temp_list)
	
	# get the maximum and minimum MIs in the dictionary
	min_mi = 100000
	max_mi = -1
	for pos_num, pos in enumerate(positions_list):
		for sec_pos_num, sec_pos in enumerate(positions_list):
			try:
				mi = mutual_information_dict[pos][sec_pos]
				array_to_plot[pos_num][sec_pos_num] = mi
				if (mi < min_mi):
					min_mi = mi
				if (mi > max_mi):
					max_mi = mi
			except KeyError:
				if (pos != sec_pos):
					mi = mutual_information_dict[sec_pos][pos]
					array_to_plot[pos_num][sec_pos_num] = mi
	# makes a figure, axis, and plots the mutual information between
	# residues i and j as a box at (i, j) with the color of the box
	# indicating the Mutual Information; here (i, j) are the #s
	# corresponding to what well-represented positions it is
	print("\nPlotting indexed relative to well-represented positions...")
	fig = plt.figure();
	ax = fig.add_subplot(111);
	#plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	this_cmap = plt.cm.Reds
	plt.imshow(np.array(array_to_plot), cmap=this_cmap, vmin=0., vmax=max_mi)
	ax.set_aspect("equal")
	sm = plt.cm.ScalarMappable(cmap=this_cmap,\
			norm=plt.normalize(vmin=0., vmax=max_mi));
	# include a bar indicating what colors correspond to what MI values
	sm._A = []
	ticks_list = []
	ticks_labels = []	
	for i in range(11):
		ticks_list.append( i*(max_mi-min_mi)/10. )
		ticks_labels.append( "{0:.2f}".format(i*(max_mi-min_mi)/10.) )
	pp = plt.colorbar(sm, ticks = ticks_list)
	pp.ax.get_yaxis().set_ticklabels(ticks_labels)

	ax.set_title("Mutual information heatmap")	
	ax.set_ylim((0, len(positions_list)+1))
	ax.set_xlim((0, len(positions_list)+1))
	ax.set_xlabel("Well-represented Position 1",fontsize=14)
	ax.set_ylabel("Well-represented Position 2",fontsize=14)
	# this will save the figure in the directory from which this script is run
	fig.savefig("mutual_information_plot_well_represented_positions.pdf")
	print("\n\tPlot saved as mutual_information_plot_well_represented_positions.pdf")
	plt.clf()

	## does the same thing, but plots according to the original position in
	## the multiple sequence alignment, not which "well-represented position"
	## it is
	array_to_plot = []
	for i in range(max_pos+1):
		temp_list = [0] * (max_pos+1)
		array_to_plot.append(temp_list)
	
	# get the maximum and minimum MIs in the dictionary
	min_mi = 100000
	max_mi = -1
	for pos_num, pos in enumerate(positions_list):
		for sec_pos_num, sec_pos in enumerate(positions_list):
			try:
				mi = mutual_information_dict[pos][sec_pos]
				array_to_plot[pos][sec_pos] = mi
				if (mi < min_mi):
					min_mi = mi
				if (mi > max_mi):
					max_mi = mi
			except KeyError:
				if (pos != sec_pos):
					mi = mutual_information_dict[sec_pos][pos]
					array_to_plot[pos][sec_pos] = mi
	# makes a figure, axis, and plots the mutual information between
	# residues i and j as a box at (i, j) with the color of the box
	# indicating the Mutual Information; here (i, j) are the #s
	# corresponding to what well-represented positions it is
	print("\nPlotting indexed relative to original MSA positions...")
	fig = plt.figure();
	ax = fig.add_subplot(111);
	#plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	this_cmap = plt.cm.Reds
	plt.imshow(np.array(array_to_plot), cmap=this_cmap, vmin=0., vmax=max_mi)
	ax.set_aspect("equal")
	sm = plt.cm.ScalarMappable(cmap=this_cmap,\
			norm=plt.normalize(vmin=0., vmax=max_mi));
	# include a bar indicating what colors correspond to what MI values
	sm._A = []
	ticks_list = []
	ticks_labels = []	
	for i in range(11):
		ticks_list.append( i*(max_mi-min_mi)/10. )
		ticks_labels.append( "{0:.2f}".format(i*(max_mi-min_mi)/10.) )
	pp = plt.colorbar(sm, ticks = ticks_list)
	pp.ax.get_yaxis().set_ticklabels(ticks_labels)

	ax.set_title("Mutual information heatmap")	
	ax.set_ylim((min_pos, max_pos + 1))
	ax.set_xlim((min_pos, max_pos + 1))
	ax.set_xlabel("Original MSA Position 1",fontsize=14)
	ax.set_ylabel("Original MSA Position 2",fontsize=14)
	# this will save the figure in the directory from which this script
	# is run
	fig.savefig("mutual_information_plot_original_MSA_positions.pdf")
	print("\n\tPlot saved as mutual_information_plot_original_MSA_positions.pdf")

if __name__ == "__main__":
	### the FASTA MSA file should be the command line argument
	aligned_fasta_file = sys.argv[1]

	### parse the FASTA MSA file to get back a dictionary which gives
	### the positions that have residues (non-gaps) in 50%
	### or more of the sequences and the residue at each of those positions
	### in all of the sequences
	###    - also returns the human sequence in the MSA for future analysis
	seqs_byposition_dict, human_sequence = \
	    parse_alignment_file(aligned_fasta_file)

	### get the frequency of a.a.s and counts of di-a.a.s - see the function
	### for the structure of how these frequencies are represented in the
	### dictionaries that are returned
	mono_freqs_dict, di_counts_dict =\
	get_mono_freqs_and_di_counts(seqs_byposition_dict)

	### get the information content at each position
	info_content_list = get_information_content(mono_freqs_dict)

	### you can uncomment the line below to plot the
	### information content at each position if you have matplotlib
	### installed; otherwise, plot info_content_list with some other
	### tool of your choice
	plot_info_content_list(info_content_list)

	### Calculate the mutual information between all pairs of
	### "well-represented positions"
	mutual_information_dict = get_MI_at_pairs_of_positions(di_counts_dict)
	
	### This will "pickle" mutual_information_dict (this is Python-speak
	### for saving the dictionary in binary format - "wb" corresponds
	### to Write in Binary mode - as mutual_information_dict.pkl)
	### so you can later load it without having to
	### run the above commands to regenerate it
	with open("mutual_information_dict.pkl", "wb") as f:
		pickle.dump(mutual_information_dict, f)

	### This will make 2 heatmap plots of the mutual information
	### dictionary, one indexed by the
	### original positions in the MSA and a condensed one indexed by the
	### "well-represented positions"
	### - only uncomment and use this if you have matplotlib installed
	plot_mutual_information(mutual_information_dict)

	### determine which residues in the human sequence correspond
	### to the 10 consecutive positions with the highest avg. MI

	### if you have previously pickled mutual_information_dict,
	### you can uncomment this to load it without having to run
	### the above functions!
	#with open("mutual_information_dict.pkl") as f:
	#	mutual_information_dict = pickle.load(f)

	### get the 10 consecutive positions with the highest avg. MI
	highest_MI_block_of_10 = get_highest_MI_block_of_10(mutual_information_dict)

	### This will print out
	### the MSA entries of the human sequence corresponding to the range
	### implied by the returned highest_MI_block_of_10
	if (len(highest_MI_block_of_10) == 10):
		start_of_MI_block = highest_MI_block_of_10[0]
		end_of_MI_block = highest_MI_block_of_10[-1]
		print("\nResidues ",start_of_MI_block,"-",end_of_MI_block,\
			" of the human sequence in the MSA are:")
		print(human_sequence[start_of_MI_block:(end_of_MI_block+1)])

		print("\nNow go figure out where these residues are in the ungapped")
		print("human protein!")
	else:
		print("Fill in get_highest_MI_block_of_10() !!")


