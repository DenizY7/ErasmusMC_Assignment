import re
import pandas as pd
import argparse

from Bio import SeqIO


def main(args):
	# Load and sort the motifs based on their length (longest to shortest)
	seq_motifs    = load_motifs(args.motifs)
	sorted_motifs = sorted(seq_motifs, key=len, reverse=True) 

	# Generate short-hand annotation of each target sequence in the fasta file generated by the retrieve_target_sequences.py script
	target_str_sequences  = extract_strs(sorted_motifs, args.fasta_file)
	# Generate short-hand annotation report as dataframe & save to file
	short_hand_sequencess = generate_STR_records(target_str_sequences, args.output)


def load_motifs(path):
	with open(path) as file:
		seq_motifs = file.readlines()
		seq_motifs = [x.strip() for x in seq_motifs] # Strip off line breaks

	return seq_motifs


def get_compl_seq(primer_seq):
	# Retrieve complement sequences of the primers

	# Define the base dictionary to retrieve complements
	base_dictionary = {
		'A': 'T',
		'T': 'A',
		'C': 'G',
		'G': 'C'
	}

	comp_strand = []
	for i in primer_seq:
		letter = base_dictionary[i]
		comp_strand.append(letter)

	compl_seq = ''.join(comp_strand)

	return compl_seq


def extract_strs(sorted_motifs, fasta_file_path):
	# Convert target sequences into short-hand annotation using STRs

	# List to collect short-hand annotation of target sequences
	target_str_sequences = []

	# Identify STRs in each target sequence using the motifs from sorted_motifs
	for record in SeqIO.parse(fasta_file_path, "fasta"):
		print(f"Converting target sequence ID {record.id} into STRs")
		target_seq = str(record.seq)
	
		# Create a list to collect all STR elements with their counts from the target sequence
		str_sequences = []

		n = 0 # Keep track of non-repetitive elements (Ns) in the sequence
		while len(target_seq) > 0:
			found_match = False

			for motif in sorted_motifs:
				compl_motif = get_compl_seq(motif) # Retrieve complementary seq of the motif

				# Motif list contains complementary seqs for 2bp motifs, thefore only use complementary motif in the search if motif is longer than 2bp
				if len(motif) > 2:
					valid_motif   = re.compile(rf"^({motif}|{compl_motif})+") # generate the regex motif for the search
				else: 
					valid_motif   = re.compile(rf"^({motif})+") # generate the regex motif for the search
				
				# Search the regex motif at the beginning of the target sequence
				matched_motif = valid_motif.match(target_seq)

				# When a match is found, retrieve the found motif and its count
				if matched_motif is not None:
					found_str = matched_motif.group()
					str_count = max((found_str.count(motif), found_str.count(compl_motif)))

					# If STR is repeated more than once, append it to the list along with its count
					if str_count > 1:
						# Check if any Ns were found, and append them to the list along with their count 
						if n > 0:
							str_sequences.append(f"[N]{n}")
							n = 0

						# Add the found motif (STR) to the list with its count
						str_sequences.append(f"[{motif}]{str_count}")

						# Remove the found motif section from the sequence, continue the search
						target_seq = target_seq[len(found_str):]
						found_match = True
						break # To continue the search from the beginning of the motif list
			
			# If no motif matches are found, remove one base (store as N), and
			# repeat the search until a new motif is found
			if not found_match:
				n += 1
				target_seq = target_seq[1:]

		# If target sequence ends with Ns, append them to the list at the end with their count
		if n > 0:
			str_sequences.append(f"[N]{n}")

		# Join the list elements to form the STR elements, and append to the main list
		target_str_seq = ''.join(str_sequences)
		target_str_sequences.append((str(record.seq), target_str_seq))

	return target_str_sequences


def generate_STR_records(target_str_sequences, output_path):
	# Put the found STRs into dataframe & save to file
	short_hand_seqs = pd.DataFrame(target_str_sequences, columns = ["target_sequence", "short_hand_annotation"])

	# Save to file
	short_hand_seqs.to_csv(f"{output_path}/short_hand_annotations.csv")

	return short_hand_seqs


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description = 'Retrieve short-hand annotations')
	parser.add_argument('-f','--fasta_file', help='path to the fasta file location', default="Output/target_sequences.fasta")
	parser.add_argument('-m','--motifs', help='path to the list of STR motifs', default="Data/Motifs.txt")
	parser.add_argument('-o', '--output', help='path to the output folder', default="Output")
	args = parser.parse_args()
	
	main(args)
