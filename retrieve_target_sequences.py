import argparse
import pandas as pd
import pickle

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def main(args):
	# Get the list of primers
	primers_df = get_primers(args.primers)

	# Set the max target seq length (1000bp)
	seq_len = args.seq_length

	#Retrieve possible target sequences from the reference genome
	found_sequences = get_list_found_seqs(primers_df, seq_len, args.reference_genome)

	# # Save the found sequences into file
	# with open(f"{args.output}/found_seqs.pickle", "wb") as handle:
	#     pickle.dump(found_sequences, handle)

	# # Load found sequences
	# with open(f"{args.output}/found_seqs.pickle", 'rb') as handle:
	# 	found_sequences = pickle.load(handle)

	# Extract the amplicons from the possible target sequences found in the genome
	found_amplicons = extract_amplicon_sequence(primers_df, found_sequences)

	# Generate amplicon report as dataframe & save to file
	all_amplicons = generate_amplicon_records(found_amplicons, args.output)
	# Create fasta file of the found amplicons & save to file
	generate_targetSeq_fasta(all_amplicons, args.output)

	
def get_primers(path):
	# Load the primers file as dataframe
	primers_df = pd.read_csv(path, sep="\t", names=["number", "primer"])
	# Get the reverse complement primer sequences, named as compl_primer
	primers_df['compl_primer'] = primers_df.apply(lambda x: get_revcompl_seq(x.primer), axis=1)

	return primers_df


def get_revcompl_seq(primer_seq):
	# Retrieve reverse complement sequences of the primers, aka, possible reverse primers
	
	# Define the base dictionary to retrieve complements
	base_dictionary = {
		'A': 'T',
		'T': 'A',
		'C': 'G',
		'G': 'C',
		'R': 'Y',
		'Y': 'R',
	}

	comp_strand = []
	for i in primer_seq:
		letter = base_dictionary[i]
		comp_strand.append(letter)

	rev_compl_seq = ''.join(reversed(comp_strand))

	return rev_compl_seq


def get_list_found_seqs(primers_df, seq_len, reference_genome_path):
	# Retrieve possible target sequences from the Reference Genome
	# Search only in one direction (forward primer + 1000bp), reverse primers will be search in these possible target sequences

	found_sequences = {}
	
	for idx, primer, compl_primer in primers_df.values:
		found_sequences[idx] = extract_sequences(primer, seq_len, reference_genome_path)

	return found_sequences


def extract_sequences(primer, seq_len, reference_genome_path):
	# Function to search reference genome for primer(s) and extract possible target sequences
	print(f"Extracting sequences for the primer: {primer}")
	
	# Create empty list to store genomic location and sequence for possible target sequence
	possible_sequences = []
	
	# Check if a primer has nucleotide in compact IUPAC code and generate all possible sequences
	target_primers = [primer]
	if "Y" in primer: # C or T
		target_primers = [primer.replace("Y", "C"), primer.replace("Y", "T")] 

	elif "R" in primer: # A or G
		target_primers = [primer.replace("R", "A"), primer.replace("R", "G")]

	# Searh each chromosome for the primer(s)
	for record in SeqIO.parse(reference_genome_path, "fasta"):
		for target_primer in target_primers:
			i = str(record.seq).find(target_primer)

			# If primer (target_primer) is found, record its genomic location and target sequence to the list
			# And, initiate another search in the same chromosome starting from that location (i+1)
			while i >=0:
				posb_targetseq = str(record.seq)[i:i + seq_len]
				possible_sequences.append((record.id, i, posb_targetseq))

				i = str(record.seq).find(target_primer, i + 1)

	return possible_sequences


def extract_amplicon_sequence(primers_df, found_sequences):

	# Create empty list to store all identified amplicons
	found_amplicons = []

	# Search each possible target sequence for non-self complement primers to identify amplicons
	for index, sequences in found_sequences.items():
		# retrieve the forward primer used to find the possible target sequence
		forward_primer = primers_df[primers_df.number == index].primer.values[0]

		# Identify last occurence of non-self complement primers in the possible target sequence
		for record_id, position, forw_sequence in sequences:
			for idx, primer, complement in primers_df[~(primers_df.number == index)].values:
				n = forw_sequence.rfind(complement, len(forward_primer)) # Start the search after forward_primer

				if n >= 0:
					found_amplicon  = forw_sequence[0:n+len(complement)] # Retrieve the amplicon sequence
					target_sequence = found_amplicon[len(forward_primer):len(found_amplicon)-len(complement)] # Strip the amplicon sequence for target sequence
					
					# Append the found amplicons to the list
					found_amplicons.append((record_id, position, index, forward_primer, complement, found_amplicon, len(found_amplicon), target_sequence))

	return found_amplicons


def generate_amplicon_records(found_amplicons, output_path):
	# Put the found amplicons into dataframe & save to file
	all_amplicons = pd.DataFrame(found_amplicons,
		columns = ["chromosome", "start", "primer_number", "forward_primer", "reverse_primer", "amplicon_sequence", "amplicon_length", "target_sequence"])

	# Save to file
	all_amplicons.to_csv(f"{output_path}/all_amplicons.csv")

	return all_amplicons


def generate_targetSeq_fasta(all_amplicons, output_path):
	# Convert the found sequences into fasta format

	# create list to store SeqRecord objects for fasta file
	sequences = []

	# Create sequence records of all found amplicons to write them into fasta file:
	for index, row in all_amplicons.iterrows():
		record = SeqRecord(
			Seq(row['target_sequence']),
			id = f"{row['chromosome']}:{row['start']}",
			description = f"forward_primer: {row['forward_primer']} reverse_primer: {row['reverse_primer']}" 
		)
		sequences.append(record)

	# Write the target sequence record into fasta file
	SeqIO.write(sequences, f"{output_path}/target_sequences.fasta", "fasta")


if __name__ == '__main__':

	parser = argparse.ArgumentParser(description = 'Retrieve target sequences')
	parser.add_argument('-r','--reference_genome', required=True, help='path to the reference genome fasta file location')
	parser.add_argument('-p','--primers', help='path to the list of primers', default="Data/Primers.txt")
	parser.add_argument('-l', '--seq_length', help='maximum length for amplicon', default=1000)
	parser.add_argument('-o', '--output', help='path to the output folder', default="Output")
	args = parser.parse_args()
	
	main(args)