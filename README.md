# Extracting PCR target sequences using primers

The repository has two main scripts that are dedicated to each assignment:  

   1. *"retrieve_target_sequences.py"* for extracting target sequences from the reference genome using the provided list of primers.
   2. *"generate_short_hand_annotation.py"* for converting the extracted target sequences into short-hand annotations using the provided motif list.
   
The Data folder contains the provided list of primers and motifs, and the Output folder contains:  

   1. *"all_amplicons.csv"* - the list of all found amplicons along with their location, primers, and sequence;
   2. *"target_sequences.fasta"* - the fasta format of all extracted target sequences; 
   3. *"short_hand_annotations.csv"* - the short-hand annotations of each target sequence; and
   4. *"found_seqs.pickle"* - pickle file of all possible amplicons identified.
   
#### NOTE: The reference genome used in this project is GRCh38.p13, which can be found [#Here](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_33/)
