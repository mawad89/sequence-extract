# sequence-extract
--------------------------------------------------------------------------------
- This project aim to extract each gene and intergenic regions from annotated  -
-                        genomes to separated fasta files                      -
--------------------------------------------------------------------------------

The project consest of three main functions (gene_extract, intron_extract, and 
primer_pick_up) in addition to one assistance function (rc)   

                          **********



(gene_extract) Function 
#####
This function determine the start and end of each gene according to
GenBank annotation and extract each gene in an indvedual fasta file
#####

input files:
1. Sequence Fasta format
2. sequence genbank format

return:
1. separated genes in fasta format
2. Number of genes at the input file  

                          **********




(intron_extract) Function 
#####
This function determine the start and end of each gene according to
    GenBank annotation and extract the whole entergenic regions
                       in a single fasta file
#####

input files:
1. Sequence Fasta format
2. sequence genbank format

return:
1. whole entergenic regions in a single fasta file

                          **********


(rc) Function

#####
The function recieve DNA sequence and return sequence reverse complementary
#####

input:
DNA sequence 

return:
sequence reverse complementary

                          **********



(primer_pick_up) Function 

#####
The function Pick up sequences between forward and reverse primers from 
                              different files                                 
#####

Function name is (primer_pick_up)

inputs:
1.Forward_primer
2.Reverse_primer
3. numer_of_fasta_files: number of files the user need to test the primers on it

return:
sequences between forward and reverse primers in fasta format file
                          **********



Run Notes:
The user should rename the tested sequence files to numbers started from 1, and put this file in the same folder with these python script.   
