README
---
To run this code, run the main file with Python along with the file path to the reference genome and paired reads file. The following is the command for this file hierarchy:

python .\Project1A.py .\project1a_10000_reference_genome.fasta .\project1a_10000_with_error_paired_reads.fasta

To change the genome or paired reads file used, simply change the file paths. 

NOTE: With the 10000 length reference genome and a kmer size of 10, the program takes around 15 minutes to run. Progress can be checked via the console output of the program, with the current read number being printed out along with the time elapsed. 

Results will be printed out to the console and formatted in mutations.txt. The results file path and file name can be easily modified within the program. 