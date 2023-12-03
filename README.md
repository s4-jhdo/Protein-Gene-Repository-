# Protein-Gene-Repository-
Hi! This guide outlines a systematic procedure for employing the tools, commands, and programs I utilized to collect the required data for the Bio 312 research paper focusing on the evolutionary analysis of a gene assigned within the context of deuterostomes.

# Creating a Directory
Firstly, create a directory where we can store all computations related to AKAP17A (or the specific protein you are working with) by executing the following command: 
```bash
mkdir /home/ec2-user/AKAP
```
This command will create a directory within the ec2-user folder in Visual Studio code (VScode). You also have the flexibility to designate a different or specific location for placing this file.

# Downloading AKAP17A gene/protein from NCBI Blast database
First, let's acquire the CDC40 gene as our query sequence. Obtain the sequence by executing the subsequent command:
```bash
ncbi-acc-download -F fasta -m protein NP_005079.2
```

This command downloads the AKAP17A gene sequence from the NCBI database and converts it into a fasta file.

# Perform the Blast Search 

To conduct a BLAST search with our query protein (CDC40), utilize the following command.
```bash
blastp -db ../allprotein.fas -query NP_005079.2.fa -outfmt 0 -max_hsps 1 -out AKAP.blastp.typical.out
```
# Conduct a BLAST Search and Specify Tabular Output Format
The subsequent command will generate a more intricate output of the earlier analysis.
```bash
blastp -db ../allprotein.fas -query NP_005079.2.fa  -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle"  -max_hsps 1 -out AKAP.blastp.detail.out
```
# Refining and Filtering the BLAST results to identify potential homologs with high scores.
We will then refine the homologs gathered through the BLAST search by incorporating only those with high-scoring matches. In this instance, the command employed involves an e-value of 1e-35, ensuring the consideration of homologs closely matching the AKAP17A gene copy identified in H. sapiens.
```bash
awk '{if ($6< 1e-35)print $1 }' globins.blastp.detail.out > AKAP.blastp.detail.filtered.out
```
To determine the number of paralogs, which was used within Table 1 of my research study, utilize the following command:
```bash
grep -o -E "^[A-Z]\.[a-z]+" AKAP.blastp.detail.filtered.out  | sort | uniq -c
```
After filtering to include high-scoring matches from the BLAST search, this command tallies the quantity of paralogs identified in each species.
