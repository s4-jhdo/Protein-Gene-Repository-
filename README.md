# Protein-Gene-Repository
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

# Gene Family Sequence Alignment 
First, we need to install some software and programs to turn text-based sequence alignment into an HTML file (with aha) and open-source package/environment management system (with conda). Then we will install a2ps to turn the HTML file into a postscript file and install yum as a default software manager. (This installation only needs to happen ONCE by typing in these commands). 
```bash
conda install -y -n base -c conda-forge aha
```
```bash
sudo yum install -y a2ps
```
```bash
pip install buddysuite
```
Then, we will install a VScode extension to view PDFS by installing (vscode-pdf by tomoki1207). 

# Aligning AKAP Gene Family Members 
Before inputting the following commands, create a folder to organize your specific protein for future use and access to bioinformatic tools. (Make sure you are in the correct folder by changing the directory if needed). Use the following commands: 
```bash
mkdir ~/labs/lab04-$MYGIT/AKAP
```
This command is used to create a folder for the AKAP gene family to be later used in a series of bioinformatic steps. 
```bash
cd ~/labs/lab04-$MYGIT/AKAP

Type in the following command to obtain the sequences that are within the BLAST output file:
```bash
seqkit grep --pattern-file ~/labs/lab03-$MYGIT/AKAP/AKAP.blastp.detail.filtered.out ~/labs/lab03-$MYGIT/allprotein.fas > ~/labs/lab04-$MYGIT/AKAP/AKAP.homologs.fas
```
Since we already downloaded the proteomes, Seqkit enables us to access them easily with the command above. 

# Performing a Global Multiple Sequence Alignment in MUSCLE: 

Utilize the following command to create a multiple alignment by using MUSCLE: 
```bash
muscle -in ~/labs/lab04-$MYGIT/AKAP/AKAP.homologs.fas -out ~/labs/lab04-$MYGIT/AKAP/AKAP.homologs.al.fas
```
Next, view the alignment in ALV with the following command: 
```bash
alv -kli  ~/labs/lab04-$MYGIT/AKAP/AKAP.homologs.al.fas | less -RS
```
We can also try the "majority" option to observe some regions where the alignment is believable or not by typing in the command: 
```bash
alv -kli --majority ~/labs/lab04-$MYGIT/AKAP/AKAP.homologs.al.fas | less -RS
```
Then, we will create this alignment into a PDF format for usage within our research investigation by typing in the following commands: 
```bash
muscle -in ~/labs/lab04-$MYGIT/AKAP/AKAP.homologs.fas -html -out ~/labs/lab04-$MYGIT/AKAP/AKAP.homologs.al.html
```
```bash
alv -ki -w 100 ~/labs/lab04-$MYGIT/AKAP/AKAP.homologs.al.fas | aha > ~/labs/lab04-$MYGIT/AKAP/AKAP.homologs.al.html
```
```bash
a2ps -r --columns=1 ~/labs/lab04-$MYGIT/AKAP/AKAP.homologs.al.html -o ~/labs/lab04-$MYGIT/AKAP/AKAP.homologs.al.ps
```
```bash
ps2pdf ~/labs/lab04-$MYGIT/AKAP/AKAP.homologs.al.ps ~/labs/lab04-$MYGIT/AKAP/AKAP.homologs.al.pdf
```
# Obtaining Information About the Alignment 
To obtain and calculate the length of the alignment use the following command: 
```bash
alignbuddy  -al  ~/labs/lab04-$MYGIT/AKAP/AKAP.homologs.al.fas
```
We can also calculate the length of the alignment after removing any columns with gaps with the following command: 
```bash
alignbuddy -trm all  ~/labs/lab04-$MYGIT/AKAP/AKAP.homologs.al.fas | alignbuddy  -al
```
The command calculates the alignment length with gaps for a more substantial alignment for the specific gene in the study. 

We can also calculate the length of the alignment after removing completely conserved positions with the following command: 
```bash
alignbuddy -dinv 'ambig' ~/labs/lab04-$MYGIT/globins/globins.homologs.al.fas | alignbuddy  -al
```
# Calculating the Average Percent Identity 
To collect the average percent identity, use the following t_coffee command: 
```bash
t_coffee -other_pg seq_reformat -in ~/labs/lab04-$MYGIT/AKAP/AKAP.homologs.al.fas -output sim
```
This command is utilized to calculate the average percent identity among all sequences in the alignment. 

Then, repeat calculating the average percent identity with alignbuddy by using the following command: 
```bash
 alignbuddy -pi ~/labs/lab04-$MYGIT/AKAP/AKAP.homologs.al.fas | awk ' (NR>2)  { for (i=2;i<=NF  ;i++){ sum+=$i;num++} }
     END{ print(100*sum/num) } 
```
# Constructing a Phylogenetic Tree for AKAP Homologs from Sequence Data
First, the software IQ-TREE will be employed to deduce the most suitable phylogenetic tree from a sequence alignment. To organize this part's bioinformatic data collection and session, first create a directory for this part and label accordingly in sequence to the previous parts using the command: 
```bash
mkdir ~/labs/lab05-$MYGIT/AKAP
cd ~/labs/lab05-$MYGIT/AKAP 
```
Mkdir is used to make a new directory, and CD is used here to change to that directory. 

In order to continue, copy the alignment from AKAP folder from lab04 into lab05 as IQTree will put the output files in the same directory as the input file. Use the following command to go to that specific directory: 
```bash
cp ~/labs/lab04-$MYGIT/AKAP/AKAP.homologs.al.fas ~/labs/lab05-$MYGIT/AKAP/AKAP.homologs.al.fas   
```
The command above is utilized to go to a specific directory that IQTree sorted in terms of files. 

Second, we will use IQ-Tree to find the maximum likehood tree estimate by utilizing the following command: 
```bash
iqtree -s ~/labs/lab05-$MYGIT/AKAP/AKAP.homologs.al.fas -bb 1000 -nt 2 
```
Initially, the software will compute the optimal amino acid substitution model and amino acid frequencies. Subsequently, it will conduct a tree search, estimating branch lengths during the process with the command above. 

# Rooting the Optimal Phylogeny
First, to root the optimal phylogeny we will be using a type of rooting named midpoint by using the following command: 
```bash
gotree reroot midpoint -i ~/labs/lab05-$MYGIT/AKAP/AKAP.homologs.al.fas.treefile -o ~/labs/lab05-$MYGIT/AKAP/AKAP.homologs.al.mid.treefile
```
The command above will provide a view of an optimal phylogenetic tree that can be used for our research. 

Second, we can view the rooted tree with the command line: 
```bash
nw_order -c n ~/labs/lab05-$MYGIT/AKAP/AKAP.homologs.al.mid.treefile  | nw_display -
```
This command is used to view the resulting ASCII Tree that can be used to support our investigation. 

Then, to make the large trees easier to view, nw_order can be used as a graphic with the command line: 
```bash
nw_order -c n ~/labs/lab05-$MYGIT/AKAP/AKAP.homologs.al.mid.treefile | nw_display -w 1000 -b 'opacity:0' -s  >  ~/labs/lab05-$MYGIT/AKAP/AKAP.homologs.al.mid.treefile.svg -
```
The nw_order initially arranges the clades by the number of descendants, a process known as ladderization, aiming to enhance the visibility of large trees. However, it's important to note that this approach may introduce certain biases in interpretation.

# Branch Lengths 
To distinguish between phylograms and cladograms, due to the need of visualizing short branch lengths use the following command line: 
```bash
nw_order -c n ~/labs/lab05-$MYGIT/AKAP/AKAP.homologs.al.mid.treefile | nw_topology - | nw_display -s  -w 1000 > ~/labs/lab05-$MYGIT/AKAP/AKAP.homologs.al.midCl.treefile.svg -
```
The default tree displayed in nw_display (and many other programs) is a phylogram, indicating that the lengths of branches correspond proportionally to the accumulated number of substitutions in the sequence along each branch. When dealing with very short branch lengths, visualizing clades on a phylogram can be challenging. Consider switching the view to a cladogram using the following command.

# Outgroup Rooting 
A secondary approach to rooting is to specify an outgroup. To express an outgroup use the following command lines: 
```bash
nw_reroot subphyla.tre Echinodermata >subphyla.echinoderm_reroot.tre
nw_display ~/labs/lab05-$MYGIT/subphyla.echinoderm_reroot.tre
```
This command allows us to reroot the tree with the branching of echinoderms as the initial event. It's important to note that although this establishes a root, the true root is the common ancestor of Echinoderms/Hemichordates and the remaining deuterostomes.

# Printing AKAP Alignment into PDF for Research 
First, we need to install latex suit and use Rscript to print the alignments faster with the following command: 
```bash
sudo yum install texlive
```
```bash
Rscript --vanilla ~/labs/lab05-$MYGIT/installMSA.R
```
```bash
Rscript --vanilla ~/labs/lab05-$MYGIT/plotMSA.R  ~/labs/lab05-$MYGIT/AKAP/AKAP.homologs.al.fas ~/labs/lab05-$MYGIT/AKAP/AKAP.homologs.al.fas.pdf
```
