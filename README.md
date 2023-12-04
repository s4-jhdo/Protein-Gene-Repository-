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
```
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
The commands above was used to print the AKAP alignment into a file for later analysis within the research paper. 

# Reconciling the AKAP Gene Family with the Species Tree
In order to continue and stay organized with bioinformatic data collection, create a directory by using the following command: 
```bash
cd ~/labs/lab06-$MYGIT/AKAP
```
The command creates a new directory to work under for proper file storage.

For the first step, a midpoint-rooted gene tree that was generated previously for the AKAP gene family. You can check if the tree is present with the following command: 
```bash
ls ~/labs/lab06-$MYGIT/AKAP/AKAP.homologs.al.mid.treefile
```
It is important to note that a midpoint rooted gene tree will be used for your own protein and gene family. To make sure you are utilizing your own protein use the following command line: 
```bash
mkdir ~/labs/lab06-$MYGIT/AKAP
```
These two commands above help view if your midpoint-rooted gene was generated correctly and to be sure that you are working under the correct directory. 

Then we would want to make a copy of this gene tree from the previous stages (in this case lab5, as it was labeled as such) into your current directory folder (in this case lab6) with the following command line: 
```bash
cp ~/labs/lab05-$MYGIT/AKAP/AKAP.homologs.al.mid.treefile ~/labs/lab06-$MYGIT/AKAP/AKAP.homologs.al.mid.treefile
```
# Reconcile the Gene and Species Tree Utilizing Notung
First, to perform the reconciliation use the following command line: 
```bash
java -jar ~/tools/Notung-3.0-beta/Notung-3.0-beta.jar -s ~/labs/lab06-$MYGIT/species.tre -g ~/labs/lab06-$MYGIT/AKAP/AKAP.homologs.al.mid.treefile --reconcile --speciestag prefix --savepng --events --outputdir ~/labs/lab06-$MYGIT/AKAP/
```
This command will generate the reconciliation tree of your own specific gene family. 

Second, view the "duplications" and "losses" columns for your specific tree with the following command: 
```bash
nw_display ~/labs/lab06-$MYGIT/species.tre
```
This command line helps focus specifically on the "Duplications" and "Losses" columns. 

Then, view the name that Notung assigned to the internal nodes that do not have formal taxonomic names with the following command:
```bash
grep NOTUNG-SPECIES-TREE ~/labs/lab06-$MYGIT/AKAP/AKAP.homologs.al.mid.treefile.reconciled | sed -e "s/^\[&&NOTUNG-SPECIES-TREE//" -e "s/\]/;/" | nw_display -
```
This command is used as a tool in order to know which lineages the internal nodes came from that did not have a formal taxonomic name. 

Finally, generate a RecPhyloXML object and view the gene-within-species tree with thirdkind by using the following commands: 
```bash
python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g ~/labs/lab06-$MYGIT/AKAP/AKAP.homologs.al.mid.treefile.reconciled --include.species
```
```bash
thirdkind -Iie -D 40 -f ~/labs/lab06-$MYGIT/AKAP/AKAP.homologs.al.mid.treefile.reconciled.xml -o  ~/labs/lab06-$MYGIT/AKAP/AKAP.homologs.al.mid.treefile.reconciled.svg
```
These commands were used to create a gene-reconciliation-with species tree reconciliation with the use of thirdkind to view the gene tree. 

# Protein Domain Prediction 
First, we need to use unaligned protein sequences from the previous bioinformatic steps. In order to do this, make a directory for the AKAP sequences and change into that specific directory with the following command line: 
```bash
mkdir ~/labs/lab08-$MYGIT/AKAP && cd ~/labs/lab08-$MYGIT/AKAP
```
```bash
cp ~/labs/lab05-$MYGIT/AKAP/AKAP.homologs.al.mid.treefile ~/labs/lab08-$MYGIT/AKAP
```
Second, make a copy of the raw unaligned sequence by using sed's substitute command and direct it into our specific folder (in this case, lab8) with the following command: 
```bash
sed 's/*//' ~/labs/lab04-$MYGIT/AKAP/AKAP.homologs.fas > ~/labs/lab08-$MYGIT/AKAP/AKAP.homologs.fas
```
Sed's substiture command is used to substitute any instance of an asterisk with nothing. The other command is used to output our AKAP protein folder into lab8 for later analysis. 

Third, download the Pfam database with the following command: 
```bash
wget -O ~/data/Pfam_LE.tar.gz ftp://ftp.ncbi.nih.gov/pub/mmdb/cdd/little_endian/Pfam_LE.tar.gz && tar xfvz ~/data/Pfam_LE.tar.gz  -C ~/data
```
This command is used to download the Pfam database which contains models that enable us to predict domains and important sites alongside functional analysis of proteins through classification. 

# Running RPS-Blast
First, to run RPS-Blast, use the following command: 
```bash
rpsblast -query ~/labs/lab08-$MYGIT/AKAP/AKAP.homologs.fas -db ~/data/Pfam -out ~/labs/lab08-$MYGIT/AKAP/AKAP.rps-blast.out  -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue .0000000001
```
This command line will give us the output of the Pfam data for the AKAP gene which will be used in later analysis and research paper. 

Second, run a script that enables the plotting of the pfam domain predictions from rps-blast next to their cognate protein on the phylogeny with the following command: 
```bash
sudo /usr/local/bin/Rscript  --vanilla ~/labs/lab08-$MYGIT/plotTreeAndDomains.r ~/labs/lab08-$MYGIT/AKAP/AKAP.homologs.al.mid.treefile
~/labs/lab08-$MYGIT/AKAP/AKAP.rps-blast.out
~/labs/lab08-$MYGIT/AKAP/AKAP.tree.rps.pdf
```
The Sudo command is used to run as a computer administrator in the case of needing to install packages. /usr/local/bin/Rscript is used as a program to allow you to run an R script from the command line without opening up the R console. --vanilla command indicates to R to not save or restore a workspace or previous settings. The tree Rps-blast output file will generate the midpoint rooted tree file for the AKAPS proteins, pfam domains from rps-blast and an output pdf file will be created. 

Now, we would want to look at the predited domains in more detail for our research paper. 
To focus on the predictions in the Pfam database use the following command line: 
```bash
mlr --inidx --ifs "\t" --opprint  cat ~/labs/lab08-$MYGIT/AKAP/AKAP.rps-blast.out | tail -n +2 | less -S
```
This command line will generate a domain-on-tree graphic.

Then, examine the distribution of Pfam domains across proteins with the following command line: 
```bash
cut -f 1 ~/labs/lab08-$MYGIT/AKAP/AKAP.rps-blast.out | sort | uniq -c
```
```bash
cut -f 6 ~/labs/lab08-$MYGIT/AKAP/AKAP.rps-blast.out | sort | uniq -c
```
```bash
awk '{a=$4-$3;print $1,'\t',a;}' ~/labs/lab08-$MYGIT/AKAP/AKAP.rps-blast.out |  sort  -k2nr
```
```bash
awk '{a=$4-$3;print $1,'\t',a;}' ~/labs/lab08-$MYGIT/AKAP/AKAP.rps-blast.out |  sort  -k2nr
```
The commands above will generate and tell us which proteins have more than one annotation, which Pfam domain annotation is commonly found, which protein has the longest annotated protein domain, and which protein has the shortest annotated protein domains respectively. 

To sort which protein has a domain annotation with the best e-value use the following command: 
```bash
cut -f 1,5 -d $'\t' ~/labs/lab08-$MYGIT/AKAP/AKAP.rps-blast.out | sort -k2,2rn -t $'\t'
```
This command is used to know which protein has a domain annotation with the best e-value while avoiding the process by hand. 



