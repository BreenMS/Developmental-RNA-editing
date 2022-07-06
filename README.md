# Spatiotemporal and genetic regulation of RNA editing across human brain development
Winston H. Cuddleston*, Xuanjia Fan*, Laura Sloofman*, Lindsay Liang*, Enrico Mossotto, Kendall Moore, Sarah Zipkowitz, Minghui Wang, Bin Zhang, Jiebiao Wang, Nenad Sestan, Bernie Devlin, Kathryn Roeder, Stephan J. Sanders, Joseph D. Buxbaum, Michael S. Breen6**

*Co-first authors: WHC, XF, LS, LL
**Correspondence to: michael.breen@mssm.edu, @breenPsychgene

Posttranscriptional RNA modifications by adenosine-to-inosine (A-to-I) editing are abundant in the brain, yet elucidating functional sites remains challenging. To bridge this gap, we expose spatiotemporal and genetically regulated A-to-I editing sites across prenatal and postnatal stages of human brain development. We identify more than 10,000 temporally regulated and spatially conserved editing sites that occur predominately in 3’UTRs and introns, as well as 37 sites that recode amino acids in protein coding regions with precise changes in editing levels across development. A massive expansion of RNA hyper-edited transcripts accumulates in the brain well into advanced age, stabilizing RNA secondary structures. These features are conserved in murine and non-human primate models of neurodevelopment. Further, we uncover thousands of cis-editing quantitative trait loci (edQTLs) with unique regulatory effects during prenatal and postnatal development. Collectively, this work offers a highly resolved atlas linking spatiotemporal variation in editing levels to genetic regulatory effects throughout distinct stages of brain maturation.  <br /> <br /> 

Here we provide all computational code, supplemental tables, and intermediate data files used to generate all results and figures in this body of work. 

# This work entails four main levels of analysis:
1. Compute an Alu Editing Index (AEI) from a STAR mapped bam file  [(RNAEditingIndexer v1.0)](https://github.com/a2iEditing/RNAEditingIndexer)<br /> 
2. Quantifying RNA editing sites from STAR mapped bam files using de novo methods [(reditools v2.0)](https://github.com/tizianoflati/reditools2.0)<br /> 
3. Quantifying RNA editing from STAR mapped bam files using a list of predefined list of sites (code provided below)<br /> 
4. Quantifying RNA hyper-editing sites from STAR unmapped fastq files [(method based on Porath et al., 2017)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1315-y)<br /> 

# 1. Compute AEI from a STAR mapped bam file:
We used already available software from the [RNAEditingIndexer GitHub account](https://github.com/a2iEditing/RNAEditingIndexer) to compute an AEI based on a mapped bam file. The method is describe in the original publication, [Nat. Methods (2019)](https://pubmed.ncbi.nlm.nih.gov/31636457/). Here, we provide an example of bash shell script that executes the AEI on one sample. Requirements and parameters are described in full in the bash script.  <br /> 
 
An example for computing AEI on human samples:
```ruby
RNAEditingIndex -d -f Aligned.sortedByCoord.out.bam -o .
--genes_expression ucscHg38GTExGeneExpression.bed.gz
--refseq ucscHg38RefSeqCurated.bed.gz
--snps ucscHg38CommonGenomicSNPs150.bed.gz
-gf ucscHg38Genome.fa
-rb ucscHg38Alu.bed.gz
--genome UserProvided  --paired_end --stranded
```
<br />  

An example for computing AEI on macaque samples:
```ruby
RNAEditingIndex -d -f Aligned.sortedByCoord.out.bam -o .
--genes_expression RheMac8_selectedregions_expression_avg.bed.gz
--refseq RheMac8_refSeq.bed.gz
--snps rheMac8_common_0.01_SNVs.bed.gz
-gf rheMac8.fa
-rb RheMac8_SINE.bed.gz
--genome UserProvided  
```
<br />  

An example for computing AEI on mouse samples:
```ruby
RNAEditingIndex -d -f Aligned.sortedByCoord.out.bam -o .
--genes_expression ucscMM10GTExGeneExpression.bed.gz
--refseq ucscMM10RefSeqCurated.bed.gz
--snps ucscMM10CommonGenomicSNPs142.bed.gz
-gf ucscMm10Genome.fa
-rb ucscMM10SINE_B1_B2.bed.gz
--genome UserProvided  --paired_end 
```
<br />  


# 2. Quantify RNA editing sites from STAR mapped bam files using de novo methods:
We used already available software from the [reditools v2.0 GitHub account](https://github.com/tizianoflati/reditools2.0) and [JACUSA2 GitHub account](https://github.com/dieterich-lab/JACUSA2) to quantify de novo RNA editing sites based on a STAR mapped bam file. The methods are describe in the original publications: [BMC Bioinformatics (2020)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03562-x) and [Genome Biology (2022)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02676-0). Here, we provide examples of how we executed reditools 2.0 and JACUSA2 on the BrainVar data set. 

<br /> 
REDITOOLS 2.0
```ruby
mpirun parallel_reditools.py -f Aligned.sortedByCoord.out.bam -r GRCh38.chrom.fa -S -s 2 -ss 5 -mrl 50 -q 10 -bq 20 -C -T 2 -m homopolymeric_sites_hg38.txt -os 5 -Z GRCh38.chrom.fa.fai -G Aligned.sortedByCoord.out.bam.cov -D Aligned.sortedByCoord.out.bam_out

```
<br />  

JACUSA2
```ruby
java -jar $JACUSA2_JAR call-1 -r Aligned.sortedByCoord.out.bam -p 10 -a D,M,Y,E:file=hg38-blacklist.v2_sort.bed:type=BED -s -m 20 -R GRCh38.chrom.fa -P RF-FIRSTSTRAND
```
<br />  



# 3. Quantify RNA editing from STAR mapped bam files using a list of predefined list of sites (based on a predefined list of sites):
It is often of interest to quantify RNA editing sites based on a user defined list of sites. Samtools mpileup has the functionality to execute this task. Here we provide two perl scripts that will achieve this task. The only requirement is installing a recent version of samtools.  <br /> 


query_known_sites.pl= excute mpileup (samtools) to query a list of known editing sites.<br />
parse_pileup_query.pl = a requirement for query_known_sites.pl<br />  
Usage: perl query_known_sites.pl [A predefined list of known editing sites] [STAR mapped bam file] [Output file name]
```ruby
perl query_known_sites.pl CNS_A2G_events.txt SampleName.bam OutputFileName.txt
```
<br />  

# Helpful data files:
CNS_A2G_events.txt = A predefined list of 166,215 A-to-I RNA editing sites detected within each cell population.<br /> 
CNS_A2G_15221edits.txt = A matrix of 15,221 RNA editing sites we detected across all three cell populations.<br /> 

<br />  

# 4. Quantifying RNA hyper-editing sites from STAR unmapped fastq files:
We used already available software from the [RNA hyper-editing GitHub account](https://github.com/hagitpt/Hyper-editing) to quantify RNA hyper-editing sites based on a STAR unmapped fastq files. The method is describe in the original publication, [Genome Biology (2017)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1315-y). We provide a detailed example of how we execute this pipeline providing all code, with minor modifications, at [Ryn Cuddleston's GitHub account](https://github.com/ryncuddleston/RNA-hyper-editing) .  <br /> 


<br />  

# Supplemental Data Tables 1-9:
Table S1. Alu editing index (AEI) and metadata across all samples used in the current study. <br /> 
Table S2. Annotation of all prenatal, postnatal and commonly edited sites.<br /> 
Table S3. Developmentally regulated RNA editing sites and annotations. <br /> 
Table S4. Developmentally regulated RNA editing sites in independent datasets. <br /> 
Table S5. RNA hyper-editing results across all datasets. <br /> 
Table S6. RNA hyper-editing, splicing , and intron retention results. <br /> 
Table S7. AEI, hyper-editing, and metadata for animal models of neurodevelopment. <br /> 
Table S8. Top temporal edQTLs. <br /> 
Table S9. GWAS-edQTL co-localization with edQTLs.  <br /> 


<br />  


# Additional data files:
