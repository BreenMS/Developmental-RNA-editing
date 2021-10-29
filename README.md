# Spatiotemporal and genetic regulation of RNA editing across human brain development
Ryn Cuddleston Laura Sloofman, Lindsay Liang, Enrico Mossotto, Xuanjia Fan, Kendall Moore, Sarah Zipkowitz, Minghui Wang, Bin Zhang, Jiebiao Wang, Nenad Sestan, Bernie Devlin, Kathryn Roeder, Joseph D. Buxbaum, Stephan J. Sanders, Michael S. Breen

Adenosine-to-inosine editing is abundant in the brain and required for neurodevelopment, yet the developmental and genetic regulation for most sites remains uncharted. We constructed a full temporal view of base-specific RNA editing in the developing cortex, from early progenitors through fully mature cells in the adult brain. RNA editing is highly developmentally regulated, with more than 10,000 selective editing sites showing temporal changes in editing levels, shifting between mid-fetal periods and infancy. A massive expansion of RNA hyper-editing sites accumulate in the cortex well into advanced age, occurring disproportionally in 3’UTRs of neurodevelopmental genes. These profiles are preserved in non-human primate and murine models of cortical development, illustrating an evolutionarily conserved selected function. Further, we uncover thousands of RNA editing quantitative trait loci (edQTLs) across development, including 6,394 edQTLs with unique regulatory effects during prenatal and postnatal periods, and nine that co-localize with risk loci for neuropsychiatric disorders. We provide a highly resolved atlas of RNA editing sites involved in brain development, linking developmental variation in editing to genetic regulatory effects and the regulation of complex traits. <br /> <br /> 
Here we provide all computational code, supplemental tables, and intermediate data files used to generate all results and figures in this body of work. 

# This work entails four main levels of analysis:
1. Compute an Alu Editing Index (AEI) from a STAR mapped bam file  [(RNAEditingIndexer v1.0)](https://github.com/a2iEditing/RNAEditingIndexer)<br /> 
2. Quantifying RNA editing sites from STAR mapped bam files using de novo methods [(reditools v2.0)](https://github.com/tizianoflati/reditools2.0)<br /> 
3. Quantifying RNA editing from STAR mapped bam files using a list of predefined list of sites (code provided below)<br /> 
4. Quantifying RNA hyper-editing sites from STAR unmapped fastq files [(method based on Porath et al., 2017)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1315-y)<br /> 

# 1. Compute AEI from a STAR mapped bam file:
We used already available software from the [RNAEditingIndexer GitHub account](https://github.com/a2iEditing/RNAEditingIndexer) to compute an AEI based on a mapped bam file. The method is describe in the original publication, [Nat. Methods (2019)](https://pubmed.ncbi.nlm.nih.gov/31636457/). Here, we provide an example of bash shell script that executes the AEI on one sample. Requirements and parameters are described in full in the bash script.  <br /> 
 
```ruby
AEI.sh
```
<br />  

# 2. Quantify RNA editing sites from STAR mapped bam files using de novo methods:
We used already available software from the [reditools v2.0 GitHub account](https://github.com/tizianoflati/reditools2.0) to quantify de novo RNA editing sites based on a STAR mapped bam file. The method is describe in the original publication, [BMC Bioinformatics (2020)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03562-x). Here, we provide an example of bash shell script that executes reditools 2.0 on one sample. Requirements and parameters are described in full in the bash script.  <br /> 

```ruby
reditools_caller.sh
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
