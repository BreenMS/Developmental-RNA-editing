# Spatiotemporal and genetic regulation of RNA editing across human brain development
Winston H. Cuddleston*, Xuanjia Fan*, Laura Sloofman*, Lindsay Liang*, Enrico Mossotto, Kendall Moore, Sarah Zipkowitz, Minghui Wang, Bin Zhang, Jiebiao Wang, Nenad Sestan, Bernie Devlin, Kathryn Roeder, Stephan J. Sanders, Joseph D. Buxbaum, Michael S. Breen**

*Co-first authors: WHC, XF, LS, LL<br /> 
**Correspondence to: michael[dot] breen [at] mssm [dot] edu, @breenPsychgene

Posttranscriptional RNA modifications by adenosine-to-inosine (A-to-I) editing are abundant in the brain, yet elucidating functional sites remains challenging. To bridge this gap, we expose spatiotemporal and genetically regulated A-to-I editing sites across prenatal and postnatal stages of human brain development. We identify more than 10,000 temporally regulated and spatially conserved editing sites that occur predominately in 3’UTRs and introns, as well as 37 sites that recode amino acids in protein coding regions with precise changes in editing levels across development. A massive expansion of RNA hyper-edited transcripts accumulates in the brain well into advanced age, stabilizing RNA secondary structures. These features are conserved in murine and non-human primate models of neurodevelopment. Further, we uncover thousands of cis-editing quantitative trait loci (edQTLs) with unique regulatory effects during prenatal and postnatal development. Collectively, this work offers a highly resolved atlas linking spatiotemporal variation in editing levels to genetic regulatory effects throughout distinct stages of brain maturation.  <br /> <br /> 

![Abstract (1)](https://user-images.githubusercontent.com/22500312/177622931-7c53720c-11cb-4787-b266-5969afbb6c71.png)

Here we described the main computational code used to generate all results and figures in this body of work.  

All RNA editing matrices and sample level RNA editing sites across the DLPFC, forebrain and hindbrain can be downloaded from [Synapse.org](https://www.synapse.org/#!Synapse:syn26434508/files/).  

# This work entails four main levels of analysis:
1. Compute an Alu Editing Index (AEI) from a STAR mapped bam file  [(RNAEditingIndexer v1.0)](https://github.com/a2iEditing/RNAEditingIndexer)<br /> 
2. Quantifying RNA editing sites from STAR mapped bam files using de novo methods [(reditools v2.0)](https://github.com/tizianoflati/reditools2.0) and [(JACUSA2)](https://github.com/dieterich-lab/JACUSA2)<br /> 
3. Quantifying RNA editing from STAR mapped bam files using a list of predefined list of sites (code provided below)<br /> 
4. Quantifying RNA hyper-editing sites from STAR unmapped fastq files [(method based on Porath et al., 2017)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1315-y)<br /> 

An overview of our analytical approach:<br /> 
<img width="500" alt="Screen Shot 2022-07-07 at 10 16 40 AM" src="https://user-images.githubusercontent.com/22500312/177795799-90edf08a-106e-44a6-a646-4569efdac555.png">

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
We used already available software from the [reditools v2.0 GitHub account](https://github.com/tizianoflati/reditools2.0) and [JACUSA2 GitHub account](https://github.com/dieterich-lab/JACUSA2) to quantify de novo RNA editing sites based on a STAR mapped bam file. The methods are describe in the original publications: [BMC Bioinformatics (2020)](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-020-03562-x) and [Genome Biology (2022)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02676-0). Here, we provide examples of how we executed reditools 2.0 and JACUSA2 on the BrainVar data set. <br /> 

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
It is often of interest to quantify RNA editing sites based on a user defined list of sites. Samtools mpileup has the functionality to execute this task. Here we provide two perl scripts that will achieve this task. The only requirement is installing a recent version of samtools. In the current study, we leveraged lists of known sites from these three resources: [REDIportal](https://academic.oup.com/nar/article/49/D1/D1012/5940507), [cellular and genetic drivers of RNA editing variability in the human brain](https://www.nature.com/articles/s41467-022-30531-0), and [an atlas of human recoding sites](https://www.nature.com/articles/s41467-022-28841-4). <br /> 


query_known_sites.pl= excute mpileup (samtools) to query a list of known editing sites.<br />
parse_pileup_query.pl = a requirement for query_known_sites.pl<br />  
Usage: perl query_known_sites.pl [A predefined list of known editing sites] [STAR mapped bam file] [Output file name]
```ruby
perl query_known_sites.pl CNS_A2G_events.txt SampleName.bam OutputFileName.txt
```
<br />  

# Helpful data files:
CNS_A2G_events.txt = A predefined list of 166,215 A-to-I RNA editing sites detected within each cell population can be found [here](https://github.com/BreenMS/RNA-editing-in-CNS-cell-types).<br /> 

<br />  

# 4. Quantifying RNA hyper-editing sites from STAR unmapped fastq files:
We used already available software from the [RNA hyper-editing GitHub account](https://github.com/hagitpt/Hyper-editing) to quantify RNA hyper-editing sites based on a STAR unmapped fastq files. The method is describe in the original publication, [Genome Biology (2017)](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1315-y). We provide a detailed example of how we execute this pipeline [here](https://github.com/ryncuddleston/RNA-hyper-editing) .  <br /> 


<br />  

# Supplemental Data Tables 1-10:
Table S1. Alu editing index across all samples in the current study. Alu editing indexes and corresponding metadata for all samples in the current study. (Code described above) <br />  
Table S2. Sample level summary of all high-quality RNA editing sites. Summary of all high-quality RNA editing sites detected and their genomic locations per sample in the DLPFC, forebrain and hindbrain. (All data available on [Synapse](https://www.synapse.org/#!Synapse:syn26434508/files/)) <br />  
Table S3. Differentially regulated RNA editing sites across development. Summary statistics for sites that are dynamically regulated across brain development in the DLPFC, forebrain and hindbrain. We also provide the following information: 1) pathway enrichment; 2) disease-gene set enrichment; 3) associations between changes in editing levels and gene expression levels; 4) validation rates in hiPSCs and aging samples.  <br />  
<br />  An example of running differential editing:

```ruby
library(limma)
library(edgeR)

GeneExprs <- read.delim("BrainVar_editing_matrix.txt", check.names=FALSE, stringsAsFactors=FALSE, row.names=1)
info <- read.delim("BrainVar_metadata.txt", check.names=FALSE, stringsAsFactors=FALSE) #Count_NPC_noOut_targets.txt

Stage = as.factor(info$Stage)
Ancestry = as.factor(info$Ancestry)
Sex = as.factor(info$Sex)
ADAR1 = (info$ADAR1)
ADAR2 = (info$ADAR2)
ADAR3 = (info$ADAR3)
RIN = (info$RIN)
Neurons = (info$Neurons)

design <- model.matrix(~0+Stage+Sex+Ancestry)
#design <- model.matrix(~0+Stage+Sex+Ancestry+Neurons) #model to correct for neuronal fractions
#design <- model.matrix(~0+Stage+Sex+Ancestry+Neurons+ADAR1+ADAR2) #model to correct for neuronal fractions, ADAR1, and ADAR2
fit <- lmFit(GeneExprs,design)
cm <-makeContrasts(DevEffect = (StagePostnatal - StagePrenatal),levels=design)
fit2 <- contrasts.fit(fit, cm)
fitDupCor <- eBayes(fit2)
topTable(fitDupCor, coef="DevEffect")
DE_sites<- topTable(fitDupCor, coef="DevEffect", n=nrow(GeneExprs))
write.table(DE_sites, "DEG_BrainVar_sites.txt", sep="\t")
```

Table S4. miRNA binding affinity predictions to 3’UTRs. miRanda analysis output for editing sites overlapping miRNA seed regions and the resulting alignments and minimum free energy calculations for edited and unedited sequences.<br />  
<br />  An example of running miRANDA:
```ruby
miranda hsa.mature.fa 3UTR_sites_targetA_unedited.fa -strict > miRanda_out_targetA_unedited.txt
miranda hsa.mature.fa 3UTR_sites_targetG_edited.fa -strict > miRanda_out_targetG_edited.txt

grep -A 1 "Scores for this hit:" miRanda_out_targetA_unedited.txt | sort | grep '>' > match_miRanda_out_targetA_unedited.txt
grep -A 1 "Scores for this hit:" miRanda_out_targetG_edited.txt | sort | grep '>' > match_miRanda_out_targetG_edited.txt

```


Table S5. Temporally regulated recoding sites. This table includes information temporally regulated recoding sites across human brain development, pathway level enrichment and validation in hiPSC models of corticogenesis.  <br />  
Table S6. Temporally regulated RNA hyper-editing events. Summary statistics on RNA hyper-editing events per gene across development as well as per sample for all transcriptome samples in the current study. <br />  
Table S7. Mechanistic investigation of RNA hyper-editing. Summary statistics for RNA hyper-editing sites that are predicted to be splice alternating via SpliceAI, occur in retained introns via SIRI.<br />  
<br />  An example of running SpliceAI:
```ruby
spliceai -I BrainVar_sites.vcf -O output.vcf -R GRCh38.chrom.fa -A grch38
```
<br />  An example of running SIRI:
```ruby
SIRI --gtf gencode.v30.primary_assembly.annotation.gtf --bam Aligned.sortedByCoord.out.bam --anchor 8 --length 100 --lib first --read P 

```

Table S8. RNA editing summary results in animal models of neurodevelopment. Summary statistics on the AEI and global RNA hyper-editing per sample in mouse and macaque models of neurodevelopment. (Code described above). <br />  
Table S9. Temporal predominate edQTLs. Summary statistics for leading editing-variant pairs (edQTLs), including constant, prenatal and postnatal predominate edQTLs. <br />  
Table S10. GWAS-edQTL co-localization. Summary statistics of all significant GWAS risk loci that colocalize with edQTLs in the current study. <br />  


