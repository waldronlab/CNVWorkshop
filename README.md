![](https://github.com/waldronlab/CNVWorkshop/workflows/.github/workflows/basic_checks.yaml/badge.svg)

# Copy number variation analysis with Bioconductor

[Workshop website](https://waldronlab.github.io/CNVWorkshop)

[Docker image](https://hub.docker.com/repository/docker/ludwigg/cnvworkshop)

# Instructor(s) name(s) and contact information

Ludwig Geistlinger, Marcel Ramos, Sehyun Oh, and Levi Waldron

CUNY School of Public Health
55 W 125th St, New York, NY 10027

Ludwig_Geistlinger@hms.harvard.edu  
Marcel.Ramos@sph.cuny.edu   
Sehyun.Oh@sph.cuny.edu   
Levi.Waldron@sph.cuny.edu   

# Workshop Description

This workshop gives an overview of Bioconductor solutions for the analysis of 
copy number variation (CNV) data. 
The workshop introduces Bioconductor core data structures for efficient 
representation, access, and manipulation of CNV data, and how to use these
containers for structured downstream analysis of CNVs and integration with gene 
expression and quantitative phenotypes. 
Participants will be provided with code and hands-on practice for a comprehensive 
set of typical analysis steps including exploratory analysis, summarizing individual 
CNV calls across a population, overlap analysis with functional genomic regions 
and regulatory elements, expression quantitative trait loci (eQTL) analysis, 
and genome-wide association analysis (GWAS) with quantitative phenotypes.
As an advanced application example, the workshop also introduces allele-specific 
absolute copy number analysis and how it is incorporated in cancer genomic analysis 
for the estimation of tumor characteristics such as tumor purity and ploidy. 

# Pre-requisites

* Basic knowledge of R syntax
* Familiarity with the SummarizedExperiment class
* Familiarity with the GenomicRanges class

* Familiarity with high-throughput genomic assays such as microarrays and 
  next-generation sequencing
* Familiarity with the biological definition of single nucleotide polymorphism 
  (SNP) and copy number variation (CNV) 

# Workshop Participation

Execution of example code and hands-on practice

# _R_ / _Bioconductor_ packages used

* [RaggedExperiment](http://bioconductor.org/packages/RaggedExperiment)
* [CNVRanger](http://bioconductor.org/packages/CNVRanger)
* [regioneR](http://bioconductor.org/packages/regioneR)
* [PureCN](http://bioconductor.org/packages/PureCN)

# Time outline

| Activity                                              | Time |
|-------------------------------------------------------|------|
| Overview                                              | 5m   |
| Data representation and manipulation                  | 20m  |
| Integrative downstream analysis (eQTL, GWAS, ...)     | 20m  |
| Allele-specific CN analysis in cancer                 | 15m  |


# Learning Goals

* Gain familiarity with elementary concepts of CNV analysis
* Learn how to efficiently represent, access, and manipulate CNV data 
  in Bioconductor data structures
* Gain familiarity with different strategies for summarizing individual CNV
  calls across a population
* Learn how to assess the significance of overlaps between CNVs and functional
  genomic regions
* Learn how carry out association analysis with gene expression and quantitative
  phenotypes
* Gain familiarity with allele-specific absolute CN analysis of cancer genomic data 
* Understand how CNVs can be experimentally detected and computationally inferred
  from SNP arrays and next-generation sequencing data
 
# Specific objectives

* Use `GRangesList` and `RaggedExperiment` to represent, access, and 
  manipulate CNV data 
* Identify recurrent CNV regions in a population,
  including density trimming, reciprocal overlap, and recurrence significance estimation
* Use the [regioneR](http://bioconductor.org/packages/regioneR) package
  to assess the significance of overlaps between CNVs and functional genomic regions
  such as genes, promoters, and enhancers.
* Carry out eQTL analysis for CNV and RNA-seq data using `GenomicRanges` / `RaggedExperiment` / `CNVRanger` architecture
* Carry out a GWAS analysis for CNV and quantitative phenotype data
* Perform estimation of tumor purity and ploidy from absolute CN analysis with 
  [PureCN](http://bioconductor.org/packages/PureCN)
  
