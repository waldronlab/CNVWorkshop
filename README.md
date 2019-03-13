# Copy number variation analysis with Bioconductor

# Instructor(s) name(s) and contact information

Ludwig Geistlinger, Marcel Ramos, Seyhun Oh, and Levi Waldron

CUNY School of Public Health
55 W 125th St, New York, NY 10027

Ludwig.Geistlinger@sph.cuny.edu
Marcel.Ramos@sph.cuny.edu
Seyhun.Oh@sph.cuny.edu
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
absolute copy number analysis and how it is incorporated in genomic cancer analysis 
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

* get familiar with elementary concepts of CNV analysis
* learn how to efficiently represent, access, and manipulate CNV data 
  in Bioconductor data structures
* get familiar with different strategies for summarizing individual CNV
  calls across a population
* learn how to assess the significance of overlaps between CNVs and functional
  genomic regions
* learn how carry out association analysis with gene expression and quantitative
  phenotypes
* get familiar with allele-specific absolute CN analysis of genomic cancer data 
 
# Specific objectives

* understand how CNVs can be experimentally detected and computationally inferred
  from SNP arrays and next-generation sequencing data
* learn how to use `GRangesList` and `RaggedExperiment` to represent, access, and 
  manipulate CNV data 
* understand different strategies for finding recurrent CNV regions in a population,
  including density trimming, reciprocal overlap, and recurrence significance estimation
* learn how to use the [regioneR](http://bioconductor.org/packages/regioneR) package
  to assess the significance of overlaps between CNVs and functional genomic regions
  such as genes, promoters, and enhancers.
* learn how to carry out eQTL analysis for CNV and RNA-seq data
* learn how to carry out a GWAS analysis for CNV and quantitative phenotype data
* learn how to estimate tumor purity and ploidy from absolute CN analysis with 
  [PureCN](http://bioconductor.org/packages/PureCN)
  
