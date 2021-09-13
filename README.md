# PEmRNA-seq-nf

A Nextflow pipeline used for quantification and quanlity control for PE mRNA-seq data.

## Execution of pipeline using Nextflow
```
git clone https://github.com/AndersenLab/PEmRNA-seq-nf.git

cd PEmRNA-seq-nf

nextflow PEmRNAseq.nf --ref=c_elegans.PRJNA13758.WS276.genome.fa.gz --vcf=WI.20200815.hard-filter.vcf.gz --gtf=c_elegans.PRJNA13758.WS276.canonical_geneset.gtf --teref=bin/cele_Dfam_157_singleline.fasta --fqs=test.tsv 
```


## Pipeliine parameters

* --fqs

We use a sample sheet as the input of sequences here, see the example in `test.tsv`.

* --ref

Reference genome from WormBase

* --vcf

VCF file from CeNDR. 
This VCF file and the reference genome are used to generate SNVs-substituted genome for each strain.

* --gtf

Unzipped GTF file from WormBase.
The GTF file and the strain-specific genome are used to generate strain-specific transcriptome.


* --teref
Transcriptome of transposons. The file in the bin folder is generated using script here: https://github.com/fansalon/TEconsensus#dfam



