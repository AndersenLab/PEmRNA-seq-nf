#!/usr/bin/env nextflow



date = new Date().format( 'yyyyMMdd' )
params.out = "PEmRNAseq-${date}"
params.email = ""



 



/*
~ ~ ~ > * FastQ path file
*/


params.fqs = "$baseDir/test.tsv"

File fq_file = new File(params.fqs)


fileinfo = Channel
    .from(fq_file.collect { it.tokenize("\t") })
    .map { strain, SM, reads_1, reads_2, pol -> [ strain, SM, file("${reads_1}"), file("${reads_2}"), pol ] }
    .into { fq_file_1;fq_file_2 }

Channel
    .from(fq_file.collect { it.tokenize("\t") })
    .map { ST, SM, reads_1, reads_2, pol -> ST  }
    .unique()
    .into { Strains; 
            Strains2}

 


/*
~ ~ ~ > * reference genome from WormBase
*/


params.ref                    = null

File reference = new File("${params.ref}")
reference_handle = reference.getAbsolutePath()




/*
~ ~ ~ > * VCF from CeNDR
*/

params.vcf                    = null

File vcf = new File("${params.vcf}")
vcf_handle = vcf.getAbsolutePath()



/*
~ ~ ~ > * reference GTF from WormBase
*/


params.gtf                    = null

File gtf = new File("${params.gtf}")
gtf_handle = gtf.getAbsolutePath()



/*
~ ~ ~ > * TE ref from Dfam, https://github.com/fansalon/TEconsensus#dfam
*/


params.teref                    = null

File tereference = new File("${params.teref}")
teref_handle = tereference.getAbsolutePath()










/* 
   ==================================
   pre trim QC
   ==================================
*/ 


 


process pre_trim_fastqc {


    errorStrategy 'retry'

    tag "${ST}_${SM}"

    cpus 8

    input:
        set ST, SM, file(reads_1), file(reads_2), pol from fq_file_2

    output:

        set file("${ST}_${SM}_${pol}_prelog_1/*.zip"), file("${ST}_${SM}_${pol}_prelog_2/*.zip") optional true into prefastqc_ch

    script:
        """

        
        mkdir -p ${ST}_${SM}_${pol}_prelog_1

        fastqc -t 4 -o ${ST}_${SM}_${pol}_prelog_1 -f fastq -q ${reads_1}

        mkdir -p ${ST}_${SM}_${pol}_prelog_2

        fastqc -t 4 -o ${ST}_${SM}_${pol}_prelog_2 -f fastq -q ${reads_2}

       
        """
}


process summary_multi_qc_pre {

    errorStrategy 'retry'

    publishDir "${params.out}/multiqc_report", mode: 'copy'

    cpus 4
    memory '64 GB'

    input:

        file('*') from prefastqc_ch.collect()


    output:

        file("multiqc_fastqc_pretrim.html") 

    script:
    
        """


        multiqc *_fastqc.zip --cl_config '{max_table_rows: 1600}' --filename multiqc_fastqc_pretrim.html --interactive


        """
}







/* 
   ==================================
   trim
   ==================================
*/ 





process fastp_Trim {


    tag "reads: ${strain}_${SM}"

    errorStrategy 'retry'

    cpus 8 
    memory '32 GB'

    input:
       set val(strain), val(SM), file(reads_1), file(reads_2), val(pol) from fq_file_1

    
    output:
       set val(strain), val(SM), file ("${strain}_${SM}_${pol}_fastp1.fq.gz"),file ("${strain}_${SM}_${pol}_fastp2.fq.gz"), val(pol) into fq_file_trim1, fq_file_trim2

 


    script:
    //
    // fastp
    //
  
        """
        fastp --thread 8 --detect_adapter_for_pe --length_required 20 --in1 ${reads_1} --in2 ${reads_2} --out1 ${strain}_${SM}_${pol}_fastp1.fq.gz --out2 ${strain}_${SM}_${pol}_fastp2.fq.gz  
    
        """
    }









/*
   ==================================
    Generate Strain specific Transcriptome file
   ==================================
*/ 






process generate_specific_trancriptome {

   
    errorStrategy 'retry'
 
    tag "${ST}"

    cpus 4
    memory '32 GB'
    
    input:

      val(ST) from Strains
      
    output:
      
      set val(ST), file("${ST}.transcriptome2.fa") into sample_transcriptome 
 


    """

        bcftools consensus -f ${reference_handle} -i "TYPE='snp'" --sample ${ST} ${vcf_handle} > ${ST}.fa

        
        gffread -w ${ST}.transcriptome.fa -g ${ST}.fa ${gtf_handle}       


        cat ${ST}.transcriptome.fa | sed 's/ CDS=.*//' | sed 's/CDS://' | sed 's/Transcript://' | sed 's/Pseudogene://' > ${ST}.transcriptome3.fa


        cat ${ST}.transcriptome.fa ${teref_handle} > ${ST}.transcriptome2.fa
       
    """

}






/* 
   ==================================
   Kallisto index & mapping
   ==================================
*/ 


process kal_index {


    errorStrategy 'retry'

    cpus 1
    memory '2 GB'

    input:
        set val(ST), file(transcriptome_file) from sample_transcriptome


    output:
        set val(ST), file ("${ST}.transcriptome.index") into transcriptome_index

    script:
        //
        // Kallisto mapper index
        //
        """
        kallisto index -i ${ST}.transcriptome.index ${transcriptome_file}
        """
}


fq_file_trim1
	.combine(transcriptome_index, by: 0)
	.into{mapping_data_set;
             print_mapping}









process kal_mapping {

    publishDir "${params.out}/kallisto", mode: 'copy', pattern: "kallisto_*"



    errorStrategy 'retry'

    tag "reads: ${strain}_${SM}_${pol}"

    cpus 4
     memory '16 GB'

    input:
       set val(strain), val(SM), file(reads_1), file(reads_2), val(pol), file(kalIndex) from mapping_data_set
       


    output:
       set val(strain), val(SM), val(pol), file ("kallisto_${strain}_${SM}_${pol}") into kallisto_out_dirs

        file("*_log.txt") into kallisto_log

   

    script:
    //
    // Kallisto tools mapper
    //
  
        """
        mkdir kallisto_${strain}_${SM}_${pol}

        kallisto quant --bootstrap 100 -i ${kalIndex} -o kallisto_${strain}_${SM}_${pol} ${reads_1} ${reads_2} --threads=4 &> ${strain}_${SM}_${pol}_kallisto_log.txt


        
    
        """
    }





 






/* 
   ==================================
   post trim QC
   ==================================
*/ 








    
process post_trim_fastqc {

 
    errorStrategy 'retry'

    tag "${ST}_${SM}"

    cpus 8

    input:
        set ST, SM, file(reads_1), file(reads_2), pol from fq_file_trim2

    output:

        set file("${ST}_${SM}_${pol}_log_1"), file("${ST}_${SM}_${pol}_log_2") 

        set file("${ST}_${SM}_${pol}_log_1/*.zip"), file("${ST}_${SM}_${pol}_log_2/*.zip") into fastqc_ch

    script:
        """

        
        mkdir -p ${ST}_${SM}_${pol}_log_1

        fastqc -t 4 -o ${ST}_${SM}_${pol}_log_1 -f fastq -q ${reads_1}

        mkdir -p ${ST}_${SM}_${pol}_log_2

        fastqc -t 4 -o ${ST}_${SM}_${pol}_log_2 -f fastq -q ${reads_2}



       
        """
}




fastqc_ch
    .mix(kallisto_log)
    .collect()
    .into{qc_data;
    qc_data_2 }


 
 

process summary_multi_qc {

    errorStrategy 'retry'

    publishDir "${params.out}/multiqc_report", mode: 'copy'

    cpus 4
    memory '64 GB'

    input:

        file('*') from qc_data


    output:

         set file("multiqc_fastqc.html"), file("multiqc_kallisto.html") 

    script:
    
        """

        multiqc *_fastqc.zip --cl_config '{max_table_rows: 1600}' --filename multiqc_fastqc.html --interactive


        multiqc *_kallisto_log.txt --cl_config '{max_table_rows: 800}' --filename multiqc_kallisto.html --interactive

        """
}



















// Pipeline work summary, there would be a email notification if email has been given.
workflow.onComplete {

    summary = """   
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    Error report: ${workflow.errorReport ?: '-'}
    Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
    """

    println summary

    // mail summary
    if (params.email) {
        ['mail', '-s', 'eQTL-nf', params.email].execute() << summary
    }


}
