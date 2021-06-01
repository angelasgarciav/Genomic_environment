params.output_dir = false

process RUN_BLAST_TO_TXT {
    tag "$sample_id"
    publishDir "${params.output_dir}/blast_files", mode: 'copy'

    input:
        tuple(val(sample_id), path(genome_file))
        file(reference_gene)

    output:
        tuple(val(sample_id), path(genome_file), path("${sample_id}_${reference_name}.txt"), emit: blast_file_tuple)
        tuple(val(sample_id), path("${sample_id}_${reference_name}.txt"), emit: blast_file)
        path("${sample_id}_${reference_name}.txt")

    script:
    reference_name = reference_gene.baseName.replaceAll(/\..+$/,'')
    """
    #!/usr/bin/env python3
    import sys
    sys.path.append("${workflow.projectDir}/modules/")
    import blast_from_co_located as blast_py
    blast_py.run_onesample_blast_to_txt("${genome_file}", "${reference_gene}", 95, '${sample_id}_${reference_name}.txt')
    """
}

process EXTRACT_CONTIG_FILE {
    tag "$sample_id"
    publishDir "${params.output_dir}/contigs", mode: 'copy'

    input:
        tuple(val(sample_id), path(genome_file), path(blast_file))

    output:
        tuple(val(sample_id), path("${sample_id}_contig.fasta"))

    script:
    """
    awk -F \$'\\t' 'FNR==2 {print \$4}' ${blast_file} | seqtk subseq ${genome_file} - > ${sample_id}_contig.fasta
    """
}

process RUN_PROKKA {
    tag "$sample_id"
    publishDir params.output_dir, mode: 'copy'
    input:
        tuple(val(sample_id), path(fastas))

    output:
        path("${sample_id}/${sample_id}.tsv")
        path("${sample_id}/${sample_id}.txt")
        tuple(val(sample_id), path("${sample_id}/${sample_id}.gff"), emit: gene_annotation)
    script:
    """
    prokka --outdir $sample_id --prefix ${sample_id} ${fastas} --metagenome
    """
}

process GFF_TO_BED {
    tag "$sample_id"

    input:
    tuple(val(sample_id), path(annotated_genome_file))
    
    output:
    tuple(val(sample_id), path("${sample_id}.bed"))

    script:
    """
    sed -n '/##FASTA/q;p' ${annotated_genome_file} > ${sample_id}.bed
    """
}

// process LOCATE_GENE_IN_BED {
//     tag "$sample_id"

//     input:
//     tuple(val(sample_id), path(bed_file))
//     tuple(val(sample_id), path(blast_file))

//     output:
//     tuple(val(sample_id), path(bed_file), val(position), val(prokka_tag))

//     script:
//     """
//     #!/usr/bin/env Rscript

//     gff_file <- read.delim("${bed_file}", header = F, comment.char = "#")
//     blast_file <- read.table("${blast_file}", sep = "\t", header = T)

//     start = blast_file\$subject_start_1
//     end = blast_file\$subject_end_1
//     range <- abs(end - start)
//     min_data = min(start, end)

//     # id closest match to start
//     position <- which.min(abs(gff_file[,2]-min_data))
//     which(abs(gff_file[,2]-min_data)==min(abs(gff_file[,2]-min_data)))

//     ${position} = position
//     ${prokka_tag} <- gff_file[position, 1]
//     """
// }

//process EXTRACT_nGENES_AROUND_xGEN {
//    tag "$sample_id"
//
//    input: 
//      tuple(val(sample_id), path(annotated_file))
//      val(gene_name)
//      file ("params.output_dir/${sample_id}/${sample_id}.tsv")
//      file (¨${sample_id}.tsv¨)
//    output:
//      path ("${sample_id}_extracted.txt")
//    script:
//    """
//      grep 
//    """
//}
