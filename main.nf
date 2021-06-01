nextflow.enable.dsl=2

params.fastas_dir = false
params.num_genes = 5
params.output_dir = false
params.reference_gene = false

include { RUN_PROKKA; EXTRACT_CONTIG_FILE; RUN_BLAST_TO_TXT; GFF_TO_BED} from './modules/process.nf'

workflow {
    if (params.fastas_dir && params.num_genes && params.output_dir) {
        fastas = Channel
        .fromPath(params.fastas_dir)
        .map{ file -> tuple (file.baseName.replaceAll(/\..+$/,''), file)}
        .ifEmpty { error "Cannot find any fastas matching ${params.fastas_dir}"}
    
    fastas.view()
    RUN_BLAST_TO_TXT(fastas, file(params.reference_gene)) 
    //RUN_BLAST_TO_TXT.out.blast_file_tuple.view()
    EXTRACT_CONTIG_FILE(RUN_BLAST_TO_TXT.out.blast_file_tuple)
    RUN_PROKKA(EXTRACT_CONTIG_FILE.out)
    RUN_PROKKA.out.gene_annotation.view()
    GFF_TO_BED(RUN_PROKKA.out.gene_annotation)
    //LOCATE_GENE_IN_BED(GFF_TO_BED.out, RUN_BLAST_TO_TXT.out.blast_file)

    } else {
        error "Please specify a fastas directory '--fastas_dir fastas/*.fasta' and an output directory '--output_dir output'" 
    }
}