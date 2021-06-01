# Script to show from blast_file and gff_file, locate the first blast file in the gff
# input1: 2-line file of one gene located through blast
# input2: annotated file.gff without ##FASTA information

locate_gene_position  <- function (annotation_file, blast_file) {
  gff_file <- read.delim(file = annotation_file, header = F, comment.char = "#")
  blast_file <- read.table(file = blast_file, sep = "\t", header = T)

  start = blast_file[,5]
  end = blast_file[,6]
  range <- abs(end - start)
  min_data = min(start, end)
  position <- which.min(abs(gff_file[,2]-min_data))
  return(position)
}

# open dir files with pattern .bed or .txt

# extract name to match files

# id position of my gene of interest

# get surrounding n number, gene names

# parse the output per sample id


# load files
gff_file <- read.delim("tests/genes/PROKKA_05182021/test.bed", header = F, comment.char = "#")
gff_file <- read.delim("output/G18001493/G18001493.bed", header = F, comment.char = "#")
blast_file <- read.table("output/blast_files/G18001493_NG_050347.txt", sep = "\t", header = T)

start = blast_file$subject_start_1
end = blast_file$subject_end_1
range <- abs(end - start)

# other possile way
min_data = min(start, end)
if (min_data == start) {
  range=(end - start)  
} else {
  range=(start - end)
}

# id closest match to start
position <- which.min(abs(gff_file[,2]-min_data))
which(abs(gff_file[,2]-min_data)==min(abs(gff_file[,2]-min_data)))

gff_gene_tag <- gff_file[position, 1]

x=c(1:10^6)
your.number=90000.67
which(abs(x-your.number)==min(abs(x-your.number)))
