# ConservationProfiler
This script creates multiple sequence alignments (MSAs) from the BLAST results of input sequences against input proteomes/genomes

# Rquirements:
-Biopython (tested on 1.70)<br>
-Numpy (tested on 1.13.1)
<br><br>
# Description:
This script creates multiple sequence alignments (MSAs) for the BLAST results of input sequences against a list of target proteomes/genomes.<br>
An alignment is created for each sequence as it is BLASTed or Supermatched against the target proteomes/genomes.
The resulting alignments are then filtered using several threshold, and an MSA is outputted using several multiple sequence aligners for the remaining sequences.<br>
The inputs:<br>
proteome_file:      Reference proteome/genome file. Must be in fasta format, can be gzipped.<br>
genes_data:         Tab delimited genes_data file, each row must be one sequence, columns must be sequence name, 1-based start (can just be
                    "start"), 1-based end (can just be "end"), and can also include the id of the sequence (for protein data), the chromosome
                    of the sequence and the strand (for nucleotide data).<br>
proteome_list:      List of all the target proteomes/genomes, the name must be similar to the names of the files or the directories in which are saved. Should have "name" at the top.<br>
                    (you can have spaces or slashes instead of underlines, but don't get too creative).<br>
database_dir:       The directory where the proteomes/genomes are saved, they should be in fasta format. The proteomes/genomes can all be saved as files
                    in one directory, or each proteome/genome can be in its own subdirectory. If you are intending to run BLAST mode, blast dbs
                    must be created a priori!<br>
workdir:            Directory where the script will create and edit files, doesn't have to be an existing directory.<br>
file_type:          Filetype for the target proteomes/genomes, the defualt is .fna.gz.<br>
exclude:            Proteomes/Genomes including any string inputted here will be excluded.<br>
source:             Proteomes or genomes? Input can be Protein/Proteome (default) or Nucleotide/Genome.<br>
mode:               What algorithm should be used for sequence alignment step. Script supports Supermatcher (default) and BLAST, but for
                    BLAST the user must create blast dbs a priori.<br>
identity_threshold: Identity percentage minimum threshold to query for alignments to enter the multiple sequence alignment. Default is 0.59.<br>
eval_threshold:     Maximum E-value threshold for query so that alignments can enter the multiple sequence alignment. Default is 0.05.<br>
gene_percentage:    Minimum fraction of BLAST/Supermatcher hit length / query length for alignments to enter the multiple sequence alignment. Default is 0.5.<br>
per_proteome_sequences:
                    How many sequences will be extracted form the BLAST/Supermatcher results per proteome/genome. The default is None
                    meaning all of them, if you want to change it, please input an integer.<br>
msa:                What multiple sequence alignment method should be used to construct the MSA, options are: Clustal Omega, ClustalW2, emma, MSUCLE, MAFFT,
                    Dialign, MSAProbs, Probcons, PRANK, T-Coffee.<br><br>

The outputs:<br>
workdir/$seq/$seq.fasta:      The sequence that was extracted for $seq from the reference proteome/genome<br>
workdir/$seq/msa.fasta:       All the resulting hits of the BLAST/supermatcher search for the sequence<br>
workdir/$seq/msa.aln:         Multiple Sequence Alignment for the sequence, format is clustal unless emma or mafft are used as the
                              algorithm in which case it will be fasta. If you want emma in clustal format just use clustalw2,
                              if you want MAFFT in clustal format I would consider MUSCLE, or just converting the resulting
                              fasta to clustal format.<br>
workdir/$seq/msa.dnd:         Dendrogram built based on the alignment<br>
workdir/$seq/seq_scores.csv:  File summarizing BLAST/Supermatcher results with the name, description, e-value, identity percentage,
                              score, sequence, length and fraction of hit length/query length for all target genomes/proteomes<br>
workdir/$seq/seq_scores_statistics.csv:
                              File with statistics about the seq scores of all alignments. You can see the distributuion of all of the measures in teh seq_scores file across all sequences.<br>
