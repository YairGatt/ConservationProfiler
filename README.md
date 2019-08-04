# ConservationProfiler
This script creates multiple sequence alignments (MSAs) from the BLAST results of input sequences against input proteomes/genomes

# Rquirements:
-Biopython (tested on 1.70)<br>
-Numpy (tested on 1.13.1)

# Description:
This script creates multiple sequence alignments (MSAs) for the BLAST results of input sequences against a list of target proteomes/genomes.<br>
An alignment is created for each sequence as it is BLASTed or Supermatched against the target proteomes/genomes.
The resulting alignments are then filtered using several threshold, and an MSA is outputted using several multiple sequence aligners for the remaining sequences.<br>
<u><b>The inputs:</b></u><br>
<b>proteome_file:</b>      Reference proteome/genome file. Must be in fasta format, can be gzipped.<br>
<b>genes_data:</b>         Tab delimited genes_data file, each row must be one sequence, columns must be sequence name, 1-based start (can just be
                    "start"), 1-based end (can just be "end"), and can also include the id of the sequence (for protein data), the chromosome
                    of the sequence and the strand (for nucleotide data).<br>
<b>proteome_list:</b>      List of all the target proteomes/genomes, the name must be similar to the names of the files or the directories in which are saved. Should have "name" at the top.<br>
                    (you can have spaces or slashes instead of underlines, but don't get too creative).<br>
<b>database_dir:</b>       The directory where the proteomes/genomes are saved, they should be in fasta format. The proteomes/genomes can all be saved as files
                    in one directory, or each proteome/genome can be in its own subdirectory. If you are intending to run BLAST mode, blast dbs
                    must be created a priori!<br>
<b>workdir:</b>            Directory where the script will create and edit files, doesn't have to be an existing directory.<br>
<b>file_type:</b>          Filetype for the target proteomes/genomes, the defualt is .fna.gz.<br>
<b>exclude:</b>            Proteomes/Genomes including any string inputted here will be excluded.<br>
<b>source:</b>             Proteomes or genomes? Input can be Protein/Proteome (default) or Nucleotide/Genome.<br>
<b>mode:</b>               What algorithm should be used for sequence alignment step. Script supports Supermatcher (default) and BLAST, but for
                    BLAST the user must create blast dbs a priori.<br>
<b>identity_threshold:</b> Identity percentage minimum threshold to query for alignments to enter the multiple sequence alignment. Default is 0.59.<br>
<b>eval_threshold:</b>     Maximum E-value threshold for query so that alignments can enter the multiple sequence alignment. Default is 0.05.<br>
<b>gene_percentage:</b>    Minimum fraction of BLAST/Supermatcher hit length / query length for alignments to enter the multiple sequence alignment. Default is 0.5.<br>
<b>per_proteome_sequences:</b>
                    How many sequences will be extracted form the BLAST/Supermatcher results per proteome/genome. The default is None
                    meaning all of them, if you want to change it, please input an integer.<br>
<b>msa:</b>                What multiple sequence alignment method should be used to construct the MSA, options are: Clustal Omega, ClustalW2, emma, MSUCLE, MAFFT,
                    Dialign, MSAProbs, Probcons, PRANK, T-Coffee.<br><br>

<b><u>The outputs:</b></u><br>
<b>workdir/$seq/$seq.fasta:</b>      The sequence that was extracted for $seq from the reference proteome/genome<br>
<b>workdir/$seq/msa.fasta:</b>       All the resulting hits of the BLAST/supermatcher search for the sequence<br>
<b>workdir/$seq/msa.aln:</b>         Multiple Sequence Alignment for the sequence, format is clustal unless emma or mafft are used as the
                              algorithm in which case it will be fasta. If you want emma in clustal format just use clustalw2,
                              if you want MAFFT in clustal format I would consider MUSCLE, or just converting the resulting
                              fasta to clustal format.<br>
<b>workdir/$seq/msa.dnd:</b>         Dendrogram built based on the alignment<br>
<b>workdir/$seq/seq_scores.csv:</b>  File summarizing BLAST/Supermatcher results with the name, description, e-value, identity percentage,
                              score, sequence, length and fraction of hit length/query length for all target genomes/proteomes<br>
<b>workdir/$seq/seq_scores_statistics.csv:</b>
                              File with statistics about the seq scores of all alignments. You can see the distributuion of all of the measures in teh seq_scores file across all sequences.<br>
