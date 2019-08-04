"""
This script creates multiple sequence alignments (MSAs) for the BLAST results of input sequences against a list of target proteomes/genomes.
An alignment is created for each sequence as it is BLASTed or Supermatched against the target proteomes/genomes,
The resulting alignments are then filtered using several threshold, and an MSA is outputted using several multiple sequence aligners for the remaining
sequences.
The inputs:
proteome_file:      Reference proteome/genome file. Must be in fasta format, can be gzipped.
genes_data:         Tab delimited genes_data file, each row must be one sequence, columns must be sequence name, 1-based start (can just be
                    "start"), 1-based end (can just be "end"), and can also include the id of the sequence (for protein data), the chromosome
                    of the sequence and the strand (for nucleotide data).
proteome_list:      List of all the target proteomes/genomes, the name must be similar to the names of the files or the directories in which are saved. Should have "name" at the top.
                    (you can have spaces or slashes instead of underlines, but don't get too creative).
database_dir:       The directory where the proteomes/genomes are saved, they should be in fasta format. The proteomes/genomes can all be saved as files
                    in one directory, or each proteome/genome can be in its own subdirectory. If you are intending to run BLAST mode, blast dbs
                    must be created a priori!
workdir:            Directory where the script will create and edit files, doesn't have to be an existing directory.
file_type:          Filetype for the target proteomes/genomes, the defualt is .fna.gz.
exclude:            Proteomes/Genomes including any string inputted here will be excluded.
source:             Proteomes or genomes? Input can be Protein/Proteome (default) or Nucleotide/Genome.
mode:               What algorithm should be used for sequence alignment step. Script supports Supermatcher (default) and BLAST, but for
                    BLAST the user must create blast dbs a priori.
identity_threshold: Identity percentage minimum threshold to query for alignments to enter the multiple sequence alignment. Default is 0.59.
eval_threshold:     Maximum E-value threshold for query so that alignments can enter the multiple sequence alignment. Default is 0.05.
gene_percentage:    Minimum fraction of BLAST/Supermatcher hit length / query length for alignments to enter the multiple sequence alignment. Default is 0.5.
per_proteome_sequences:
                    How many sequences will be extracted form the BLAST/Supermatcher results per proteome/genome. The default is None
                    meaning all of them, if you want to change it, please input an integer.
msa:                What multiple sequence alignment method should be used to construct the MSA, options are: Clustal Omega, ClustalW2, emma, MSUCLE, MAFFT,
                    Dialign, MSAProbs, Probcons, PRANK, T-Coffee.

The outputs:
workdir/$seq/$seq.fasta:      The sequence that was extracted for $seq from the reference proteome/genome
workdir/$seq/msa.fasta:       All the resulting hits of the BLAST/supermatcher search for the sequence
workdir/$seq/msa.aln:         Multiple Sequence Alignment for the sequence, format is clustal unless emma or mafft are used as the
                              algorithm in which case it will be fasta. If you want emma in clustal format just use clustalw2,
                              if you want MAFFT in clustal format I would consider MUSCLE, or just converting the resulting
                              fasta to clustal format.
workdir/$seq/msa.dnd:         Dendrogram built based on the alignment
workdir/$seq/seq_scores.csv:  File summarizing BLAST/Supermatcher results with the name, description, e-value, identity percentage,
                              score, sequence, length and fraction of hit length/query length for all target genomes/proteomes
workdir/$seq/seq_scores_statistics.csv:
                              File with statistics about the seq scores of all alignments. You can see the distributuion of all of the measures in teh seq_scores file across all sequences.
"""

#general
import sys, os
import optparse
import csv
import errno
import glob
import gzip as gz
import numpy as np
from multiprocessing import Pool
from itertools import repeat
import __future__
import warnings
import traceback
#Biopython stuff
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO
from Bio.Application import ApplicationError
from Utils.Emboss import AlignIO
#BLAST stuff
from Bio.Blast.Applications import NcbiblastpCommandline, NcbiblastnCommandline
from Bio.Blast import NCBIXML
#Multiple alignment stuff
from Bio.Align.Applications import ClustalOmegaCommandline, MuscleCommandline, MafftCommandline, TCoffeeCommandline, ClustalwCommandline, PrankCommandline, MSAProbsCommandline, ProbconsCommandline, DialignCommandline
from Utils.Emboss.EmbossCommands import EmmaCommandline, SuperMatcherCommandline
from Exceptions.NoMatchForSeqException import NoMatchForSeqException

def process_command_line(argv):
    """
    Return a 2-tuple: (settings object, args list).
    `argv` is a list of arguments, or `None` for ``sys.argv[1:]``.
    """
    if argv is None:
        argv = sys.argv[1:]

    # initialize the parser object:
    parser = optparse.OptionParser(
        formatter=optparse.TitledHelpFormatter(width=78),
        add_help_option=None)

    parser.add_option(
        '-p' , '--proteome_file',
        help='Path to the reference proteome file.')
    
    parser.add_option(
        '-d' , '--genes_data',
        help='Path to file containing the data of the sequences to be processed.')

    parser.add_option(
        '-l' , '--proteome_list',
        help='Path to file containing all the proteomes against which the sequences will be searched.')

    parser.add_option(
        '-D' , '--database_dir', default="/home/hosts/disk20/microbial" ,
        help='Path to directory holding the proteomes.')

    parser.add_option(
        '-w' , '--workdir', default="./temp/",
        help='Path where the script can create and edit files.')
    
    parser.add_option(
        '-t' , '--file_type', default=".fna.gz",
        help='File type of proteome/genome files.')
    
    parser.add_option(
        '-x' , '--exclude', default=None,
        help='Genome/proteome files including string will be excluded.')

    parser.add_option(
        '-s' , '--source', default="Protein",
        help='Proteomes or genomes? Input can be Protein/Proteome or Nucleotide/Genome')
    
    parser.add_option(
        '-m' , '--mode', default="BLAST",
        help='What algorithm should be used for sequence alignment. script supports Supermatcher and BLAST, but for BALST the user must supply genomes/proteomes with blast dbs.')
    
    parser.add_option(
        '-i' , '--identity_threshold', type="float", default=0.59,
        help='Minimum identity percentage to query for significant sequences.')
    
    parser.add_option(
        '-e' , '--eval_threshold', type="float", default=0.05,
        help='Only for BLAST mode. Minimum E-value against query for significant sequences.')

    parser.add_option(
        '-c' , '--gene_percentage', type="float", default=0.5,
        help='Only for BLAST mode. What percentage of the length of the query is the minimum length of the BLAST/Supermatcher hit for significant sequences.')

    parser.add_option(
        '-g' , '--per_proteome_sequences', type="int", default=None,
        help='How many sequences to keep from each proteome/genome, default is None, meaning all of them. Input is integer.')
    
    parser.add_option(
        '-a' , '--msa', default="clustalo",
        help='What multiple sequence alignment method should be used, options are: Clustal Omega, ClustalW2, emma, MSUCLE, MAFFT, Dialign, MSAProbs, Probcons, PRANK, T-Coffee.')

    parser.add_option(    #customized description; put --help last
        '-h', '--help', action='help',
        help='Show this help message and exit.')

    settings, args = parser.parse_args(argv)

    #check number of arguments, verify values, etc.:
    if args:
        parser.error('Program takes no command-line arguments; '
                     '"%s" ignored.' % (args,))

    #raise input errors if necessary
    if not settings.proteome_file:
        parser.error("Missing proteome_file.")

    if not settings.genes_data:
        parser.error("Missing genes_data file.")

    if not settings.proteome_list:
        parser.error("Missing proteome_list file.")
    #return
    return settings, args

def extract_genes_data(genes_data_file):
    """
    This function extracts the data of the sequences from the tab delimited file, it can deal with different formats, but the file
    must have one sequence per row, and each sequence must have 1-based start and end coordinates, which can also just be "start"
    or "end". The sequence will later be extracted from the relevant id/chr sequence in the reference genome (if none is
    inputted the first one will be used arbitrarily) from the start coordinate to the end coordinate according to the
    chosen strand (for nucleotide data).
    """
    #initialize
    gene_data_dict = {}
    #open file
    with open(genes_data_file, "rb") as fl:
        #read as csv
        reader = csv.DictReader(fl, delimiter="\t")
        #iterate through rows
        for row in reader:
            #if protein
            if "id" in row.keys(): seq_specification_field = row["id"] #id field was given and specifies the relevant protein
            #if nucleotide
            elif "chr" in row.keys(): seq_specification_field = row["chr"] #chr field was given and specifies the relevant chromosome
            elif "chromosome" in row.keys(): seq_specification_field = row["chromosome"] #chromosome field was given and specifies the relevant chromosome
            else: seq_specification_field = None
            #get start coordinate
            try:
                start = min(int(row["min_start"]),int(row["max_end"]))
            except ValueError:
                start = "start"
            #get end coordinate
            try:
                end = max(int(row["max_end"]),int(row["min_start"]))
            except ValueError:
                end = "end"
            #assess format
            if len(row) == 5: #data includes id/chr and also strand
                gene_data_dict[row["name"]] = (start, end, seq_specification_field, row["strand"])
            elif len(row) == 4: #data includes only id/chr or only strand
                if seq_specification_field:
                    gene_data_dict[row["name"]] = (start, end, seq_specification_field)
                else:
                    gene_data_dict[row["name"]] = (start, end, row["strand"])
            elif len(row) == 3: #data includs only start and end coordinates
                gene_data_dict[row["name"]] = (start, end)
            else:
                raise BaseException("genes_data file must contain a name, a start and an end for each row")
    #return
    return gene_data_dict

def extract_seq(data, proteome_seq_list):
    """
    This function extracts the sequence of sequences from the reference proteome/genome.
    """
    #Default is no id/chr and no strand
    seq_specification_field = None
    strand = None
    #coordinates
    start, end = data[0:2]
    #process end
    if end == "end": end = None
    #process start
    if start == "start": start = 1
    #get data format
    try: #Data is in the form start,end, seq_spec, strand
        if data[3] in [-1,1,"-","+","plus","minus", "True", "False"]: #checking for strand at final position
            strand = data[3]
            seq_specification_field = data[2]
    except IndexError: #Data has only 2 or 3 fields
        try:
            if data[2] in [0,1,"-","+","plus","minus", "True", "False"]: #checking for strand at final position
                strand = data[2] #there is only a strand and no id/chr
            else: #final field is seq_spec
                seq_specification_field = data[2]
        except IndexError: #Data has only 2 fields which are start and end coordinates
            pass
    #get sequence
    if seq_specification_field != None:
        #iterate
        for prot in proteome_seq_list:
            #seearch for match for seq specification field
            if seq_specification_field in prot.id:
                #get seq
                seq = prot.seq[start-1:end]
    #get sequence from first record
    elif seq_specification_field == None:
        prot = proteome_seq_list[0]
        seq = prot.seq[start-1:end]
    #flip if strand minus
    if strand in (-1,"-","minus","False"):
        seq = seq.reverse_complement()
    #return
    return seq
    
def mkdir(path):
    #Create directory and don't crash if it already exists
    try:
        os.mkdir(path)
    except OSError as exc:
        if exc.errno != errno.EEXIST:
            raise exc
        pass

def make_fasta_per_gene(path, name, seq):
    """
    Create a fasta file for sequence so it can be used as a query file for BALST/Supermatcher
    """
    #name
    seq_file_name = os.path.join(path, "%s.fasta" % name)
    #record
    record = SeqRecord(seq, id=name, name=name)
    #write
    SeqIO.write(record, seq_file_name, "fasta")
    #return
    return seq_file_name

def extract_proteome_list(proteome_list_file):
    """
    Get a list of proteomes/genomes from the proteome_list file
    """
    with open(proteome_list_file, "rb") as fl:
        try:
            reader = csv.DictReader(fl, delimiter="\t")
            proteome_name_list = [row["name"].replace(" ","_").replace("/","_") for row in reader]
        except KeyError:
            proteome_list_file = [i.strip() for i in fl.readlines()]
    #return
    return proteome_name_list

def parse_supermatcher_result(file_name, number):
    """
    Get identity percentage and supermatcher score values from the supermatcher result file, number signifies the number
    of the alignment (first,second, etc.) and it must be a 0-based int.
    """
    #initialize
    identity_percentage = None
    score = None
    #open file
    with open(file_name, "rb") as fl:
        #get first line
        line = fl.readline()
        #initialize lists
        identity_percentage_list = []
        score_list = []
        #iterate through lines
        while (line != ""):
            #get identity
            if line.startswith("# Identity:"):
                #this linejust converts the identity format to a proper float
                identity_percentage_list.append(eval(compile(line.split("# Identity:")[1].split(" (")[0].strip(), '<string>', 'eval', __future__.division.compiler_flag)))
            #get score
            if line.startswith("# Score: "):
                #Get score number as float
                score_list.append(float(line.replace("# Score: ", "")))
            #next line
            line = fl.readline()
    #return
    return score_list[number], identity_percentage_list[number] #return score and identity percentage of number # alignment

def directed_local_alignment(proteome_file, ref_seq_file, proteome_name, workdir, mode, source, strand, per_proteome_sequences):
    """
    This function does the main lifting, running Supermatcher or BLAST on your proteomes/genomes.
    """    
    #open result file
    alignment_file = os.path.join(workdir, "curr_alignment.aln")
    #If Supermatcher was chosen as the alignment algorithm
    if mode.lower() == "supermatcher":
        #Unlike BLAST, Supermatcher only searches one strand, so we create a temp file to use as the supermatcher
        #bsequence, and we can write into it either the positive or negative strand depending on the "strand" parameter
        records = open_proteome(proteome_file)
        if not records: return []
        if not strand: records = [i.reverse_complement(id=True, name=True, description=True, features=True, annotations=True, letter_annotations=True, dbxrefs=True) for i in records]
        #write file
        proteome_fasta_file = os.path.join(workdir, "curr_proteome.fasta")
        SeqIO.write(records, proteome_fasta_file, "fasta")
        #Build the matching command
        if source.lower() == "protein" or source.lower() == "proteome": matrix = "EBLOSUM62" #AA matrix
        elif source.lower() == "nucleotide" or source.lower() == "genome": matrix = "EDNAFULL"
        #run supermatcher
        cmd = SuperMatcherCommandline(asequence=ref_seq_file, 
                          bsequence=proteome_fasta_file, 
                          gapopen=10, 
                          gapextend=0.5,
                          datafile=matrix,
                          outfile=alignment_file)
        #Excecute the command
        stdout, stderr = cmd()
        #Parse the resulting alignments
        alignments=[]
        try:
            #Create list of MultipleSeqAlignment objects representing the Supermatcher results (list may be empty)
            align_seq_list = list(AlignIO.parse(alignment_file, "amir_emboss")) #List of MultipleSeqAlignment objects
            #Iterate through list, only the first alignment will be used if per_proteome_sequences==1 and all alignments
            #will be used if per_proteome_sequences==None
            for number,alignment in enumerate(align_seq_list[0:per_proteome_sequences]):
                #get the alignent
                align_seq = alignment[1] #SeqRecord objects, [0] is query and [1] is sbjct
                #remove gaps
                align_seq._set_seq(align_seq.seq.ungap("-"))
                #get name
                if (not per_proteome_sequences or per_proteome_sequences > 1) and len(align_seq_list) > 1:
                        usable_name = proteome_name + "_" + str(number)
                else:
                        usable_name = proteome_name
                align_seq.name = usable_name
                align_seq.id = usable_name
                #finalize parsing
                score, identity_percentage = parse_supermatcher_result(alignment_file, number)
                #add to list
                alignments.append((score, align_seq, identity_percentage,0,1)) #1 is given arbitrarily as gene_percentage
                #alignments is a list of tuples with each element being (score, align_seq, identity_percentage)
            #detlete temp files
            os.remove(alignment_file)
            os.remove(proteome_fasta_file)
            #return
            return alignments
        except ValueError or IndexError:
            raise NoMatchForSeqException(proteome_fasta_file, ref_seq_file)
    #If BLAST was chosen as the alignment algorithm
    elif mode.upper() == "BLAST":
        #if file is called XXXX.file_type.gz, db name should be XXXX        
        if proteome_file.endswith(".gz"): db_name_temp = ".".join(proteome_file.split(".")[:-2])
        #if file is called XXXX.file_type, db name should be XXXX
        else: db_name_temp = ".".join(proteome_file.split(".")[:-1])
        #define dir
        directory = os.path.dirname(proteome_file)
        directory_files = os.listdir(directory)
        #iterate through files and find database file
        for dir_file in directory_files:
                #determine type of die
                if os.path.basename(db_name_temp) in dir_file and ".nhr" in dir_file: db_name = directory+"/"+dir_file.split(".nhr")[0]
                elif os.path.basename(db_name_temp) in dir_file and ".phr" in dir_file: db_name = directory+"/"+dir_file.split(".phr")[0]
        #open file
        records = open_proteome(proteome_file)
        if not records: return []
        #Build matching command
        if source.lower() == "protein" or source.lower() == "proteome":
            cmd =  NcbiblastpCommandline(query=ref_seq_file, db=db_name, out=alignment_file, outfmt=5)
        elif source.lower() == "nucleotide" or source.lower() == "genome":
            cmd =  NcbiblastnCommandline(query=ref_seq_file, db=db_name, out=alignment_file, outfmt=5, task="blastn")
            #cmd =  NcbiblastnCommandline(query=ref_seq_file, db=db_name, out=alignment_file, outfmt=5)
        #Execute command
        stdout, stderr = cmd()
        #Open result
        try:
            result_handle = open(alignment_file)
            blast_record = list(NCBIXML.parse(result_handle))[0] #BLAST record object
            result_handle.close()
        except ValueError:
            raise NoMatchForSeqException(proteome_file, ref_seq_file)
        #Parse resulting alignments
        alignments = []
        try:
            #Iterate through list, only the first alignment will be used if per_proteome_sequences==1 and all alignments
            #will be used if per_proteome_sequences==None
            for number, alignment in enumerate(blast_record.alignments[0:per_proteome_sequences]):
                hsp = alignment.hsps[0] #HSP contains all the details about the alignment
                sequence = hsp.sbjct
                score = hsp.score
                evalue = hsp.expect
                identities = hsp.identities
                query_length = blast_record.query_letters
                align_length = hsp.align_length
                #calculate percentages
                identity_percentage = float(identities)/align_length
                if (not per_proteome_sequences or per_proteome_sequences > 1) and len(blast_record.alignments) > 1: name = proteome_name+"_"+str(number)
                else: name = proteome_name
                #length percentage
                percentage = float(align_length)/query_length
                #convert to SeqRecord0
                align_seq = SeqRecord(Seq(sequence,IUPAC.protein),id=name, name=name, description=name)
                align_seq._set_seq(align_seq.seq.ungap("-"))  #Remove the gaps
                #score, identity_percentage = parse_blast_result(alignment_file)
                alignments.append((score,align_seq,identity_percentage,evalue, percentage))
                #alignments is a list of tuples with each element being (score, align_seq, identity_percentage, evalue, length_percentage)
        except IndexError: #If alignments are empty, doesn't actually do anything
            sequence = ""
            score = 0
            identity_percentage = 0
        #remove temp file
        os.remove(alignment_file)
        #return
        return alignments
    #not BLAST or Supermatcher!
    else:
        raise BaseException("Only Supermatcher and BLAST modes are currently supported.")

def extract_aligned_seqs(proteome_name_list, ref_seq_file, database_dir, workdir, file_type, exclude, mode, source, per_proteome_sequences):
    """
    Function running directed_local_alignment for all proteomes/genomes and creates all the results
    """
    #define filetype properly
    if file_type[0] != ".": file_type = "." + file_type #file_type should start with period
    #initialzie
    best_aligned_seq_list = []
    #iterate through proteomes
    for proteome_name in proteome_name_list:
        #get the actual proteomes
        match_list = []
        #If the each proteome/genome is found in a subdirectory, we will use this to get all files of the file_type in the
        #subdirectory whose name is the name of the proteome/genome    
        if os.path.isdir("%s/%s" % (database_dir, proteome_name)): separator="/*"
        #If each proteome/genome is a single file we will use this to get all the flies of the file_type whose name is the
        #name of the proteome/genome    
        elif os.path.isfile("%s/%s%s" % (database_dir, proteome_name, file_type)): separator=""
        else: raise BaseException("No directory or file %s/%s of filetype %s" % (database_dir, proteome_name, file_type))
        #Supermatcher only searches the given strand, so we have to search each genome twice if it's nucleotide data
        if mode.lower() == "supermatcher" and (source.lower() == "nucleotide" or source.lower() == "genome"): strands=[True,False]
        #BLAST searches both strands automatically
        else: strands=[True]        
        #Go over the proteome/genome and plasmids
        #iterate through files and run!
        for proteome_file in glob.glob("%s/%s%s%s" % (database_dir, proteome_name, separator, file_type)):
            #check for exclusion
            if exclude != None and exclude in proteome_file: continue
            #iterate through strands
            for strand in strands:
                #run!
                try:
                    #match list includes the top per_proteome_sequences alignments from each plasmid/genome/strand
                    match_list += directed_local_alignment(proteome_file,
                                           ref_seq_file,
                                           proteome_name,
                                           workdir, mode,
                                           source, strand,
                                           per_proteome_sequences)
                except NoMatchForSeqException:
                    continue
        #sort list
        try:
            #sort match_list by supermatcher/blast score in descending order
            match_list.sort(key=lambda x: x[0], reverse=True)
            #keep the first per_proteome_sequences alignments
            best_aligned_seq_list+=match_list[0:per_proteome_sequences]
        except IndexError:
            #if there are less alignments available than per_proteome_sequences, just keep all the alignments
            best_aligned_seq_list += match_list
    #return
    return best_aligned_seq_list

def filter_insegnificant_seqs(aligned_seq_list, identity_threshold, eval_threshold, gene_percentage):
    """
    This function filters alignments with insufficient identity percentage as determined by the identity_threshold. evalue and length percentage.
    FUTURE: Add the score/best_possible score measure.
    We calculate lots of stuff about the score, and then don't even use it!
    """
    #initialize results list
    signif_seq_list = []
    #iterate through alignments
    for score, record, identity_percentage, evalue, percentage in aligned_seq_list:
        #if score > median + 2 * std:
        #if score >= lower_bound and score <= upper_bound:
        #see if it passes all thresholds
        if identity_percentage > identity_threshold and evalue < eval_threshold and percentage > gene_percentage: #identity_percentage should be float, if it isn't add float()
            signif_seq_list.append(record)
    #return
    return signif_seq_list #list of sequences which were found to be significant

def run_multiple_sequence_alignment(records, workdir, msa):
    """
    This runs the MSA, user can choose between emma, clustalw (old and busted), clustal omega (recommended for proteins and
    also uses HMM), MUSCLE or MAFFT (recommended for nucleotide data, and MUSCLE should be pretty fast), T-Coffee
    (good for distantly related sequences).
    FUTURE: Add more iterative methods to improve runtime? Add HMMER? HHpred is also quite fast
    """
    #get filename for fasta file
    sequence_list_file = os.path.join(workdir, "msa.fasta")
    #write sequences
    SeqIO.write(records, sequence_list_file, "fasta")
    #prepare filenames for MSA output
    outfile = os.path.join(workdir, "msa.aln")
    treefile  =os.path.join(workdir, "msa.dnd")
    #Prepare command line according to chosen algorithm
    if msa.lower() == "emma": #output is fasta
        print "Aligning by emma"
        cmd = EmmaCommandline(sequence=sequence_list_file, outseq=outfile, dendoutfile=treefile)
    elif msa.lower() == "clustalo" or msa.lower() == "clustal_omega" or msa.lower() == "clustal-omega":
        print "Aligning by Clustal Omega"
        cmd = ClustalOmegaCommandline(infile=sequence_list_file, outfile=outfile, verbose=True, auto=True,
                                        guidetree_out=treefile, outfmt="clu", force=True)
    elif msa.lower() == "t-coffee" or msa.lower() == "t_coffee": #should output tree file automatically
        print "Aligning by T-Coffeee"
        cmd = TCoffeeCommandline(infile=sequence_list_file,output="clustalw", outfile=outfile)
    elif msa.lower() == "muscle":
        print "Aligning by MUSCLE"
        #cmd = MuscleCommandline(input=sequence_list_file, out=outfile, tree2=treefile, clw=True)
        cmd = MuscleCommandline(input=sequence_list_file, out=outfile, tree2=treefile)
    elif msa.lower() == "mafft": #probably gonna save tree as input.tree
        print "Aligning by MAFFT"
        cmd = MafftCommandline(input=sequence_list_file, clustalout=True, treeout=True)
    elif msa.lower() == "clustalw" or msa.lower() == "clustalw2":
        print "Aligning by ClustalW2"
        cmd = ClustalwCommandline("clustalw", infile=sequence_list_file, outfile=outfile, tree=True, newtree=treefile)
    elif msa.lower() == "prank": #output is fasta, tree will be outputted to .dnd file?
        print "Aligning by PRANK"
        cmd = PrankCommandline(d=sequence_list_file, o=outfile, f=8, showtree=True, noxml=True)
    elif msa.lower() == "msaprobs": #doesn't use a guide tree
        print "Aligning by MSAprobs"
        cmd = MSAProbsCommandline(infile=sequence_list_file, outfile=outfile, clustalw=True)
    elif msa.lower() == "probcons":
        print "Aligning by ProbCons"
        cmd = ProbconsCommandline(input=sequence_list_file, clustalw=True)
    elif msa.lower() == "dialign": #phylip tree should be created automatically, names are a mystery?
        print "Aligning by Dialign"
        cmd = DialignCommandline(input=sequence_list_file, cw=True, fn=outfile)
    else:
        raise BaseException("Only Multiple Sequence Alignment algorithms currently supported are emma, clustalo, t_coffee, muscle and mafft")
    #Execute the command
    stdout, stderr = cmd()
    #For algorithms that don't have an option to save ouptut to file, capture the stdout
    if msa.lower() == "mafft" or msa.lower() == "probcons":
        with open(outfile, "w") as handle:
            handle.write(stdout)

def export_aligned_sequences(gene_path, aligned_seq_list):
    """
    Write the alignment results of Supermatcher/BLAST to a tab delimited file and add statistics
    """
    #open results file
    with open(os.path.join(gene_path, "seq_scores.csv"), "wb") as fl:
        #initialize write object
        writer = csv.writer(fl, delimiter="\t")
        #iterate through alignments for gene
        for score, seq, identity_percentage, evalue, percentage in aligned_seq_list:
            #write
            writer.writerow([seq.name, seq.description, evalue, identity_percentage, score, seq.seq, len(seq.seq), percentage])
    #Calculating statistics about the evalue, score and identity percentage of all alignments and writing them to file
    #Initialize score list
    SCORE_INDEX = 0 #Index number for the score in the alignment tuple, which is (score,aligned_sequence,identity_percentage)
    score_list = [entry[SCORE_INDEX] for entry in aligned_seq_list]
    #Initialize identity list
    IDENTITY_INDEX = 2 #Index number for the score in the alignment tuple, which is (score,aligned_sequence,identity_percentage)
    identity_list = [entry[IDENTITY_INDEX] for entry in aligned_seq_list]
    #Initialize evalue list
    EVALUE_INDEX = 3 #Index number for the score in the alignment tuple, which is (score,aligned_sequence,identity_percentage)
    evalue_list = [entry[EVALUE_INDEX] for entry in aligned_seq_list]
    #Initialize percetnage list
    PERCENTAGE_INDEX = 4 #Index number for the score in the alignment tuple, which is (score,aligned_sequence,identity_percentage)
    percentage_list = [entry[PERCENTAGE_INDEX] for entry in aligned_seq_list]
    #if no results
    if score_list == [] or identity_list == []: #there are no successful alignments at all even before filtration
        raise BaseException("No successful alignment was produced for any proteome/genome")
    #Write to statistics file
    with open(os.path.join(gene_path, "seq_scores_statistics.csv"), "wb") as fl:
        #initialize writer
        writer = csv.writer(fl, delimiter="\t")
        #write header
        writer.writerow(["Statistic","Minimum","First_quantile","Median","Third_quantile","maximum","Standard_deviation"])
        #write scores
        writer.writerow(["Scores:",min(score_list),np.percentile(score_list, 25),np.median(score_list),
                 np.percentile(score_list, 75),max(score_list), np.std(score_list)])
        #rite identity percentages
        writer.writerow(["Identity:",min(identity_list),np.percentile(identity_list, 25),
                 np.median(identity_list),np.percentile(identity_list, 75),max(identity_list), np.std(identity_list)])
        #write e-values
        writer.writerow(["E-value:",min(evalue_list),np.percentile(evalue_list, 25),
                 np.median(evalue_list),np.percentile(evalue_list, 75),max(evalue_list), np.std(evalue_list)])
        #write length percentages
        writer.writerow(["Percentage:",min(percentage_list),np.percentile(percentage_list, 25),
                 np.median(percentage_list),np.percentile(percentage_list, 75),max(percentage_list), np.std(percentage_list)])

def open_proteome(proteome_file):
    """
    Open a proteome or genome file whether it's gzipped or not, and output its contents as a list of SeqRecord objects (one
    object for each sequence in the file)
    """
    #initialize
    proteome_seq_list = []
    #determine if gzipped
    if proteome_file[-3:] == ".gz":
        #open
        with gz.open(proteome_file) as fl:
            #iterate
            for prot in SeqIO.parse(fl,"fasta"):
                #add to list
                proteome_seq_list.append(prot)
    else:
        #open
        with open(proteome_file) as fl:
            #iterate
            for prot in SeqIO.parse(fl,"fasta"):
                #add to list
                proteome_seq_list.append(prot)
    #return
    return proteome_seq_list

def run(args):
    """
    Run the conservation profiler pipeline
    """
    #initialize from arguements
    proteome_seq_list, proteome_name_list, workdir, database_dir, file_type, exclude, mode, identity_threshold, eval_threshold, gene_percentage, source, per_proteome_sequences, msa, name, data = args
    #create results path for sequence
    gene_path = os.path.join(workdir, name)
    mkdir(gene_path)
    #run everything
    try:
        #Get the reference sequence
        ref_seq = extract_seq(data, proteome_seq_list)
        #Create the reference seq file
        ref_seq_file = make_fasta_per_gene(gene_path, name, ref_seq)
        #Extract aligned sequences from other proteomes
        aligned_seq_list = extract_aligned_seqs(proteome_name_list, 
                                ref_seq_file, 
                                database_dir, 
                                gene_path,
                                #subfolders,
                                file_type,
                                exclude,
                                mode,
                                source,
                                per_proteome_sequences)
        #Save the sequences and their scores into a file
        export_aligned_sequences(gene_path, aligned_seq_list)
        #Use only significant scores, remove ones that don't pass the threshold
        signif_seq_list = filter_insegnificant_seqs(aligned_seq_list, identity_threshold, eval_threshold, gene_percentage)
        #create MSA
        run_multiple_sequence_alignment(signif_seq_list, gene_path, msa)
        #END OF THREAD
        print "finished processing: '%s'" % name

    except ApplicationError as ex:
            print "[Error] failed to run multiple alignment in dir:", workdir
            print ex
            return

    except BaseException as ex:
        print "[Error] General error with processing %s" % name
        print ex
        traceback.print_exc()


def rate_genes(proteome_seq_list,
               gene_data_dict, 
               proteome_name_list,
               workdir, 
               database_dir,
               #subfolders,
               file_type,
               exclude,
               mode,
               identity_threshold,
               eval_threshold,
               gene_percentage,
               source,
               per_proteome_sequences,
               msa):
    """
    Create processes running the alignment of the different sequences.
    """
    #initialize
    pool = Pool(processes=4)
    #get names and data of sequences
    name_list, data_list = zip(*[(name, data) for name, data in gene_data_dict.items()])
    #arguements for processes
    args = zip(repeat(proteome_seq_list, len(gene_data_dict)),
               repeat(proteome_name_list, len(gene_data_dict)),
               repeat(workdir, len(gene_data_dict)),
               repeat(database_dir, len(gene_data_dict)),
               #repeat(subfolders, len(gene_data_dict)),
               repeat(file_type, len(gene_data_dict)),
               repeat(exclude, len(gene_data_dict)),
               repeat(mode, len(gene_data_dict)),
               repeat(identity_threshold, len(gene_data_dict)),
               repeat(eval_threshold, len(gene_data_dict)),
               repeat(gene_percentage, len(gene_data_dict)),
               repeat(source, len(gene_data_dict)),
               repeat(per_proteome_sequences, len(gene_data_dict)),
               repeat(msa, len(gene_data_dict)),
               name_list,
               data_list)
    #run
    results = pool.map(run, args)

def main(argv=None):
    #Parse command line
    try:
        settings, args = process_command_line(argv)
    except BaseException as ex:
        print ex
        return 1
    #Create directory
    mkdir(settings.workdir)
    #Get the target proteomes/genomes name list
    proteome_name_list = extract_proteome_list(settings.proteome_list)
    #Extract reference proteome/genome
    proteome_seq_list = []
    if os.path.isfile(settings.proteome_file): #check if the reference proteome is one file
        proteome_seq_list += open_proteome(settings.proteome_file) #if so, get all its fasta sequences
    elif os.path.isdir(settings.proteome_file):
        for proteome_file in glob.glob("%s/*%s" % (settings.proteome_file,settings.file_type)):
            proteome_seq_list += open_proteome(proteome_file)
    else:
        raise BaseException("Reference Proteome/Genome must be either a file or a directory with genomes")
    #Extract genes names to positions - data is (strand, start, end)
    gene_data_dict = extract_genes_data(settings.genes_data)
    #run everything
    rate_genes(proteome_seq_list,
               gene_data_dict,
               proteome_name_list,
               settings.workdir,
               settings.database_dir,
               #settings.subfolders,
               settings.file_type,
               settings.exclude,
               settings.mode,
               settings.identity_threshold,
               settings.eval_threshold,
               settings.gene_percentage,
               settings.source,
               settings.per_proteome_sequences,
               settings.msa)
    #THE END
    return 0

if __name__ == "__main__":
    exit(main())
