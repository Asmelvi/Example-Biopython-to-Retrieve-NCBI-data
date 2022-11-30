from Bio import SeqIO
from Bio import Entrez
from Bio.Align import AlignInfo
from Bio import AlignIO
from Bio import Align

# PART 1

def get20protfasta(seqname):

    Entrez.email = "ans@gmail.com"
    ids = Entrez.esearch(db="protein", term=seqname, retmax="20") #Although retmax is 20 by default
    record = Entrez.read(ids)
    id_list = record["IdList"]
    handle = Entrez.efetch(db="protein", id=id_list, rettype="fasta")
    fastas = handle.read()
    return(fastas)

prot_fastas = get20protfasta("mdm2[gene]")
print(prot_fastas)

txt_file_with_fastas = open("WHERE YOU WANT TO SAVE YOUR DATA - YOU ARE TOLD TO SAVE IT AS SEARCHRESULTS.FASTA", "w")
txt_file_with_fastas.write(prot_fastas)
txt_file_with_fastas.close()


# PART 2

def get_consensus_seq(file):
    alignment = AlignIO.read(file, "clustal")
    alignment_summary = AlignInfo.SummaryInfo(alignment)
    consensus = alignment_summary.dumb_consensus(threshold=0)
    return(consensus)

filename = "PATH OF THE CLUSTALW ALN OR THE CLUSTALW ACCURATE ALN"
consensusseq = get_consensus_seq(filename)
print(consensusseq)


# PART 3

#For this one we are getting the functions we already created. There is no need since we can just use the files we have, but since the assignment is about doing 1 step after another, we can simply run it all.

def get20protfasta(seqname):

    Entrez.email = "ans@gmail.com"
    ids = Entrez.esearch(db="protein", term=seqname, retmax="20") #Although retmax is 20 by default
    record = Entrez.read(ids)
    id_list = record["IdList"]
    handle = Entrez.efetch(db="protein", id=id_list, rettype="fasta")
    fastas = handle.read()
    return(fastas)

def get_consensus_seq(file):
    alignment = AlignIO.read(file, "clustal")
    alignment_summary = AlignInfo.SummaryInfo(alignment)
    consensus = alignment_summary.dumb_consensus(threshold=0)
    return(consensus)


#First of all, we get the sequences from those 20 proteins and save that into a fasta file (prot_fastas.txt). We have already done this step before so there is no need, you can simply use the next one and delete this part.

#prot_fastas = get20protfasta("mdm2[gene]")
#txt_file_with_fastas = open("C:/Users/Avalon/Desktop/prot_fastas.txt", "w")
#txt_file_with_fastas.write(prot_fastas)
#txt_file_with_fastas.close()

#Now we open "prot_fastas.txt" and we parse it to get the info.

protein_sequences_all_info = SeqIO.parse("PATH TO THE PLACE WHERE YOU SAVED YOUR SEARCHRESULTS.FASTA", "fasta")

"""for protein_seq in protein_sequences_all_info:
    print(protein_seq.seq)"""                          #You can print this to see only the sequences. This is what we will be using.


#Obviously, for this excercise, we also need the consensus sequence we have generated. We just have to do what we did before. //// Another option would be "copypasting" the consensus seq we got in the Part II.
#--> In this part we DO need the clustalw.aln file we generated before.

filename = "PATH WHERE YOU SAVED YOUR ALIGNMENT FILE"
consensusseq = get_consensus_seq(filename)


#Now what we need to do is performing a Pairwise Alignment between all samples and the consensus seq. That is: Seq1 vs Consensus; Seq2 vs Consensus; Seq3 vs Consensus... and so on. (We have 20 fasta seqs)
#Then we need to calculate the score.
#And then we need to organize them: First by score number (the higher the better) and then, if 2 sequences have the same score, the order them by their name.

#---First we initialize the aligner:

aligner = Align.PairwiseAligner(match_score=1.0)

#Now we align the sequences 1 by 1 against the consensus and we get the score:
#In this assignment it is asked to get the "Score or Common residues", the "gapped sequences" and the "Description of the sequence or ID in this case"

sequences_data ={}
header =[]

for protein_seq in protein_sequences_all_info:
    seq_score = aligner.align(consensusseq, protein_seq.seq).score
    seq_alignment = aligner.align(consensusseq, protein_seq.seq)[0]
    prot_id = protein_seq.id
    #Now we save the seq and the id into a small list
    header = [seq_score, prot_id]
    #And we save that info in a dictionary with the protein name
    sequences_data[str(seq_alignment)] = header

#Now we sort the sequences by score (and name)
sorted_seqs = dict(sorted(sequences_data.items(), key=lambda kv:kv[1], reverse=True))

#And finally we print the results. Keep in mind that we have saved each alignment as "key" and the "ID and Score" as values.
#And we just need to print the 10 best matches, so:

counter = 1

for seq in sorted_seqs.keys():
    head_title = str(sorted_seqs[seq][1]) + " ---> " + str(sorted_seqs[seq][0])
    alignment_shown = seq
    print(head_title)
    print(alignment_shown)
    if counter >= 10:
        break
    counter+=1

#Note, there is no need to place the consensus seq vs consensus seq at the top of the results. It is obvious that it is going give a 100% match.
