from structure import *
from collections import*

def validateseq(dna_seq):
    '''check whether given sequence is a dna or not'''
    temp_seq = dna_seq.upper()
    for nuc in temp_seq:
        if nuc not in nucleotide:
            return False
    return temp_seq

def countNucFreq(seq):
    '''count the nculeotide frequency'''
    tempDNAdict={"A":0, "C": 0, "G": 0, "T": 0}
    for nuc in seq:
        tempDNAdict[nuc]+=1
    return tempDNAdict

def transcription(seq):
    '''Transcription process: DNA--->RNA'''
    return seq.replace("T", "U")

def reverse_transcription(seq):
    '''rersing the trancription process: RNA--->DNA '''
    return seq.replace("U", "T")

def complement_DNA(seq):
    '''Complement of given DNA Seq'''
    #return ''.join([DNA_complement[nuc] for nuc in seq])
    mapping=str.maketrans('ATCG', 'TAGC')
    return seq.translate(mapping)

def reversecomplement_DNA(seq):
    '''Reverse complement of given DNA seq'''
    #return ''.join([DNA_complement[nuc] for nuc in seq])[::-1]
    mapping=str.maketrans('ATCG', 'TAGC')
    return seq.translate(mapping)[::-1]

def gc_content(seq, k=5):
    """ To know the percentage of G & C in given DNA/RNA seq """
    return round(((((seq.count('C'))+seq.count('G'))/ len(seq))*100),k)

def gc_content_subseq(seq, k=20):
     """ To know the percentage of G & C for length (default=20) in given DNA/RNA seq"""
     res=[]
     for i in range(0, len(seq)-k + 1, k):
        subseq= seq[i:i+k]
        res.append(gc_content(subseq))
     return res

def translation(seq, init_pos= 0):
    '''Translate DNA Sequence into Amino acid sequence'''
    return [DNA_Codons[seq[pos:pos + 3]]for pos in range(init_pos, len(seq)-2,3)]

def translation_rna(seq, init_pos= 0):
    ''''Translate RNA Sequence into Amino Acid Sequence'''
    return [RNA_Codons[seq[pos:pos + 3]]for pos in range(init_pos, len(seq)-2,3)]

def codon_freq(seq, aminoacid):
    '''Provide the frequency of each codon representating the particular amino acid'''
    templist= []
    for i in range(0,len(seq)-2,3):
        if DNA_Codons[seq[i:i+3]]==aminoacid:
            templist.append(seq[i:i+3])

    freqDict= dict(Counter(templist))
    totalweight= sum(freqDict.values())

    for seq in freqDict:
        freqDict[seq]= round(freqDict[seq] / totalweight,3)
    return freqDict

def reading_frames(seq):
    '''Generate the six reading frames of DNA sequence, including Reverse Complement'''
    frames=[]
    frames.append(translation(seq,0))
    frames.append(translation(seq,1))
    frames.append(translation(seq,2))
    frames.append(translation(reversecomplement_DNA(seq),0))
    frames.append(translation(reversecomplement_DNA(seq),1))
    frames.append(translation(reversecomplement_DNA(seq),2))
    return frames

def protein_from_reading_frame(aa_seq):
    '''Compute all possible proteins in an aminoacid seq and return a list of possibles'''
    current_prot= []
    proteins= []
    for aa in aa_seq:
        if aa =="_":
            # STOP accumulating amino acid if '_' STOPwas found
            if current_prot:
                for p in current_prot:
                    proteins.append(p)
                current_prot= []
        else:
            #START accumulating amino acid If 'M' START was found
            if aa == 'M':
                current_prot.append("")
            for i in range(len(current_prot)):
                current_prot[i]+= aa
    return proteins

def all_proteins_from_orfs(seq, startReadPos=0, endReadPos=0, ordered=False):
    """Compute all possible proteins for all open reading frames"""
    """Protine Search DB: https://www.ncbi.nlm.nih.gov/nuccore/NM_001185097.2"""
    """API can be used to pull protein info"""
    if endReadPos > startReadPos:
        rfs = reading_frames(seq[startReadPos: endReadPos])
    else:
        rfs = reading_frames(seq)

    res = []
    for rf in rfs:
        prots = protein_from_reading_frame(rf)
        for p in prots:
           res.append(p)

    if ordered:
        return sorted(res, key=len, reverse=True)
    return res

