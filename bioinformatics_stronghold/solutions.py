from structures import Nucleotides, DNA_ReverseComplement, RNA_Codons, Monoisotopic_Masses
from utilities import readFile

# Counting DNA Nucleotides 
def countNuc(seq):
    freqDict = {"A": 0, "C": 0, "G": 0, "T": 0}
    for nuc in seq:
        freqDict[nuc] += 1
    results = " ".join([str(val) for key, val in freqDict.items()])
    return results

# Transcribing DNA into RNA 
def DNA2RNA(seq):
    return seq.replace("T", "U")

# Complementing a Strand of DNA 
def reverse_complement(seq):
    reverseSeq = seq[::-1]
    reverseComplement = "".join([DNA_ReverseComplement[nuc] for nuc in reverseSeq])
    return reverseComplement

# Counting Point Mutations 
def hamming_distance(seq1, seq2):
    count = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            count += 1
    return count

# Computing GC Content
def compute_GC(filePath):
    FASTAFile = readFile(filePath)
    FASTADict = {}
    FASTALabel = ""
    for line in FASTAFile:
        if ">" in line:
            FASTALabel = line
            FASTALabel = FASTALabel.split(">")[1]
            FASTADict[FASTALabel] = ""
        else:
            FASTADict[FASTALabel] += line
    
    RESULTDict = {}
    for key, value in FASTADict.items():
        C_count = value.count("C")
        G_count = value.count("G")
        result = (C_count + G_count) / len(value) * 100
        RESULTDict[key] = round(result, 6)
        
    maxGCKey = max(RESULTDict, key=RESULTDict.get)
    maxGCValue = RESULTDict[maxGCKey]
    return maxGCKey, maxGCValue

# Translating RNA into Protein
def RNA_to_protein(seq):
    k = 3   # length of a single codon
    protein = ""
    for i in range(0, len(seq) - k, k):
        codon = seq[i: i+k]
        aminoacid = RNA_Codons[codon]
        protein += aminoacid
    return protein


# Calculating Protein Mass
def calculate_mass(protein):
    total_mass = 0
    for aminoacid in protein:
        mass = Monoisotopic_Masses[aminoacid]
        total_mass += mass
    return round(total_mass, 3)


# Finding a Motif in DNA
def find_pattern(seq, pattern):
    locations = []
    k = len(pattern)
    for i in range(0, len(seq)-k):
        window = seq[i:i+k]
        if window == pattern:
            locations.append(i+1)
    return locations


# Rabbits and Recurrence Relations
def find_num_rabbits(months, offsprings):
    parent, child = 1, 1
    for iter in range(months - 1):
        child, parent = parent, parent + (child * offsprings)
    return child
