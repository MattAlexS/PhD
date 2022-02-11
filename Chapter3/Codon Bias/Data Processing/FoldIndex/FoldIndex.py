#translate to DNA to protein
#calculate fold index
#add sliding window capability

import argparse

parser = argparse.ArgumentParser()
parser.add_argument("infile", help="Name of input file in FASTA format")
parser.add_argument("outfile", help="Name of output file, format is .csv")
parser.add_argument("-d", "--dna", help="Read DNA. Default: Reads protein sequence", action = "store_true")
parser.add_argument("-r", "--rna", help="Read RNA. Default: Reads protein sequence", action = "store_true")
parser.add_argument("-w", "--window", help="Desired Window size for sliding scale perent folded", type=int)
args = parser.parse_args()

aminos = ["I", "V", "L", "F", "C", "M", "A", "G", "T", "W", "S", "Y", "P", "H", "E", "Q", "D", "N", "K", "R"]

sample = "MVIIIRRRRRRRRRRCVQDRIMCAMIRRRRRRISHTHQSSQSIRRRRRRRRRRRRIIIVMHAISSACTILHCYDRIK"

def translate(dna):
    translation = {
        'AAA':'K','AAC':'N','AAG':'K','AAT':'N','ACA':'T','ACC':'T','ACG':'T','ACT':'T',
        'AGA':'R','AGC':'S','AGG':'R','AGT':'S','ATA':'I','ATC':'I','NTG':'M','ATG':'M',
        'ATT':'I','CAA':'Q','CAC':'H','CAG':'Q','CAT':'H','CCA':'P','CCC':'P','CCG':'P',
        'CCT':'P','CGA':'R','CGC':'R','CGG':'R','CGT':'R','CTA':'L','CTC':'L','CTG':'L',
        'CTT':'L','GAA':'E','GAC':'D','GAG':'E','GAT':'D','GCA':'A','GCC':'A','GCG':'A',
        'GCT':'A','GGA':'G','GGC':'G','GGG':'G','GGT':'G','GTA':'V','GTC':'V','GTG':'V',
        'GTT':'V','TAC':'Y','TAT':'Y','TCA':'S','TCC':'S','TCG':'S','TCT':'S','TGC':'C',
        'TGG':'W','TGT':'C','TTA':'L','TTC':'F','TTG':'L','TTT':'F'
        }
    peptide = []
    for base in range(len(str(dna))):
        if base % 3 == 0:               
            if dna[base:base + 3] in translation.keys():
                peptide.append(translation[dna[base:base+3]])
    protein = ''.join(peptide)
    return protein


def foldIndex(sequence):
    hydrophobicity = 0
    charge = 0
    kyte_doo = {
        "I":1.0000000000000000,
	"V":0.9666666666666667,
	"L":0.9222222222222222,
	"F":0.8111111111111111,
	"C":0.7777777777777778,
    	"M":0.7111111111111111,
	"A":0.7000000000000000,
    	"G":0.4555555555555556,
    	"T":0.4222222222222222,
    	"W":0.4000000000000000,
    	"S":0.4111111111111111,
    	"Y":0.3555555555555556,
    	"P":0.3222222222222222,
    	"H":0.1444444444444444,
    	"E":0.1111111111111111,
    	"Q":0.1111111111111111,
    	"D":0.1111111111111111,
    	"N":0.1111111111111111,
    	"K":0.0666666666666667,
    	"R":0.0000000000000000,
    }
    for residue in sequence:
        if residue == 'D':
            charge -= 1
        elif residue == 'E':
            charge -= 1
        elif residue == 'K':
            charge += 1
        elif residue == 'R':
            charge += 1
        hydrophobicity = hydrophobicity + kyte_doo[residue]
    fold = 2.785 * hydrophobicity/len(sequence) - abs(charge)/len(sequence) - 1.151
    return fold

def foldPercent(sequence, winsize):
    hydrophobicity = 0
    charge = 0
    hydrophobicitywin = []
    chargewin = []
    count = 0
    total = 0
    kyte_doo = {
        "I":1.0000000000000000,
	"V":0.9666666666666667,
	"L":0.9222222222222222,
	"F":0.8111111111111111,
	"C":0.7777777777777778,
    	"M":0.7111111111111111,
	"A":0.7000000000000000,
    	"G":0.4555555555555556,
    	"T":0.4222222222222222,
    	"W":0.4000000000000000,
    	"S":0.4111111111111111,
    	"Y":0.3555555555555556,
    	"P":0.3222222222222222,
    	"H":0.1444444444444444,
    	"E":0.1111111111111111,
    	"Q":0.1111111111111111,
    	"D":0.1111111111111111,
    	"N":0.1111111111111111,
    	"K":0.0666666666666667,
    	"R":0.0000000000000000,
    }
    if winsize >= len(sequence):
        for residue in sequence:
            if residue == 'D':
                charge -= 1
            elif residue == 'E':
                charge -= 1
            elif residue == 'K':
                charge += 1
            elif residue == 'R':
                charge += 1
            hydrophobicity = hydrophobicity + kyte_doo[residue]
        if (2.785 * hydrophobicity/len(sequence) - abs(charge)/len(sequence) - 1.151) >= 0:
            return '100'
        else:
            return '0'
    else:
        for i in range(winsize):
            if sequence[i] == 'D':
                chargewin.append(-1)
            elif sequence[i] == 'E':
                chargewin.append(-1)
            elif sequence[i] == 'K':
                chargewin.append(1)
            elif sequence[i] == 'R':
                chargewin.append(1)
            else:
                chargewin.append(0)
            hydrophobicitywin.append(kyte_doo[sequence[i]])
        charge = sum(chargewin)
        hydrophobicity = sum(hydrophobicitywin)
        for i in range(len(sequence) - winsize):
            if sequence[winsize + i] == 'D':
                chargewin.append(-1)
                charge -= 1
            elif sequence[winsize + i] == 'E':
                chargewin.append(-1)
                charge -= 1
            elif sequence[winsize + i] == 'K':
                chargewin.append(1)
                charge += 1
            elif sequence[winsize + i] == 'R':
                chargewin.append(1)
                charge += 1
            else:
                chargewin.append(0)
            charge -= chargewin.pop(0)
            hydrophobicitywin.append(kyte_doo[sequence[i]])
            hydrophobicity += kyte_doo[sequence[winsize + i]] - hydrophobicitywin.pop(0)
            if (2.785 * hydrophobicity/winsize - abs(charge)/winsize - 1.151) >= 0:
                count += 1
                total += 1
            else:
                total += 1
        return count/total * 100
                 
headers = []
sequences = []
foldindex = []
percentfold = []
winsize = args.window
with open(args.infile, 'r') as file:
    data = file.readlines()
    for line in data:
        line = line.strip()
        str(line)
        if line[0] == '>':
            headers.append(str(line.split(' ')[0][1:]))
        else:
            if args.dna:
                sequences.append(translate(line))
            elif args.rna:
                sequences.append(translate(line.upper().replace('U','T')))
            else:
                sequences.append(line)

for seq in sequences:
    foldindex.append(str(foldIndex(str(seq))))
    percentfold.append(str(foldPercent(str(seq), winsize)))

"""
    if len(seq) < winsize:
        if foldIndex(seq) >= 0:
            percentfold.append('100')
        else:
            percentfold.append('0')
    else:
        count = 0
        total = 0
        position = 0
        while position <= (len(seq) - winsize):
            if foldIndex(seq[position:position+winsize]) >= 0:
                count += 1
                total += 1
                position += 1
            else:
                total += 1
                position += 1
        percentfold.append(str(count/total * 100))
"""

with open(args.outfile, 'w') as file:
    print ("ID,FoldIndex,PrecentFolded", file = file)
    for i in range(len(sequences)):
        print (headers[i] + ',' + foldindex[i] + ',' + percentfold[i], file=file)
        
