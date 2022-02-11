#Extracting codon usage information

with open("Human101one_line.fa", "r") as file:
    data = file.readlines()

#gene_sym = []
#gene_ID = []
uniprotID = []
seqs =[]
codon_freq = []
alphabet = ['A','C','G','T']



for i, line in enumerate(data):
    temp = line.strip()
    if temp[0] == '>':
        #gene_sym.append(temp.split(' ')[6].split(':')[1])
        #gene_ID.append(temp.split(' ')[3].split(':')[1])
        uniprotID.append(temp[1:])
        #chrom.append(temp.split(' ')[2].split(':')[2])
        #strand.append(temp.split(' ')[2].split(':')[5])
    else:
        lencheck = len(temp) % 3
        if lencheck != 0:
            extend = (3 - lencheck)
            temp = temp + "X"*extend
        unknown = 0
        total = 0
        seqs.append(temp)
        codons = {
            'AAA':0,
            'AAC':0,
            'AAG':0,
            'AAT':0,
            'ACA':0,
            'ACC':0,
            'ACG':0,
            'ACT':0,
            'AGA':0,
            'AGC':0,
            'AGG':0,
            'AGT':0,
            'ATA':0,
            'ATC':0,
            'ATG':0,
            'ATT':0,
            'CAA':0,
            'CAC':0,
            'CAG':0,
            'CAT':0,
            'CCA':0,
            'CCC':0,
            'CCG':0,
            'CCT':0,
            'CGA':0,
            'CGC':0,
            'CGG':0,
            'CGT':0,
            'CTA':0,
            'CTC':0,
            'CTG':0,
            'CTT':0,
            'GAA':0,
            'GAC':0,
            'GAG':0,
            'GAT':0,
            'GCA':0,
            'GCC':0,
            'GCG':0,
            'GCT':0,
            'GGA':0,
            'GGC':0,
            'GGG':0,
            'GGT':0,
            'GTA':0,
            'GTC':0,
            'GTG':0,
            'GTT':0,
            'TAA':0,
            'TAC':0,
            'TAG':0,
            'TAT':0,
            'TCA':0,
            'TCC':0,
            'TCG':0,
            'TCT':0,
            'TGA':0,
            'TGC':0,
            'TGG':0,
            'TGT':0,
            'TTA':0,
            'TTC':0,
            'TTG':0,
            'TTT':0
            }
        cursor = 0
        while cursor <= len(temp)-2:
            if temp[cursor] in alphabet and temp[cursor + 1] in alphabet and temp[cursor + 2] in alphabet:
                codons[temp[cursor:cursor + 3]] += 1
                total += 1
                cursor += 3
            else:
                unknown += 1
                cursor += 3
        aminos = list(codons.keys())
        aminos.sort()
        usage = []
        for entry in aminos:
            usage.append(str(codons[entry]))
        usage.append(str(unknown))
        usage.append(str(total))
        usage = ','.join(usage)
        codon_freq.append(usage)
"""
with open("Human101.codon_freq.csv", 'w') as file:
    print("ID,Sequence,AAA(K),AAC(N),AAG(K),AAT(N),ACA(T),ACC(T),ACG(T),ACT(T),AGA(R),AGC(S),AGG(R),AGT(S),ATA(I),ATC(I),ATG(M),ATT(I),CAA(Q),CAC(H),CAG(Q),CAT(H),CCA(P),CCC(P),CCG(P),CCT(P),CGA(R),CGC(R),CGG(R),CGT(R),CTA(L),CTC(L),CTG(L),CTT(L),GAA(E),GAC(D),GAG(E),GAT(D),GCA(A),GCC(A),GCG(A),GCT(A),GGA(G),GGC(G),GGG(G),GGT(G),GTA(V),GTC(V),GTG(V),GTT(V),TAA(Stop),TAC(Y),TAG(Stop),TAT(Y),TCA(S),TCC(S),TCG(S),TCT(S),TGA(Stop),TGC(C),TGG(W),TGT(C),TTA(L),TTC(F),TTG(L),TTT(F),Unknown(J),Total", file = file)

    for i in range(len(uniprotID)):
        #print(gene_sym[i] + ',' + gene_ID[i] + ',' + chrom[i] + ',' + strand[i] + ',' + codon_freq[i])
        print(uniprotID[i] + ',' + 'ATG' + ',' + codon_freq[i], file = file)
"""









          





          
