fasta = {}

with open("GRCh38.p13Human101.fa", "r") as file:
    data = file.readlines()
    for line in data:
        temp = line.strip()
        if len(temp) > 0:
            if temp[0] == ">":
                acc = temp
                fasta[acc] = []
            else:
                fasta[acc].append(temp)
        
            

with open("Human101one_line.fa", "w") as file:
    for seq in fasta.keys():
        print(seq, file = file)
        print("".join(fasta[seq]), file = file)
