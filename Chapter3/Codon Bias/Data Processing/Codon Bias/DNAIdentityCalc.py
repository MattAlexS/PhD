fasta = {}
with open('Human_dna.fa', 'r') as file:
    data = file.readlines()
    for i in range(len(data)):
        temp = data[i].strip()
        if temp[0] == ">":
            symbol = temp.split('.')[0]
            fasta[symbol[1:]] = data[i+1].strip()

print("Read in")
cutoff = 75

compares = []
with open('Human_blast.csv', 'r') as file:
    data = file.readlines()
    for i in data:
        temp = i.strip().split(',')
        if float(temp[2]) >= cutoff:                  
            compares.append(temp)

print("Assembled")

output = []
for i in compares:
    first = fasta[i[0]]
    second = fasta[i[1]]
    if first != second:
        protcentage = i[2]
        start1 = (int(i[3]) - 1) * 3
        start2 = (int(i[5]) - 1) * 3
        matches = 0
        total = 0
        buffer1 = 0
        buffer2 = 0
        length = 0
        for x in range(len(i[7])):
            if i[7][x] != "-" and i[8][x] != "-":
                if i[7][x] == i[8][x]:
                    length += 1
                    cursor1 = start1 + (buffer1 + x)*3
                    cursor2 = start2 + (buffer2 + x)*3
                    if first[cursor1:cursor1 + 3] == second[cursor2:cursor2 + 3]:
                        matches += 1
                        total += 1
                    else:
                        total += 1
            else:
                if i[7][x] == "-":
                    buffer2 += 1
                else:
                    buffer1 += 1
        output.append([i[0],i[1],float(protcentage),(matches/total)*100,(float(protcentage)-((matches/total)*100)), length])

print("Assessed")

output.sort(key = lambda  x: x[4])

print("Sorted")

minlength = 50

x=[]
y=[]
final = []

for i in output:
    if i[5] >= minlength:
        x.append(float(i[2]))
        y.append(float(i[3]))
        final.append(i)

#Only allow unique pairs
unique = {}
for i in range(len(output)):
    if output[i][5] >= 50:
        ids = [output[i][0], output[i][1]]
        ids = sorted(ids)
        unique[','.join(ids)] = [output[i][2], output[i][3], output[i][4], output[i][5]]
    
        


        



import numpy as np
import matplotlib.pyplot as plt

plt.scatter(x,y)
#plt.xlim(0,100)
plt.ylim(0,100)
plt.xlabel("Protein Identity (%)")
plt.ylabel("Percentage of Identical Codons at Identical AA (%)")
plt.title("Percentage of Matching Aligned AAs with Matching Codon Used")
plt.show()



