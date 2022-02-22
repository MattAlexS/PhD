import pandas as pd

#Read in Tissue Gene Expression Data
te = {}

with open('TE.csv', 'r') as file:
    data = file.readlines()
    tissues = data[0].strip().split(',')[1:]
    for line in data[1:]:
        temp = line.strip().split(',')
        exp = []
        gene = temp[0]
        for i in temp[1:]:
            exp.append(float(i))
        te[gene] = exp

elevated = []
enhanced = []
enriched = []

for gene in te.keys():
    temp = te[gene].copy()
    mean = sum(te[gene])/len(te[gene])
    m = max(te[gene])
    temp.remove(m)
    m2 = max(temp)
    for i in range(32):
        if te[gene][i] == m:
            elevated.append([gene, tissues[i]])
        if te[gene][i] > 5*mean:
            enhanced.append([gene, tissues[i]])
        if te[gene][i] > 5*m2:
            enriched.append([gene, tissues[i]])
    
ele = pd.DataFrame(elevated, columns = ["Gene ID", "Tissue"])
enh = pd.DataFrame(enhanced, columns = ["Gene ID", "Tissue"])
enr = pd.DataFrame(enriched, columns = ["Gene ID", "Tissue"])

cb = pd.read_csv("HumanTEcodonbias.csv")

fullCopies = list(cb['Gene ID'])
fullUniq = list(set(cb['Gene ID']))
fullCounts = []
for i in fullUniq:
    fullCounts.append(0)

for gene in fullCopies:
    fullCounts[fullUniq.index(gene)] += 1

eleU = list(set(ele['Gene ID']))
enhU = list(set(enh['Gene ID']))
enrU = list(set(enr['Gene ID']))

eleCounts = []
enhCounts = []
enrCounts = []

for i in range(len(fullUniq)):
    if fullUniq[i] in eleU:
        eleCounts.append(fullCounts[i])
    if fullUniq[i] in enhU:
        enhCounts.append(fullCounts[i])
    if fullUniq[i] in enrU:
        enrCounts.append(fullCounts[i])

print(sum(fullCounts)/len(fullCounts))
print(sum(eleCounts)/len(eleCounts))
print(sum(enhCounts)/len(enhCounts))
print(sum(enrCounts)/len(enrCounts))
    
ele.to_csv("ElevatedGenes.csv", index = False)
enh.to_csv("EnhancedGenes.csv", index = False)
enr.to_csv("EnrichedGenes.csv", index = False)
    

