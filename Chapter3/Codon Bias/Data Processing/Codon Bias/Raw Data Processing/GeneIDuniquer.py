biasData = {}
with open('Human101codon_bias.csv', 'r') as file:
    data = file.readlines()
    for i in data:
        temp = i.strip().split(',')
        biasData[temp[0]] = temp[1:]

with open('Human101.UNIQ.codon_bias.csv', 'w') as file:
    for i in biasData.keys():
        print(i + ',' + ','.join(biasData[i]), file = file)
