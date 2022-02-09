##Evolutionary point packer

import numpy as np
from HyperPlaneSampler import hyper_uniform
import random
import os

aminos = {'K':2,'N':2,'T':4,'R':6,'S':6,'I':3,'Q':2,'H':2,'P':4,'L':6,'E':2,'D':2,'A':4,'G':4,'V':4,'Stop':3,'Y':2,'C':2,'F':2}
aminosNonStop = {'K':2,'N':2,'T':4,'R':6,'S':6,'I':3,'Q':2,'H':2,'P':4,'L':6,'E':2,'D':2,'A':4,'G':4,'V':4,'Y':2,'C':2,'F':2}

codons = ['K','N','K','N','T','T','T','T','R','S','R','S','I','I','I','Q','H','Q','H','P','P','P','P','R','R','R','R','L','L','L','L','E','D','E','D','A','A','A','A','G','G','G','G','V','V','V','V','Stop','Y','Stop','Y','S','S','S','S','Stop','C','C','L','F','L','F']
codonsNonStop = ['K','N','K','N','T','T','T','T','R','S','R','S','I','I','I','Q','H','Q','H','P','P','P','P','R','R','R','R','L','L','L','L','E','D','E','D','A','A','A','A','G','G','G','G','V','V','V','V','Y','Y','S','S','S','S','C','C','L','F','L','F']

#Parameters
ecosize = 100
startsize = 30
time = 150
mpg = 50
rmr = 5
minDist = [4.5,4.0,3.75,3.5,3.25,3.00,2.75,2.5,2.25,2.00]
random.seed(611)
dimension = [aminosNonStop, codonsNonStop]


#hyperplane toggle
"""
def rand_gen(dimension):
    return(hyper_uniform(dimension))
"""
def dist(p1, p2):
    out = np.linalg.norm(p1-p2)
    return out

def ConwayCross(minDist, *argv):
    master = []
    child = []
    for arg in argv:
        for i in arg:
            master.append(i)
    random.shuffle(master)
    child.append(master.pop(0))
    for i in master:
        eligibility = True
        counter = 0
        while eligibility == True and counter < len(child):
            if dist(i,child[counter]) < float(minDist):
                eligibility = False
            else:
                counter += 1
        if eligibility == True:
            child.append(i)
    return(child)
        

#hypersimplicial complex toggle
def rand_gen(dimension):
    product = []
    holding = {}
    planes = dimension[0]
    layout = dimension[1]
    for i in planes.keys():
        holding[i] = hyper_uniform(planes[i])
    for i in layout:
        product.append(holding[i].pop())
    return(product)


os.chdir("/Users/matthew/Documents/Codon Usage Project/Clean Attempt/Point Packings/")
for d in minDist:
    print(d)
    filename =  "CodonAnchorsMinD" + str(minDist) + ".csv"
    #initialize
    ecosystem = []
    for i in range(ecosize):
        mom = []
        dad = []
        for x in range(startsize):
            mom.append(rand_gen(dimension))
            dad.append(rand_gen(dimension))
        ecosystem.append(ConwayCross(minDist, mom, dad))

    print('Intialized')

    #optimize
    for year in range(time):
        random.shuffle(ecosystem)
        newgen = []
        for mating in range(mpg):
            mom = ecosystem[random.randint(0,ecosize-1)]
            dad = ecosystem[random.randint(0,ecosize-1)]
            mutation = []
            for i in range(rmr):
                mutation.append(rand_gen(dimension))
            newgen.append(ConwayCross(minDist, mom, dad, mutation))
        for child in newgen:
            for individual in range(len(ecosystem)):
                if len(child) > len(ecosystem[individual]):
                    ecosystem[individual] = child
                    break

    print(len(ecosystem[0]))

    #write to file

    with open(filename, "w") as file:
        for i in ecosystem[0]:
            temp = []
            for d in i:
                temp.append(str(d))
            line = ",".join(temp)
            print(line, file = file)
