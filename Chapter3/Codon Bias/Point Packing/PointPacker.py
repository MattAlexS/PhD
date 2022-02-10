##Evolutionary Point Packer

import numpy as np
from HyperPlaneSampler import hyper_uniform
import random
import os

aminos = {'K':2,'N':2,'T':4,'R':6,'S':6,'I':3,'Q':2,'H':2,'P':4,'L':6,'E':2,'D':2,'A':4,'G':4,'V':4,'Stop':3,'Y':2,'C':2,'F':2}
aminosNonStop = {'K':2,'N':2,'T':4,'R':6,'S':6,'I':3,'Q':2,'H':2,'P':4,'L':6,'E':2,'D':2,'A':4,'G':4,'V':4,'Y':2,'C':2,'F':2}

codons = ['K','N','K','N','T','T','T','T','R','S','R','S','I','I','I','Q','H','Q','H','P','P','P','P','R','R','R','R','L','L','L','L','E','D','E','D','A','A','A','A','G','G','G','G','V','V','V','V','Stop','Y','Stop','Y','S','S','S','S','Stop','C','C','L','F','L','F']
codonsNonStop = ['K','N','K','N','T','T','T','T','R','S','R','S','I','I','I','Q','H','Q','H','P','P','P','P','R','R','R','R','L','L','L','L','E','D','E','D','A','A','A','A','G','G','G','G','V','V','V','V','Y','Y','S','S','S','S','C','C','L','F','L','F']

#Parameters
ecosize = 100       #maximum population Size
startsize = 30      #number of random points used to initialize a packing
time = 150          #number of generations
mpg = 50            #number of mating events each generation
rmr = 5             #mutation rate
distList = [4.5,4.0,3.75,3.5,3.25,3.00,2.75,2.5,2.25,2.00] #List of minimum distances to use, creates a new packing and file for each.
random.seed(611)

#The following lines can be used to toggle the domain to be packed
#Ensure dimension is defined once
dimension = [aminosNonStop, codonsNonStop]  #Hypersimplicial complex without STOP codons
#dimension = [amino,codons]                 #Hypersimplicial complex including STOP codons
#dimension = 6                              #Hyperplane in given number of dimensions (for single amino acid packings)
#dimension = [4,-4.0,4.0]                   #Generic Space (4 dimensions from -4.0 to 4.0)

#rand_gen is the point generation function for generating a random point
#Ensure only one of these function variants is active.
#Use this in conjunction with the dimension to correctly specify the space you wish to pack
#Hyperplane toggle
"""
def rand_gen(dimension):
    return np.asarray(hyper_uniform(dimension))
"""

#Hypersimplicial Complex toggle
def rand_gen(dimension):
    product = []
    holding = {}
    planes = dimension[0]
    layout = dimension[1]
    for i in planes.keys():
        holding[i] = hyper_uniform(planes[i])
    for i in layout:
        product.append(holding[i].pop())
    return np.asarray(product)

#Generic Space toggle
"""
def rand_gen(dimension):
    point = []
    for i in range(dimension[0]):
        point.append(random.uniform(dimension[1],dimension[2]))
    return np.asarray(point)
"""

#Euclidean Distance
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
        

os.chdir("/Users/matthew/Documents/Codon Usage Project/Clean Attempt/Point Packings/")
for d in distList:
    print(d)
    filename =  "CodonNonStopAnchorsMinD" + str(d)
    #initialize
    ecosystem = []
    for i in range(ecosize):
        mom = []
        dad = []
        for x in range(startsize):
            mom.append(rand_gen(dimension))
            dad.append(rand_gen(dimension))
        ecosystem.append(ConwayCross(d, mom, dad))

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
            newgen.append(ConwayCross(d, mom, dad, mutation))
        for child in newgen:
            for individual in range(len(ecosystem)):
                if len(child) > len(ecosystem[individual]):
                    ecosystem[individual] = child
                    break

    print(len(ecosystem[0]))
    filename = filename + "("  + str(len(ecosystem[0])) + ").csv"

    #write to file
    with open(filename, "w") as file:
        for point in ecosystem[0]:
            temp = []
            for dim in point:
                temp.append(str(dim))
            print(",".join(temp), file = file)
