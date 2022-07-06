##Evolutionary Point Packer

import numpy as np
from HyperPlaneSampler import hyper_uniform
import random
import os
import matplotlib.pyplot as plt

#Parameters
ecosize = 100       #maximum population Size
startsize = 2      #number of random points used to initialize a packing
time = 80        #number of generations
mpg = 50            #number of mating events each generation
rmr = 2             #mutation rate

np.random.seed(611)


dimension = 3       #dimension of space to pack
d = 0.1            #minimum distance

#rand_gen is the point generation function for generating a random point
#Use this in conjunction with the dimension to correctly specify the space you wish to pack

#Original generation method

def ori_gen(dimension):
    point = np.random.uniform(0,1,dimension)
    point = point/np.sum(point)
    return point



def sim_gen(dimension):
    return np.asarray(hyper_uniform(dimension))


#Euclidean Distance
def dist(p1, p2):
    out = np.linalg.norm(p1-p2)
    return out

def fittest(population):
    m = 0
    for p in population:
        if len(p) > m:
            m = len(p)
    return m


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
        
#Target file location for point packings

gen_mets = ['original','simplex']

results = {
    'original':[],
    'simplex' :[]
    }

for g in gen_mets:
    if g == 'original':
        rand_gen = ori_gen
    else:
        rand_gen = sim_gen
    #initialize
    ecosystem = []
    for i in range(ecosize):
        mom = []
        dad = []
        for x in range(startsize):
            mom.append(rand_gen(dimension))
            dad.append(rand_gen(dimension))
        ecosystem.append(ConwayCross(d, mom, dad))
    results[g].append(fittest(ecosystem))
    

    print('Intialized')

    #optimize
    for year in range(time-1):
        print(year)
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
        results[g].append(fittest(ecosystem))
        
plt.plot(results['original'], color = 'orange', label = 'original')
plt.plot(results['simplex'], color = 'blue', label = 'simplex')
plt.title('3d Simplex Packing Generation Method Comparison')
plt.ylabel('Largest Packing Size')
plt.xlabel('Generation')
plt.legend()
plt.show()

  
