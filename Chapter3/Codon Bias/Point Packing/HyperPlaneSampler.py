## Willms sampler function

import numpy as np

def hyper_uniform(dimension):
    point = []
    plane_sum = 1.0
    n = dimension - 1
    for i in range(n):
        roll = np.random.ranf() # np.random.randint(2) #for random 0s and 1s
        value = plane_sum*(1-(1-roll)**(1/(n-i)))
        plane_sum -= value
        point.append(value)
    point.append(plane_sum)
    np.random.shuffle(point)
    return(point)

##Generate random points and save to a file
"""
with open("3DhpSample.csv", "w") as file:
    for i in range(100000):
        point = hyper_uniform(3)
        line = ','.join(point)
        print(line, file = file)

"""       
