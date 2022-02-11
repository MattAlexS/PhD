import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

fold = pd.read_csv("HumanFolds.csv", header = None)

plt.hist(fold[1], 40)


plt.xlabel('FoldIndex')
plt.ylabel('Frequency')
plt.title('Histogram of FoldIndex Values for Human Coding Sequences')
#plt.yscale('log'). #Toggle Logarithmic Y axis
plt.grid(True)
plt.show()
