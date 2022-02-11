import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

fold = pd.read_csv("HumanFolds.csv", header = None)

plt.hist(fold[1], 40)


plt.xlabel('FoldIndex')
plt.ylabel('Frequency')
plt.title('Histogram of FoldIndex Values for Human Coding Sequences')
#plt.text(60, .025, r'$\mu=100,\ \sigma=15$')
#plt.xlim(-1.5, 1.5)
plt.ylim(0, 6000)
plt.grid(True)
plt.show()
