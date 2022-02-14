import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt

packing = pd.read_csv("CodonNonStopAnchorsMinD2.54(33)Header.csv", header = 0)

ax = sns.heatmap(packing, vmin=0, vmax=1, cmap='viridis', xticklabels=True, yticklabels=True)
plt.xlabel("Codon")
plt.ylabel("Anchor Number")
plt.title("Packing Bias Visualization")
plt.show()
