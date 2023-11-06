import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# read csv

df = pd.read_csv('secreted_dles_per_tax.csv')

# select top 20
df = df.iloc[:20]

# Plot the top 20 species as a bar plot with seaborn
sns.set(style="whitegrid")
sns.set_color_codes("pastel")
sns.barplot(x="Lineage", y="Mean_DLE_freq", data=df,
            label="DLEs", color="b")

# Add a legend and informative axis label
plt.xlabel("Taxonomic Lineage")
plt.ylabel("Mean DLE Frequency")
plt.title("Mean DLE Frequency per Taxonomic Lineage")
plt.xticks(rotation=90)
plt.tight_layout()
plt.savefig('species_figure.png')

