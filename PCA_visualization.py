import matplotlib

import matplotlib.pyplot as plt

from sklearn import datasets
from sklearn.decomposition import PCA
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis

# https://scikit-learn.org/stable/auto_examples/decomposition/plot_pca_vs_lda.html#sphx-glr-auto-examples-decomposition-plot-pca-vs-lda-py
# https://stackoverflow.com/questions/47006617/finding-max-min-value-of-individual-columns
# https://stackoverflow.com/questions/22408237/named-colors-in-matplotlib

df = pd.read_csv("data/abalone.data", delimiter=" ", header=None)

result_df = pd.read_csv("output.txt", delimiter=" ", header=None)

df.iloc[0]

min_max_collection = result_df.agg(['min', 'max'])


y = result_df[0]

df = df.drop(0, axis=1)
df

pca = PCA(n_components=2)

X_r = pca.fit(df).transform(df)

colors = ['blue', 'green', 'red', 'cyan', 'magenta', 'yellow', 'black', "navy", "turquoise", "darkorange", 
          'aquamarine', 'mediumseagreen', 'lime']
    
for i in range(0, 10):
    for j in range(0, len(colors)):
        colors.append(colors[j])
    
    
K = min_max_collection[0][1] + 1

real_colors = colors[0 : K]

all_clusters = []

for i in range(0, K):
    all_clusters.append(i)


for color, i, target_name in zip(colors, all_clusters, all_clusters):
    plt.scatter(
        X_r[y == i, 0], X_r[y == i, 1], color=color, alpha=0.8, lw=lw, label=target_name
    )
plt.legend(loc="best", shadow=False, scatterpoints=1)
plt.title("PCA of Abalones dataset")
