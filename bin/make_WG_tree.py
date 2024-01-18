import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import json
from collections import defaultdict
from os import path
from sys import argv, exit
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree
from scipy.spatial.distance import pdist, squareform
from sklearn.preprocessing import StandardScaler

def make_heatmap(sample_names, distance_matrix):

    fig, ax = plt.subplots(figsize=(16, 16))
    plt.imshow(distance_matrix, cmap='viridis', interpolation='nearest')
    plt.colorbar()

    # Adding labels to the rows and columns
    plt.xticks(np.arange(len(sample_names)), sample_names, rotation=45, ha='right')
    plt.yticks(np.arange(len(sample_names)), sample_names)
    plt.title('Heatmap of wgSNP Distance')

    plt.savefig("./snp_heatmap.png")


def count_snps(matrix_data):

    counts_dict = defaultdict(dict)

    for i in matrix_data.to_string().split("\n")[1:]:
        line = [j for j in i.split(" ") if j != ""]
        counts_dict[line[0]]["num_snps"] = line[1:].count("1")
        counts_dict[line[0]]["unique_snps"] = []

    for name, values in matrix_data.items():
        linetups = []
        for i in values.to_string().split("\n"):
            line = [j for j in i.split(" ") if j != ""]
            linetups.append(line)

        snp_n = [l[1] for l in linetups].count("1")

        if snp_n == 1:
            si  = [l[1] for l in linetups].index("1")
            counts_dict[linetups[si][0]]["unique_snps"].append(name)

    with open("allele_stats.json", "w") as fout:
        json.dump(counts_dict, fout, indent=4)


def pca(distance_matrix, matrix_data):

    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    distance_matrix_standardized = StandardScaler().fit_transform(distance_matrix)

    # Instantiate PCA

    if len(matrix_data.index) < 20:
        n = len(matrix_data.index)
    else:
        n = 20
    pca = PCA(n_components=n)

    # Fit and transform the data
    pca_result = pca.fit_transform(distance_matrix_standardized)

    # Create a new DataFrame with PCA results for the first 2 components
    pca_df = pd.DataFrame(data=pca_result[:, :2], columns=['PC1', 'PC2'], index=matrix_data.index)

    # Plotting
    fig, ax = plt.subplots(figsize=(10, 10))

    # PCA scatter plot
    scatter = ax.scatter(pca_df['PC1'], pca_df['PC2'], edgecolors='k', alpha=0.7)

    # Annotate points with sample names
    for sample, (pc1, pc2) in zip(pca_df.index, pca_result[:, :2]):
        ax.text(pc1, pc2, sample)

    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.set_title('PCA WG Distance')
    ax.grid(True)

    # Add scree plot as a smaller subplot in the corner
    scree_ax = fig.add_axes([0.69, 0.67, 0.2, 0.2])
    explained_variance_ratio = pca.explained_variance_ratio_
    scree_ax.plot(range(1, len(explained_variance_ratio) + 1), np.cumsum(explained_variance_ratio), marker='o', linestyle='--')
    scree_ax.set_xlabel('N. PCs')
    scree_ax.set_ylabel('Cum. Variance')

    scree_ax.axvline(x=2, color='red', linestyle='--')  # Draw a vertical line at x=2

    # Label the intersection point on the y-axis
    y_intersect = np.cumsum(explained_variance_ratio)[1]
    scree_ax.text(2, y_intersect, f'(2,{y_intersect:.2f})', color='red', verticalalignment='bottom')

    scree_ax.grid(True)

    plt.savefig("./pca.png")


def make_tree(Z, header_set, prefix, style="ggplot"):
    """ Plots a hierarchical clustering tree. """

    # plt.tight_layout()
    plt.style.use(style)

    plt.figure("HCD", figsize=[16,10])
    plt.title(f"WG Hierarchical Clustering Dendrogram (ward)")
    plt.ylabel(f"Distance")

    ddata = dendrogram(Z, labels=header_set, leaf_rotation=45)

    ## Add distances at branch roots
    for i, d, c in zip(ddata['icoord'], ddata['dcoord'], ddata['color_list']):
        x = 0.5 * sum(i[1:3])
        y = d[1]
        plt.plot(x, y, 'o', c=c)
        plt.annotate("%.3g" % y, (x, y), xytext=(0, -5),
                     textcoords='offset points',
                     va='top', ha='center')

    plt.savefig(f"./{prefix}.png")

def make_empty_out():
    outfiles = ["NO_SNP_TREE.png", "snp_tree.nwk"]

    for f in outfiles:
        with open(f, 'w') as fout:
            print(f"None", file=fout)
    
    gen_empty_plot("Cannot be generated with <3 samples.", "allele_matrix.png")
    gen_empty_plot("Cannot be generated with <3 samples.", "pca.png")
    gen_empty_plot("Cannot be generated with <3 samples.", "snp_heatmap.png")

def gen_empty_plot(text, out):
        save_path = out
        # Create an empty plot
        fig, ax = plt.subplots()

        # Add the text to the plot
        ax.text(0.5, 0.5, text, ha='center', va='center', fontsize=16, color='red')

        # Remove axis labels and ticks
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel('')
        ax.set_ylabel('')

        # Save the plot as a PNG image
        plt.savefig(save_path, bbox_inches='tight', pad_inches=0)
        return save_path

def main(args):
    prefix = path.basename(args[1].split(".csv")[0])
    matrix_data = pd.read_csv(args[1], index_col=0)

    count_snps(matrix_data)

    ## Sample names
    sample_names = list(matrix_data.index)

    if len(sample_names) < 3:
        print(f"Fewer than 3 samples present in the allele matrix. Exiting...")
        make_empty_out()
        exit(0)

    ## Create a distance matrix based on the allele presence/absence matrix
    distance_matrix = np.zeros((len(sample_names), len(sample_names)))

    for i, (index_i, row_i) in enumerate(matrix_data.iterrows()):
        for j, (index_j, row_j) in enumerate(matrix_data.iterrows()):
            d = pdist([list(row_i), list(row_j)])
            distance_matrix[i, j] = d

    make_heatmap(sample_names, distance_matrix)

    pca(distance_matrix, matrix_data)

    z = linkage(squareform(distance_matrix, force='tovector'), 'ward')

    tree = to_tree(z, False)

    make_tree(z, sample_names, prefix)

main(argv)
