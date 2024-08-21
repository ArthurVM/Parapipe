import sys
import io
from os import path
import pandas as pd
import numpy as np
from collections import defaultdict
import json

import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree
from scipy.spatial.distance import squareform

from Ardal import Ardal


def make_empty_out():
    outfiles = ["NO_SNP_TREE.png", "snp_tree.nwk"]

    for f in outfiles:
        with open(f, 'w') as fout:
            print(f"None", file=fout)
    
    gen_empty_plot("Cannot be generated with <3 samples.", "allele_matrix_tree.png")
    gen_empty_plot("Cannot be generated with <3 samples.", "allele_matrix_pca.png")
    gen_empty_plot("Cannot be generated with <3 samples.", "snp_heatmap.png")


def gen_empty_plot(text, out):
    """ Generate an empty plot showing some text.
    """
    save_path = out
    fig, ax = plt.subplots()

    ax.text(0.5, 0.5, text, ha='center', va='center', fontsize=16, color='red')

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')

    plt.savefig(save_path, bbox_inches='tight', pad_inches=0)
    return save_path


def getNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick


def dist2newick(data):

    sample_names = list(data.index)
    # Create a list of tuples for linkage information
    z = linkage(squareform(data, force='tovector'), 'ward')
    tree = to_tree(z, False)
    nwk = getNewick(tree, "", tree.dist, sample_names)

    return nwk


def pca(distance_matrix, prefix):

    distance_matrix_standardized = StandardScaler().fit_transform(distance_matrix)

    if len(distance_matrix.index) < 20:
        n = len(distance_matrix.index)
    else:
        n = 20
    pca = PCA(n_components=n)

    pca_result = pca.fit_transform(distance_matrix_standardized)

    pca_df = pd.DataFrame(data=pca_result[:, :2], columns=['PC1', 'PC2'], index=distance_matrix.index)

    fig, ax = plt.subplots(figsize=(10, 10))    
    scatter = ax.scatter(pca_df['PC1'], pca_df['PC2'], edgecolors=None, s=100, alpha=0.7, color="grey")

    for sample, (pc1, pc2) in zip(pca_df.index, pca_result[:, :2]):
        ax.text(pc1, pc2, sample)

    ax.set_xlabel('PC1')
    ax.set_ylabel('PC2')
    ax.grid(True)

    ## Add scree plot
    scree_ax = fig.add_axes([0.69, 0.67, 0.2, 0.2])
    explained_variance_ratio = pca.explained_variance_ratio_
    scree_ax.plot(range(1, len(explained_variance_ratio) + 1), np.cumsum(explained_variance_ratio), marker='o', linestyle='--')
    scree_ax.set_xlabel('N. PCs')
    scree_ax.set_ylabel('Cum. Variance')

    scree_ax.axvline(x=2, color='red', linestyle='--')  # Draw a vertical line at x=2

    ## label intersection
    y_intersect = np.cumsum(explained_variance_ratio)[1]
    scree_ax.text(2, y_intersect, f'(2,{y_intersect:.2f})', color='red', verticalalignment='bottom')

    scree_ax.grid(True)

    plt.savefig(f"{prefix}_pca.png", dpi=100)
    # plt.show()


def ete3PrettyTree(nwk, prefix):
    from ete3 import PhyloTree, TreeStyle, NodeStyle, TextFace, AttrFace, add_face_to_node, TreeFace, RectFace
                   
    tree = PhyloTree(nwk)
    tree.ladderize()

    ts = TreeStyle()
    ts.show_branch_length = False
    ts.show_branch_support = False
    ts.show_leaf_name = False

    # Show the tree with the defined style
    tree.render(f"wgSNP_tree.png", dpi=1000, w=1000, tree_style=ts)

    # tree.show(tree_style=ts)


def basicTree(nwk, prefix):
    from Bio import Phylo
    import matplotlib.pyplot as plt

    plt.style.use('ggplot')
    
    tree = Phylo.read(io.StringIO(nwk), "newick")
    tree.ladderize()

    ## Set figure size for A4 dimensions (210mm x 297mm)
    fig, ax = plt.subplots(figsize=(8.27, 11.69))  # A4 size in inches (1 inch = 25.4mm)

    Phylo.draw(tree, do_show=False, axes=ax)

    ax.axis('off')  # Turn off axes

    plt.savefig(f"{prefix}_tree.png", format='png', dpi=100)  # Adjust dpi as needed
    plt.close()


def makeHeatmap(distance_matrix):

    sample_names = distance_matrix.index

    fig, ax = plt.subplots(figsize=(16, 16))
    plt.imshow(distance_matrix, cmap='viridis', interpolation='nearest')
    plt.colorbar()

    # Adding labels to the rows and columns
    plt.xticks(np.arange(len(sample_names)), sample_names, rotation=45, ha='right')
    plt.yticks(np.arange(len(sample_names)), sample_names)
    plt.title('Heatmap of wgSNP Distance')

    plt.savefig("./snp_heatmap.png")


def countSNPs(ard_mat, sample_ids):

    counts_dict = defaultdict(dict)
    allale_mapping = ard_mat.toDict()

    for id in sample_ids:
        unique_alleles = ard_mat.unique([id])
        total_alleles = allale_mapping[id]
        counts_dict[id]["total_snps"] = len(total_alleles)
        counts_dict[id]["unique_snps"] = list(unique_alleles)

    with open("allele_stats.json", 'w') as fout:
        json.dump(counts_dict, fout, indent=4)


def main(args):
    prefix = path.basename(args[1].split(".csv")[0])
    allele_csv = pd.read_csv(args[1], index_col=0)

    ## make ardal object
    ard_mat = Ardal("./allele_matrix.csv")
    wg_D = ard_mat.pairwise(metric="absolute", chunk_size=50)
    wg_D.to_csv(f"absolute_distance.csv")

    print("pairwise finished")

    ## make plots
    ## check there are enough samples to plot
    if len(wg_D.index) < 3:
        make_empty_out()
    else:
        wg_nwk = dist2newick(wg_D)
        basicTree(wg_nwk, prefix)
        pca(wg_D, prefix)
        makeHeatmap(wg_D)

    ## produce snp report JSON
    countSNPs(ard_mat, wg_D.index)


main(sys.argv)
