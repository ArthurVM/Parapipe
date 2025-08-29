import sys
import io
from os import path
import pandas as pd
import numpy as np
from collections import defaultdict
import json
from io import StringIO
from Bio import Phylo
from plotly.offline import download_plotlyjs, init_notebook_mode,  iplot
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.cluster.hierarchy import dendrogram, linkage, to_tree
from scipy.spatial.distance import squareform

from plotting_funcs import *


def make_empty_out():
    outfiles = ["NO_SNP_TREE.png", "snp_tree.nwk"]

    for f in outfiles:
        with open(f, 'w') as fout:
            print(f"None", file=fout)
    
    return gen_empty_plot("Cannot be generated with <3 samples.")


def gen_empty_plot(text):
    """ Generate an empty plot showing some text.
    """
    save_path = out
    fig, ax = plt.subplots()

    ax.text(0.5, 0.5, text, ha='center', va='center', fontsize=16, color='red')

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')

    return fig


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


def basicTree(data_df, prefix):
    from Bio import Phylo
    import matplotlib.pyplot as plt

    nwk = dist2newick(data_df)

    plt.style.use('ggplot')
    
    tree = Phylo.read(io.StringIO(nwk), "newick")
    tree.ladderize()

    ## Set figure size for A4 dimensions (210mm x 297mm)
    fig, ax = plt.subplots(figsize=(8.27, 11.69))  # A4 size in inches (1 inch = 25.4mm)

    Phylo.draw(tree, do_show=False, axes=ax)

    ax.axis('off')  # Turn off axes

    plt.savefig(f"{prefix}_tree.png", format='png', dpi=100)  # Adjust dpi as needed
    plt.close()

    return f"{prefix}_tree.png"


def get_x_coordinates(tree):
    # Associates to  each clade a x-coord.
    # returns a dict {clade: x-coord}, i.e the key is a clade, and x-coord its value
    
    xcoords = tree.depths()
    # tree.depth() maps tree clades to depths (by branch length).
    # returns a dict {clade: depth} where clade runs over all Clade instances of the tree,
    # and depth is the distance from root to clade
    
    # If there are no branch lengths, assign unit branch lengths
    if not max(xcoords.values()):
        xcoords = tree.depths(unit_branch_lengths=True)
    return xcoords


def get_y_coordinates(tree, dist=1.3):
    # y-coordinates are   multiple of dist (i*dist below); 
    # dist: vertical distance between two consecutive leafs; it is chosen such that to get a tree of 
    # reasonable height 
    # returns  a dict {clade: y-coord}
        
    maxheight = tree.count_terminals()#Counts the number of tree leafs.
  
    ycoords = dict((leaf, maxheight - i*dist) for i, leaf in enumerate(reversed(tree.get_terminals())))
    def calc_row(clade):
            for subclade in clade:
                if subclade not in ycoords:
                    calc_row(subclade)
            ycoords[clade] = (ycoords[clade.clades[0]] +
                              ycoords[clade.clades[-1]]) / 2

    if tree.root.clades:
        calc_row(tree.root)
    return ycoords


def get_clade_lines(orientation='horizontal', y_curr=0, x_start=0, x_curr=0, y_bot=0, y_top=0,
                    line_color='rgb(25,25,25)', line_width=0.5):
    # define a Plotly shape of type 'line', for each branch
    
    branch_line = dict(type= 'line',
                       layer='below',
                       line=dict(color=line_color, 
                                 width=line_width)
                     )
    if orientation == 'horizontal':
        branch_line.update(x0=x_start,
                           y0=y_curr,
                           x1=x_curr,
                           y1=y_curr)
    elif orientation == 'vertical':
        branch_line.update(x0=x_curr,
                           y0=y_bot,
                           x1=x_curr,
                           y1=y_top)
    else:
        raise ValueError("Line type can be 'horizontal' or 'vertical'")
       
    return branch_line    
        

def draw_clade(clade, x_start, x_coords, y_coords, line_shapes,  line_color='rgb(15,15,15)', line_width=1):
    # defines recursively  the tree  lines (branches), starting from the argument clade
    
    x_curr = x_coords[clade]
    y_curr = y_coords[clade]
   
    # Draw a horizontal line 
    branch_line=get_clade_lines(orientation='horizontal', y_curr=y_curr, x_start=x_start, x_curr=x_curr,  
                                line_color=line_color, line_width=line_width)
   
    line_shapes.append(branch_line)
   
    if clade.clades:
        # Draw a vertical line connecting all children
        y_top = y_coords[clade.clades[0]]
        y_bot = y_coords[clade.clades[-1]]
       
        line_shapes.append(get_clade_lines(orientation='vertical', x_curr=x_curr, y_bot=y_bot, y_top=y_top,
                                           line_color=line_color, line_width=line_width))
       
        # Draw descendants
        for child in clade:
            draw_clade(child, x_curr, x_coords, y_coords, line_shapes)


def make_itree(dist_df, run_samples):

    ## initialise the tree phylo object
    tree_data = dist2newick(dist_df)
    handle = StringIO(tree_data)
    tree = Phylo.read(handle, 'newick')

    ## get node coordinates
    x_coords = get_x_coordinates(tree)
    y_coords = get_y_coordinates(tree)

    line_shapes = [] 
    draw_clade(tree.root, 0, x_coords, y_coords, line_shapes, line_color='rgb(25,25,25)', line_width=1)

    my_tree_clades = x_coords.keys()
    X = [] # list of nodes x-coordinates
    Y = [] # list of nodes y-coords
    text = [] #list of text to be displayed on hover over nodes

    for cl in my_tree_clades:
        X.append(x_coords[cl])
        Y.append(y_coords[cl])
        text.append(cl.name)

    intermediate_node_color = 'rgb(100,100,100)'

    colour = [intermediate_node_color]*len(X)

    for sample in dist_df.index:
        i=text.index(sample)
        if sample in run_samples:
            text[i]=text[i]+f" run sample"
            colour[i] = 'red'
        else:
            text[i]=text[i]+f" non-run sample"
            colour[i] = 'grey'

    nodes = dict(type='scatter',
             x=X,
             y=Y,
             mode='markers',
             marker=dict(color=colour, 
                         size=10),
             opacity=1.0,
             text=text,
             hoverinfo='text')
    
    layout=dict(title='wgSNP tree',
            font=dict(family='Balto',size=14),
            width=750,
            height=1500,
            autosize=False,
            showlegend=False,
            xaxis=dict(showline=True,  
                       zeroline=False,
                       showgrid=False,
                       ticklen=4,          
                       showticklabels=True,
                       title='branch length'),
            yaxis=dict(visible=False), 
            hovermode='closest',
            plot_bgcolor='rgb(250,250,250)',
            margin=dict(l=10),
            shapes=line_shapes # lines for tree branches
           )
    
    fig = go.FigureWidget(data=[nodes], layout=layout)

    return fig


def make_WG_plots(dist_df, prefix, run_samples):

    ## handle instances where there are fewer than 3 samples
    if len(dist_df.index) < 3:
        print("fewer than 3 samples detected. Skipping phylo plotting.")
        empty_plot = gen_empty_plot("Too few in this run for run tree.", "./empty_plot.png")
        wgsnp_tree_prefix, wgsnp_tree_b64 = embed_image(empty_plot)

        phylo_plots_wrapper = f"""
        <div class="plots-container">
            <div class="plot-wrapper">
                <img src="data:image/{wgsnp_tree_prefix};base64,{wgsnp_tree_b64}">
            </div>
        </div>"""

        return phylo_plots_wrapper


    ## handle instances where there are fewer than 3 run samples (and therefore no basic tree is produced)
    ## but 3 or more samples when merged with the input database
    elif len(run_samples) < 3 and len(dist_df.index) >= 3:
        print("fewer than run 3 samples detected but 3 or more present including database. Skipping run tree.")
        ## create empty basic tree plot
        empty_plot = gen_empty_plot("Too few samples for phylogenetic analysis.", "./empty_plot.png")
        wgsnp_tree_prefix, wgsnp_tree_b64 = embed_image(empty_plot)

    ## handle instances where there are 3 or more run samples
    else:
        print("Sufficient numbers of samples present. Building all plots.")
        run_dist_dif = dist_df[run_samples].T[run_samples]
        tree_path = basicTree(run_dist_dif, prefix)
        wgsnp_tree_prefix, wgsnp_tree_b64 = embed_image(tree_path)
                                                        
    itree_fig = make_itree(dist_df, run_samples)
    pca_fig = create_PCA_3D(dist_df, run_samples, fig_height=1000)
    run_dist_df = dist_df[run_samples].T[run_samples]
    hm_fig = create_distance_heatmap(run_dist_df.to_numpy(), row_labels=list(run_dist_df.index), col_labels=list(run_dist_df.index)[::-1])

    phylo_plots_wrapper = f"""
    <div class="plots-container">
        <div class="plot-wrapper">
            {itree_fig.to_html(full_html=False)}
        </div>
        <div class="plot-wrapper">
            <img src="data:image/{wgsnp_tree_prefix};base64,{wgsnp_tree_b64}">
        </div>
    </div>

    <div class="plots-container">
        <div class="plot-wrapper">
            {pca_fig.to_html(full_html=False)}
        </div>
        <div class="plot-wrapper">
            {hm_fig.to_html(full_html=False)}
        </div>
    </div>"""

    return phylo_plots_wrapper


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
