import os
import io
import sys
import json
import shutil
import numpy as np
import pandas as pd
import calendar
import argparse
from collections import defaultdict
from datetime import datetime
from io import StringIO
from Bio import Phylo
from matplotlib import pyplot as plt
from matplotlib import rcParams
from matplotlib.table import Table
import matplotlib.colors as mcolors
from pdf import PDF

rcParams['axes.spines.top'] = False
rcParams['axes.spines.right'] = False


def parse_json(json_file):

    with open(json_file, 'r') as j:
        report_json = json.loads(j.read())

    return report_json


def pair_jpegs(jpeg_dir):
    """ takes the directory containing coverage plots and pairs according to chromosome
    """
    s_jpegs = [[os.path.join(jpeg_dir, f)] for f in sorted(os.listdir(jpeg_dir))]
    # return zip(s_jpegs[::2], s_jpegs[1::2])
    return s_jpegs


def parse_tree(tree_file):

    with open(tree_file, 'r') as fin:
        return fin.read().strip("\n")


def make_plot_dir():
    PLOT_DIR = 'plots'

    # Delete folder if exists and create it again
    try:
        shutil.rmtree(PLOT_DIR)
        os.mkdir(PLOT_DIR)
    except FileNotFoundError:
        os.mkdir(PLOT_DIR)


def make_gg_plot(gg_dict):
    gg_df = pd.DataFrame([ [int(k), float(v)] for k, v in json_data["mapping"]["mapping_stats"]["gg_array"].items() ])

    plt.style.use('ggplot')
    GG_f, GG_p = plt.subplots(1, figsize=(8, 4), dpi=100)

    window_array = gg_df.loc[:,0].tolist()
    w = window_array[1]-window_array[0]

    coldata = gg_df.loc[:,1].tolist()

    norm_array = [x+1-np.max(coldata) for x in coldata]
    m_area = np.trapz(np.array(coldata), dx=w)/np.max(window_array)
    m_areaN = np.trapz(norm_array, dx=w)/np.max(window_array)

    GG_p.plot(window_array, coldata, linestyle="-", linewidth=1, label=f"GG area={np.round(m_area, 3)}\nNormalised area={np.round(m_areaN, 3)}")

    GG_p.set_xlim([0, np.max(window_array)+1])

    GG_p.set(ylabel="$G$", xlabel="$W$")

    plt.title("Gini Granularity Curve")

    GG_p.legend(ncol=3)
    plt.tight_layout()

    figout = "./plots/gg_curve.jpeg"
    plt.savefig(figout)

    return figout


def basicTree(nwk, format, prefix):

    plt.style.use('ggplot')
    
    tree = Phylo.read(io.StringIO(nwk), format)
    tree.ladderize()

    ## Set figure size for A4 dimensions (210mm x 297mm)
    fig, ax = plt.subplots(figsize=(8.27, 11.69))  # A4 size in inches (1 inch = 25.4mm)

    Phylo.draw(tree, do_show=False, axes=ax)

    ax.axis('off')  # Turn off axes

    plt.savefig(f"{prefix}_tree.png", format='png', dpi=100)  # Adjust dpi as needed
    plt.close()

    return f"{prefix}_tree.png"


def build_tree(alignments, treedata):

    genetree = PhyloTree(treedata)
    # genetree.link_to_alignment(alignments)
    return genetree, TreeStyle()


def make_map_table(json_data):
    plt.style.use('ggplot')

    array_names = ["coverage_breadth_hist", "gg_array"]

    ordered_list = [["ID", "Mean DOC", "Median DOC", "BOC>=5x", "GG area", "Norm GG area", "Av. Qual", "Av. Insert Size", "Het. k", "Bound. Prob", "SNPs: total (unique)"]]

    for sample, report in json_data.items():
        mapping_stats = report["mapping"]["mapping_stats"]
        al_stats = report["mapping"]["allele_report"]

        ## handle instances where no MOI analysis was carried out due to insufficient SNP heterozygosity
        if report["heterozygosity"]["k"] == None:
            het_k = "NA"
            het_bp = "NA"
        else:
            het_k = report["heterozygosity"]["k"]
            het_bp = np.round(report["heterozygosity"]["boundary_prob"], 3)

        ordered_list.append([ sample,\
            np.round(mapping_stats["mean_depth_of_coverage"], 1), \
            np.round(mapping_stats["median_depth_of_coverage"], 1), \
            mapping_stats["coverage_breadth_hist"]["5"], \
            np.round(mapping_stats["GG_area"], 3), \
            np.round(mapping_stats["nGG_area"], 3), \
            np.round(mapping_stats["average_quality"], 1), \
            mapping_stats["insert_size_average"], \
            het_k, \
            het_bp, \
            f"{al_stats['total_snps']} ({len(al_stats['unique_snps'])})"])

    fig, ax = plt.subplots(1,1, figsize=(12, 12), dpi=100)

    fig.patch.set_visible(False)
    ax.axis('off')
    ax.axis('tight')

    tab = ax.table(ordered_list, colLabels=None, loc='center')

    tab.auto_set_font_size(False)
    tab.auto_set_column_width(col=list(range(len(ordered_list))))

    fig.tight_layout()

    figout = "./plots/mapstats.jpeg"
    plt.savefig(figout)

    return figout


def make_st_table(json_data):
    plt.style.use('ggplot')

    figs = {}

    st_stats = defaultdict(dict)
    st_warnings = defaultdict(dict)

    for sample, report in json_data.items():
        for scheme, data in report["typing"].items():
            st_stats[scheme][sample] = data
            # st_warnings[scheme][sample] = report["st_stats"]["warnings"][scheme]

    for scheme, data in st_stats.items():

        loci = [locus for locus, cov in data[next(iter(data))]["locus_coverage"].items()]
        table = [["Sample", *loci]]

        for sample, cov_data in data.items():
            row = [sample, *[f"{np.round(float(locus_cov['BOC']), 3)} ({np.round(float(locus_cov['DOC']), 3)}x)" for locus, locus_cov in cov_data["locus_coverage"].items()]]
            table.append(row)

        fig, ax = plt.subplots(1,1, figsize=(12,12), dpi=100)

        fig.patch.set_visible(False)
        ax.axis('off')
        ax.axis('tight')

        tab = ax.table(table, colLabels=None, loc='center')

        # Iterate over rows and columns to set cell colors based on st_warnings
        for i, row in enumerate(table):
            for j, val in enumerate(row):
                if i > 0 and j > 0:  # Skip the header row and column
                    cell_color = get_cell_color(row[0], table[0][j], st_warnings[scheme])
                    tab[(i, j)].set_facecolor(cell_color)

        tab.auto_set_font_size(False)
        tab.auto_set_column_width(col=list(range(len(table))))

        fig.tight_layout()

        figout = f"./plots/{scheme}_stats.jpeg"
        plt.savefig(figout)

        figs[scheme] = figout

    return figs


def get_cell_color(sample, locus, st_warnings):
    if sample in st_warnings and locus in st_warnings[sample]:
        return mcolors.to_rgba('red', alpha=0.5)  # Light red with alpha (transparency)
    else:
        return 'white'


def gen_empty_plot(text):
    save_path = figout = f"./plots/SNP_phylo.jpeg"
    # Create an empty plot
    fig, ax = plt.subplots()

    # Add the text to the plot
    text = "Coverage too low for inclusion in SNP tree"
    ax.text(0.5, 0.5, text, ha='center', va='center', fontsize=16, color='red')

    # Remove axis labels and ticks
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel('')
    ax.set_ylabel('')

    # Save the plot as a PNG image
    plt.savefig(save_path, bbox_inches='tight', pad_inches=0)
    return save_path
    

def make_plots(json_data):

    make_plot_dir()
    map_tab = make_map_table(json_data)
    
    head_data = json_data[next(iter(json_data))]
    ## check to see if a typing scheme was defined in the first json
    if head_data["typing"] != None:
        st_tab = make_st_table(json_data)
        # gg_curve = make_gg_plot(json_data)

        st_plots = []
        for scheme, typing_data in head_data["typing"].items():
                st_plot = basicTree(typing_data["phylo"]["contree"], "newick", scheme)
                st_plots.append([scheme, st_plot])
    else:
        st_plots = None
        st_tab = None

    return map_tab, st_plots, st_tab


def make_moi_table(moi_stats):
    plt.style.use('ggplot')

    ordered_list = []

    ## manually insert fields into the list to maintain some meaningful order
    ## mapping fields
    ordered_list.append(["k", str(moi_stats["k"])])
    ordered_list.append(["Fit Score", str(moi_stats["fit_score"])])
    ordered_list.append(["Boundary Probability", str(moi_stats["boundary_prob"])])

    # df = pd.DataFrame(ordered_list)

    # hide axes
    fig, ax = plt.subplots(figsize=(8, 8), dpi=100)

    # hide axes
    fig.patch.set_visible(False)
    ax.axis('off')
    ax.axis('tight')

    ax.table(ordered_list, colLabels=None, loc='center')

    fig.tight_layout()

    figout = "./plots/moistats.jpeg"
    plt.savefig(figout)

    return figout


def parse_args(argv):

    parser = argparse.ArgumentParser()

    parser.add_argument('script_path', action='store', help=argparse.SUPPRESS)
    parser.add_argument('--env', action='store', nargs=1, help='JSON containing run environment info.')
    parser.add_argument('--report', action='store', nargs="*", help='JSON reports produced by Parapipe.')
    parser.add_argument('--moi_json', action='store', nargs="*", help='MOI JSON reports.')
    # parser.add_argument('--png', action='store', nargs="*", help='PNG plots.')

    args = parser.parse_args(argv)

    return args


def main(argv):

    args = parse_args(argv)

    env = parse_json(args.env[0])

    reports = {}
    for rep in args.report:
        prefix = os.path.basename(rep).split("_report.json")[0]
        reports[prefix] = parse_json(rep)

    for rep in args.moi_json:
        prefix = os.path.basename(rep).split(".heterozygosity.json")[0]
        reports[prefix]["heterozygosity"] = parse_json(rep)
        
    print(json.dumps(reports, indent=4))

    map_tab, st_plots, st_tab = make_plots(reports)

    pdf = PDF()

    pdf.print_page([map_tab])
    pdf.write_title(f"Parapipe Run")

    pdf.print_page(["allele_matrix_tree.png"])
    pdf.write_subheader("wgSNP tree")

    pdf.print_page(["allele_matrix_pca.png"])
    pdf.write_subheader("wgSNP PCA")

    pdf.print_page(["./snp_heatmap.png"])
    pdf.write_subheader("wgSNP Heatmap")

    if st_plots != None and st_tab != None:
        for scheme, st_plot in st_plots:
            print(scheme)
            pdf.print_page([st_plot], big_small=True)
            pdf.write_subheader(f"{scheme} Tree")
            pdf.print_page([st_tab[scheme]], big_small=True)
            pdf.write_subheader(f"{scheme} Loci Coverage")

    pdf.output(f'Parapipe_report.pdf', 'F')


if __name__=="__main__":
    main(sys.argv)
