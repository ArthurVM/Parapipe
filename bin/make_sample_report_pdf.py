import os
import io
import sys
import json
import shutil
import numpy as np
import pandas as pd
import calendar
import argparse
from datetime import datetime
from io import StringIO
from fpdf import FPDF
from Bio import Phylo
from matplotlib import pyplot as plt
from matplotlib import rcParams
from pdf import PDF

rcParams['axes.spines.top'] = False
rcParams['axes.spines.right'] = False

def parse_json(json_file):
    print(json_file)
    with open(json_file, 'r') as j:
        report_json = json.loads(j.read())
    return report_json

def pair_jpegs(jpeg_dir):
    s_jpegs = [[os.path.join(jpeg_dir, f)] for f in sorted(os.listdir(jpeg_dir))]
    return s_jpegs

def make_plots(json_data):
    def make_plot_dir():
        PLOT_DIR = 'plots'
        try:
            shutil.rmtree(PLOT_DIR)
            os.mkdir(PLOT_DIR)
        except FileNotFoundError:
            os.mkdir(PLOT_DIR)

    def make_cov_plot(cov_array):
        if cov_array:
            plt.style.use('ggplot')
            data = [b for c, b in cov_array.items()]
            x_pos = [i for i, _ in enumerate(data)]
            x_labels = [c for c, b in cov_array.items()]

            fig = plt.figure(figsize=(8, 3.5), dpi=100)

            plt.bar(x_pos, data, color='blue')
            plt.xticks(x_pos, x_labels)
            plt.xlabel("Depth of coverage")
            plt.ylabel("% Genome covered")
            plt.title("Genome covered to depth")
            plt.tight_layout()

            figout = "./plots/covplot.jpeg"
            plt.savefig(figout)
            return figout
        else:
            return None

    def make_gg_plot(gg_dict):
        if gg_dict:
            gg_df = pd.DataFrame([[int(k), float(v)] for k, v in gg_dict.items()])
            plt.style.use('ggplot')
            GG_f, GG_p = plt.subplots(1, figsize=(8, 4), dpi=100)

            window_array = gg_df.loc[:, 0].tolist()
            w = window_array[1] - window_array[0]
            coldata = gg_df.loc[:, 1].tolist()
            norm_array = [x + 1 - np.max(coldata) for x in coldata]
            m_area = np.trapz(np.array(coldata), dx=w) / np.max(window_array)
            m_areaN = np.trapz(norm_array, dx=w) / np.max(window_array)

            GG_p.plot(window_array, coldata, linestyle="-", linewidth=1, label=f"GG area={np.round(m_area, 3)}\nNormalised area={np.round(m_areaN, 3)}")

            GG_p.set_xlim([0, np.max(window_array) + 1])

            GG_p.set(ylabel="$G$", xlabel="$W$")

            plt.title("Gini Granularity Curve")
            GG_p.legend(ncol=3)
            plt.tight_layout()

            figout = "./plots/gg_curve.jpeg"
            plt.savefig(figout)
            return figout
        else:
            return None

    def make_phylo_plot(treedata, format, prefix):

        plt.style.use('ggplot')
        
        tree = Phylo.read(io.StringIO(nwk), "newick")
        tree.ladderize()

        ## Set figure size for A4 dimensions (210mm x 297mm)
        fig, ax = plt.subplots(figsize=(8.27, 11.69))  # A4 size in inches (1 inch = 25.4mm)

        Phylo.draw(tree, do_show=False, axes=ax)

        ax.axis('off')  # Turn off axes

        plt.savefig(f"{prefix}_tree.png", format='png', dpi=100)  # Adjust dpi as needed
        plt.close()

    def make_map_table(mapping_stats, allele_report):
        if mapping_stats and allele_report:
            plt.style.use('ggplot')
            array_names = ["coverage_breadth_hist", "gg_array"]
            ordered_list = []

            ordered_list.append(["Genome Size", mapping_stats["genome_size"]])
            ordered_list.append(["Mean Depth of Coverage", np.round(mapping_stats["mean_depth_of_coverage"], 3)])
            ordered_list.append(["Median Depth of Coverage", np.round(mapping_stats["median_depth_of_coverage"], 3)])
            ordered_list.append(["GG area", np.round(mapping_stats["GG_area"], 3)])
            ordered_list.append(["Normalised GG area", np.round(mapping_stats["nGG_area"], 3)])
            ordered_list.append(["SNPs", allele_report["total_snps"]])
            ordered_list.append(["unique SNPs", len(allele_report["unique_snps"])])

            ordered_list.append(["Bases mapped (cigar)", int(mapping_stats["bases_mapped_cigar"])])
            ordered_list.append(["Average Quality", int(mapping_stats["average_quality"])])
            ordered_list.append(["Bases Trimmed", int(mapping_stats["bases_trimmed"])])
            ordered_list.append(["Error Rate", np.round(mapping_stats["error_rate"], 3)])
            ordered_list.append(["Average Insert Size", np.round(mapping_stats["insert_size_average"], 1)])
            ordered_list.append(["Insert Size StdDev", np.round(mapping_stats["insert_size_standard_deviation"], 1)])
            ordered_list.append(["Inward Oriented Pairs", int(mapping_stats["inward_oriented_pairs"])])
            ordered_list.append(["Outward Oriented Pairs", int(mapping_stats["outward_oriented_pairs"])])

            fig, ax = plt.subplots(figsize=(8, 4.5), dpi=100)
            fig.patch.set_visible(False)
            ax.axis('off')
            ax.axis('tight')

            ax.table(ordered_list, colLabels=None, loc='center')
            fig.tight_layout()

            figout = "./plots/mapstats.jpeg"
            plt.savefig(figout)
            return figout
        else:
            return None

    def make_st_table(st_stats):
        if st_stats:
            figs = {}

            for scheme, loci in st_stats.items():
                locus_cov = loci["locus_coverage"]
                ordered_list = [["Locus", "DOC", "BOC"]]
                for locus, cov in locus_cov.items():
                    ordered_list.append([locus, np.round(float(cov["DOC"]), 3), np.round(float(cov["BOC"]), 3)])

                fig, ax = plt.subplots(figsize=(8, 4.5), dpi=100)
                fig.patch.set_visible(False)
                ax.axis('off')
                ax.axis('tight')

                ax.table(ordered_list, colLabels=None, loc='center')
                fig.tight_layout()

                figout = f"./plots/{scheme}_stats.jpeg"
                plt.savefig(figout)
                figs[scheme] = figout
            return figs
        else:
            return {}

    def gen_empty_plot(text):
        save_path = figout = f"./plots/SNP_phylo.jpeg"
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, text, ha='center', va='center', fontsize=16, color='red')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel('')
        ax.set_ylabel('')
        plt.savefig(save_path, bbox_inches='tight', pad_inches=0)
        return save_path

    make_plot_dir()
    cov_plot = make_cov_plot(json_data["mapping"]["mapping_stats"]["coverage_breadth_hist"])
    map_tab = make_map_table(json_data["mapping"]["mapping_stats"], json_data["mapping"]["allele_report"])

    ## check typing schemes were defined
    if json_data["typing"] != None:
        st_tab = make_st_table(json_data["typing"])

        st_plots = []
        for st, data in json_data.get("sequence_phylo", {}).items():
            if data and data.get("contree") != "NA":
                st_plot = make_phylo_plot(data["contree"], "newick", st)
                st_plots.append([st, st_plot])

    ## else escape in silico molecular typing reporting
    else:
        st_tab = None
        st_plots = None

    gg_curve = make_gg_plot(json_data["mapping"]["mapping_stats"].get("gg_array", {}))
    
    if json_data["mapping"]["mapping_stats"]["coverage_breadth_hist"]["5"] < 95.0:
        snp_phylo_plot = gen_empty_plot("Coverage too low for inclusion in SNP tree")
    else:
        snp_phylo_plot = "./allele_matrix_tree.png"

    return map_tab, cov_plot, gg_curve, snp_phylo_plot, st_plots, st_tab

def make_moi_table(moi_stats):
    if moi_stats:
        plt.style.use('ggplot')
        ordered_list = []

        if moi_stats.get("fit_score") is not None or moi_stats.get("boundary_prob") is not None:
            if moi_stats["fit_score"] >= 0.9 and moi_stats["boundary_prob"] >= 0.9:
                het_bool = "Yes"
            else:
                het_bool = "No"
        else:
            het_bool = "No"

        ordered_list.append(["k", str(moi_stats.get("k", ""))])
        ordered_list.append(["Fit Score", str(moi_stats.get("fit_score", ""))])
        ordered_list.append(["Boundary Probability", str(moi_stats.get("boundary_prob", ""))])
        ordered_list.append(["Suspected Heterozygosity", het_bool])

        fig, ax = plt.subplots(figsize=(8, 4.5), dpi=100)
        fig.patch.set_visible(False)
        ax.axis('off')
        ax.axis('tight')
        ax.table(ordered_list, colLabels=None, loc='center')
        fig.tight_layout()

        figout = "./plots/moistats.jpeg"
        plt.savefig(figout)
        return figout
    else:
        return None

def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('script_path', action='store', help=argparse.SUPPRESS)
    parser.add_argument('--env', action='store', nargs=1, help='JSON containing run environment info.')
    parser.add_argument('--report', action='store', nargs=1, help='JSON report produced by Parapipe.')
    parser.add_argument('--id', action='store', nargs=1, help='Sample ID.')
    parser.add_argument('--moi_json', action='store', nargs=1, help='MOI JSON report.')
    parser.add_argument('--moi_plots', action='store', nargs=2, help='MOI plots.')
    parser.add_argument('--png', action='store', nargs="*", help='SNP PNG trees.')
    args = parser.parse_args(argv)
    return args

def main(argv):
    args = parse_args(argv)
    env = parse_json(args.env[0])
    report = parse_json(args.report[0])
    moi_stats = parse_json(args.moi_json[0])
    
    map_tab, cov_plot, gg_curve, snp_phylo_plot, st_plots, st_tab = make_plots(report)
    moi_tab = make_moi_table(moi_stats)
    
    pdf = PDF()
    pdf.print_page([map_tab, cov_plot, gg_curve])
    pdf.write_title(args.id[0])
    pdf.write_subheader("wgSNP tree")
    pdf.print_page(["./allele_matrix_pca.png"])
    pdf.write_subheader("wgSNP PCA")
    pdf.print_page([snp_phylo_plot])
    pdf.write_subheader("wgSNP Tree")
    
    ## check typing scheme analysis was carried out and print plots
    if st_plots != None and st_tab != None:
        
        for scheme, st_plot in st_plots:
            pdf.print_page([st_plot, st_tab[scheme]], big_small=True)
            pdf.write_subheader(f"{scheme}")
    
    pdf.print_page([moi_tab, moi_stats.get("cluster_fig_path", "")], small_big=True)
    pdf.write_subheader("Heterozygosity")
    pdf.print_page([moi_stats.get("fit_fig_path", "")])
    pdf.output(f'{args.id[0]}_report.pdf', 'F')

if __name__ == "__main__":
    main(sys.argv)
