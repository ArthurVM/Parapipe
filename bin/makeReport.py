import os
import sys
import json
import shutil
import numpy as np
import pandas as pd
import calendar
from datetime import datetime
from io import StringIO
from fpdf import FPDF
from Bio import Phylo
from matplotlib import pyplot as plt
from matplotlib import rcParams

rcParams['axes.spines.top'] = False
rcParams['axes.spines.right'] = False

class PDF(FPDF):
    def __init__(self):
        super().__init__()
        self.WIDTH = 210
        self.HEIGHT = 297

    def write_subheader(self, text):
        """ write subheader text to pdf
        """
        self.set_font('Courier', 'B', 15)
        self.cell(self.WIDTH - 80)
        self.cell(-1, -1, text, 0, 0, 'R')
        self.ln(20)

    def write_title(self, sample_ID):

        self.set_font('Courier', 'B', 18)
        self.cell(self.WIDTH - 80)
        self.cell(-20, -10, f"Report for {sample_ID}", 0, 0, 'R')

    def header(self):
        # Custom logo and positioning
        # Create an `assets` folder and put any wide and short image inside
        # Name the image `logo.png`
        scriptdir = os.path.dirname(os.path.realpath(__file__))
        self.image(f'{scriptdir}/../resources/phw_logo.jpeg', 10, 8, 60)
        self.set_font('Courier', 'B', 11)
        self.cell(self.WIDTH - 80)
        self.cell(60, 1, 'Parapipe report', 0, 0, 'R')
        self.ln(20)

    def footer(self):
        # Page numbers in the footer
        self.set_y(-15)
        self.set_font('Courier', 'I', 8)
        self.set_text_color(128)
        self.cell(0, 10, 'Page ' + str(self.page_no()), 0, 0, 'C')

    def page_body(self, images):
        # Determine how many plots there are per page and set positions
        # and margins accordingly
        if len(images) == 3:
            self.image(images[0], 15, 25, self.WIDTH - 30)
            self.image(images[1], 15, self.WIDTH / 2 + 5, self.WIDTH - 30)
            self.image(images[2], 15, self.WIDTH / 2 + 90, self.WIDTH - 30)
        elif len(images) == 2:
            self.image(images[0], 15, 25, self.WIDTH - 30)
            self.image(images[1], 15, self.WIDTH / 2 + 50, self.WIDTH - 30)
        else:
            self.image(images[0], 15, 40, self.WIDTH - 30)

    def print_page(self, images):
        # Generates the report
        self.add_page()
        self.page_body(images)

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

def make_plots(json_data):

    def make_plot_dir():
        PLOT_DIR = 'plots'

        # Delete folder if exists and create it again
        try:
            shutil.rmtree(PLOT_DIR)
            os.mkdir(PLOT_DIR)
        except FileNotFoundError:
            os.mkdir(PLOT_DIR)

    def make_cov_plot(cov_array):
        """ make the coverage plot
        """
        plt.style.use('ggplot')

        data = [b for c, b in cov_array.items()]
        x_pos = [i for i, _ in enumerate(data)]
        x_labels = [c for c, b in cov_array.items()]

        fig = plt.figure(figsize=(8, 3.5), dpi=100)

        plt.bar(x_pos, data, color='blue' )
        plt.xticks(x_pos, x_labels)
        plt.xlabel("Depth of coverage")
        plt.ylabel("% Genome covered")
        plt.title("Genome covered to depth")

        plt.tight_layout()

        figout = "./plots/covplot.jpeg"
        plt.savefig(figout)

        return figout

    def make_gg_plot(gg_dict):
        gg_df = pd.DataFrame([ [int(k), float(v)] for k, v in json_data["mapping_stats"]["gg_array"].items() ])

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

    def make_phylo_plot(treedata):

        plt.style.use('ggplot')

        handle = StringIO(treedata)
        tree = Phylo.read(handle, json_data["phylo"]["format"])
        tree.ladderize()

        fig = plt.figure(figsize=(15, 22), dpi=100)
        axes = fig.add_subplot(1, 1, 1)
        Phylo.draw(tree, axes=axes, do_show=False)

        plt.title("SNP Distance")

        figout = "./plots/phylo.jpeg"
        plt.savefig(figout)

        return figout

    def make_map_table(mapping_stats):
        plt.style.use('ggplot')

        array_names = ["coverage_breadth", "gg_array"]

        ordered_list = []

        ## manually insert fields into the list to maintain some meaningful order
        ## mapping fields
        ordered_list.append(["Genome Size", mapping_stats["genome_size"]])
        ordered_list.append(["Mean Depth of Coverage", np.round(mapping_stats["mean_depth_of_coverage"],3)])
        ordered_list.append(["Median Depth of Coverage", np.round(mapping_stats["median_depth_of_coverage"],3)])
        ordered_list.append(["GG area", np.round(mapping_stats["GG_area"],3)])
        ordered_list.append(["Normalised GG area", np.round(mapping_stats["nGG_area"],3)])

        ## quality fields
        ordered_list.append(["Bases mapped (cigar)", int(mapping_stats["bases_mapped_cigar"])])
        ordered_list.append(["Average Quality", int(mapping_stats["average_quality"])])
        ordered_list.append(["Bases Trimmed", int(mapping_stats["bases_trimmed"])])
        ordered_list.append(["Error Rate", np.round(mapping_stats["error_rate"],3)])
        ordered_list.append(["Average Insert Size", np.round(mapping_stats["insert_size_average"],1)])
        ordered_list.append(["Insert Size StdDev", np.round(mapping_stats["insert_size_standard_deviation"],1)])
        ordered_list.append(["Inward Oriented Pairs", int(mapping_stats["inward_oriented_pairs"])])
        ordered_list.append(["Outward Oriented Pairs", int(mapping_stats["outward_oriented_pairs"])])

        # df = pd.DataFrame(ordered_list)

        # hide axes
        fig, ax = plt.subplots(figsize=(8, 4.5), dpi=100)

        # hide axes
        fig.patch.set_visible(False)
        ax.axis('off')
        ax.axis('tight')

        ax.table(ordered_list, colLabels=None, loc='center')

        fig.tight_layout()

        figout = "./plots/mapstats.jpeg"
        plt.savefig(figout)

        return figout

    make_plot_dir()
    cov_plot = make_cov_plot(json_data["mapping_stats"]["coverage_breadth"])
    phylo_plot = make_phylo_plot(json_data["phylo"]["tree"])
    map_tab = make_map_table(json_data["mapping_stats"])
    gg_curve = make_gg_plot(json_data["mapping_stats"]["gg_array"])

    return [[map_tab, cov_plot, gg_curve]], [[phylo_plot]]

def main(json_file, jpeg_dir, sample_ID):
    json_data = parse_json(json_file)
    paired_covplots = pair_jpegs(jpeg_dir)

    plots, phylo_plots = make_plots(json_data)

    pdf = PDF()

    for i, elem in enumerate(plots):
        pdf.print_page(elem)
        if i == 0:
            pdf.write_title(os.path.basename(json_file).split("_SNPreport")[0])

    for elem in phylo_plots:
        pdf.print_page(elem)
        pdf.write_subheader("Phylogenetics")

    for i, elem in enumerate(paired_covplots):

        chr = os.path.basename(elem[0]).split("_chrom")[0]
        pdf.print_page(elem)
        pdf.write_subheader(f"Coverage plots for chromosome {chr}")

    pdf.output(f'{sample_ID}_report.pdf', 'F')

if __name__=="__main__":
    # main("/home/amorris/BioInf/Parapipe/C_hom_SMALL_ass/UKH37_4/UKH37_4_SNPreport.json", \
    #     "/home/amorris/BioInf/Parapipe/report_WD/chromplots/")
    main(sys.argv[1], sys.argv[2], sys.argv[3])
