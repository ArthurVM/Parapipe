import os
import sys
import json
from collections import defaultdict
from datetime import datetime
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import plotly.express as px
from sklearn.decomposition import PCA
import base64
import re


def content_link(section_name, section_id, chunk_class="nav-l1"):
        return f"""<li>
            <a href="#{section_id}" class="{chunk_class}">{section_name}</a>
        </li>"""


def embed_image(image_path):
  """
  Embeds an image into an HTML file using Base64 encoding.

  Args:
    html_file: Path to the HTML file.
    image_path: Path to the image file.
  """
  with open(image_path, "rb") as image_file:
    image_data = image_file.read()
  encoded_data = base64.b64encode(image_data).decode()

  return [image_path.split(".")[-1], encoded_data]

  # Read and modify HTML content
  with open(html_file, "r+") as f:
    html_content = f.read()
    new_content = html_content.replace(
        'src="' + image_path + '"',  # Replace the original src attribute
        f'src="data:image/{image_path.split(".")[-1]};base64,{encoded_data}"'
    )
    f.seek(0)
    f.write(new_content)
    f.truncate()


def chunks(lst, n):
  """Yields successive n-sized chunks from lst."""
  for i in range(0, len(lst), n):
    yield lst[i:i + n]
    

def get_dt():
    now = datetime.now()

    formatted_datetime = now.strftime("%Y-%m-%d %H:%M:%S")
    return formatted_datetime


def parse_json(json_file):
    print(json_file)
    with open(json_file, 'r') as j:
        report_json = json.loads(j.read())
    return report_json


def create_distance_heatmap(distance_data, row_labels=None, col_labels=None, fig_height=1000):
  """
  Generates a Plotly heatmap from a distance matrix (list or pandas DataFrame).

  Args:
      distance_data (list or pd.DataFrame): The distance matrix data.
      row_labels (list, optional): List of labels for rows (default: None).
      col_labels (list, optional): List of labels for columns (default: None).

  Returns:
      plotly.graph_objects.Figure: The interactive heatmap figure.
  """

  ## Handle pandas DataFrame input
  if isinstance(distance_data, pd.DataFrame):
    distance_matrix = distance_data.to_numpy()
    if row_labels is None:
      row_labels = list(distance_data.columns)
    if col_labels is None:
      col_labels = row_labels.copy()

  ## Handle list-based input (assuming square matrix)
  elif isinstance(distance_data, np.ndarray):
    if not all(len(row) == len(distance_data) for row in distance_data):
      raise ValueError("Distance matrix must be a square matrix (numpy array).")
    distance_matrix = distance_data
    if row_labels is None:
      row_labels = [f"Sample {i+1}" for i in range(len(distance_data))]
    if col_labels is None:
      col_labels = row_labels.copy()
  else:
    raise TypeError("Distance data must be a numpy array or pandas DataFrame.")

  ## Create trace for the heatmap
  trace = go.Heatmap(
      z=distance_matrix,
      colorbar=dict(title="Distance"),
      showscale=True
  )

  ## Add row and column labels (if provided)
  if row_labels:
    trace.update(x=row_labels)
  if col_labels:
    trace.update(y=col_labels[::-1])  # Reverse column labels for bottom-up display

  ## Configure layout
  layout = go.Layout(
      title="wgSNP Distance Heatmap",
  )

  ## Create and return the figure
  fig = go.Figure(data=[trace], layout=layout)
  fig.update_layout(height=fig_height, width=fig_height)

  return fig


def create_PCA(data_df, run_samples, fig_height=750):
   
    pca = PCA(n_components=2)  # Extract 2 principal components
    components = pca.fit_transform(data_df)

    ## Create scatter plot for PCA components
    trace = go.Scatter(
        x=components[:, 0],
        y=components[:, 1],
        mode='markers',
        text=data_df.index,  ## Display sample labels on hover
        marker=dict(size=[15 if i in run_samples else 7.5 for i in data_df.index], \
                    color=['red' if i in run_samples else 'grey' for i in data_df.index], \
                    opacity=[0.9 if i in run_samples else 0.4 for i in data_df.index])
    )

    ## Configure layout
    layout = go.Layout(
        title="wgSNP PCA Plot",
        xaxis_title="PC1",
        yaxis_title="PC2",
        legend_title_text="Samples"
    )

    fig = go.Figure(data=[trace], layout=layout)
    fig.update_layout(height=fig_height, width=fig_height)

    return fig


def create_PCA_3D(data_df, run_samples, fig_height=750):
  
  pca = PCA(n_components=3)  # Extract 3 principal components
  components = pca.fit_transform(data_df)

  ## Create scatter plot for PCA components (3D)
  trace = go.Scatter3d(
        x=components[:, 0],
        y=components[:, 1],
        z=components[:, 2],
        mode='markers',
        text=data_df.index,
        marker=dict(size=[15 if i in run_samples else 7.5 for i in data_df.index], \
                    color=['red' if i in run_samples else 'grey' for i in data_df.index],
                    line=dict(width=0, color='white'),
                    symbol=['circle' if i in run_samples else 'square-open' for i in data_df.index])
        )

  ## Configure layout for 3D scene
  layout = go.Layout(
      title="wgSNP PCA Plot (3PCs)",
      scene = dict(
          xaxis_title="PC1",
          yaxis_title="PC2",
          zaxis_title="PC3"
      ),
      legend_title_text="Samples"
  )

  fig = go.Figure(data=[trace], layout=layout)
  fig.update_layout(height=fig_height)  # Adjust for 3D scene (width not required)

  return fig


def make_html_moi_plot(prefix, path):
    return f"""
    <div class="img-column">
      <p class="img-prefix">{prefix}</p>
      {path}
    </div>
    """


def make_html_moi_plot(prefix, path):
    # return f"""
    # <div class="plot-container">\
    #   <p class="plot-prefix">{prefix}</p>\
    #   <img src="{path}" alt="{prefix} class="plot-image">\
    # </div>
    # """
    plot_prefix, plot_b64 = embed_image(path)
    return f"""
    <div id="{prefix}">
    <p class="plot-prefix">{prefix}</p>\
    <img src="data:image/{plot_prefix};base64,{plot_b64}">
    """ 


def make_stacked_bar(df, fig_height):

    colors = ['royalblue', 'skyblue']  # Adjust colors as desired
    # Define a list to store traces (plots)
    traces = []

    # Create traces (one for each category) and stack them
    for i, col in enumerate(df[["non_unique_snps", "unique_snps"]].columns):  # Skip the first column (total_snps)
        traces.append(
            go.Bar(
                x=df.index,  # Sample IDs as x-axis labels
                y=df[col],  # Data for the current category (non_unique_snps or unique_snps)
                name=col,  # Category name as legend entry
                marker=dict(color=colors[i])  # Assign color from the color list
            )
        )

    # Configure plot layout
    layout = go.Layout(
        title="Stacked Bar Plot of SNPs per Sample",
        xaxis_title="Sample ID",
        yaxis_title="Number of SNPs",
        barmode='stack'  # Set bars to stack on top of each other
    )

    # Create the figure with traces and layout
    fig = go.Figure(data=traces, layout=layout)
    fig.update_layout(height=fig_height)

    return fig


def _natural_sort_key(value: str):
    """Helper to sort chromosome labels naturally (e.g. 1,2,10)."""
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", str(value))]


def load_chrom_lengths_from_fasta(fasta_path: str) -> dict:
    """
    Read a FASTA and return a dict of chromosome -> length.
    """
    lengths = {}
    chrom = None
    length = 0
    with open(fasta_path, "r") as fin:
        for line in fin:
            if line.startswith(">"):
                if chrom is not None:
                    lengths[chrom] = length
                chrom = line[1:].split()[0]
                length = 0
            else:
                length += len(line.strip())
    if chrom is not None:
        lengths[chrom] = length
    return lengths


def plot_maf_across_genome(
    vcf_obj,
    sample_id: str,
    out_path,
    min_depth: int = 5,
    chrom_gap: int = 0,
    max_points: int = 100_000,
    chrom_lengths: dict | None = None,
):
    """
    Scatter plot of minor allele frequency across the genome for a single sample.
    - vcf_obj: VCF instance providing allele depths.
    - sample_id: sample to plot.
    - out_path: file path to save the PNG.
    - min_depth: minimum total depth per site to include.
    - chrom_gap: spacing inserted between concatenated chromosomes (bp); 0 for no gap.
    - max_points: cap plotted points by downsampling.
    - chrom_lengths: optional mapping of chromosome -> length to enforce consistent spacing.
    """
    counts = vcf_obj.getAlleleCounts(sample_id)
    if counts.empty:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "No allele counts available", ha="center", va="center")
        ax.axis("off")
        fig.savefig(out_path, bbox_inches="tight")
        plt.close(fig)
        return out_path

    counts["depth"] = counts["ref_depth"] + counts["alt_depth"]
    counts = counts[counts["depth"] >= min_depth]
    if counts.empty:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, f"No sites with depth >= {min_depth}", ha="center", va="center")
        ax.axis("off")
        fig.savefig(out_path, bbox_inches="tight")
        plt.close(fig)
        return out_path

    counts["maf"] = counts[["ref_freq", "alt_freq"]].min(axis=1)
    counts = counts.sort_values(["chromosome", "position"])

    if len(counts) > max_points:
        stride = int(np.ceil(len(counts) / max_points))
        counts = counts.iloc[::stride]

    if chrom_lengths:
        chroms = sorted(chrom_lengths.keys(), key=_natural_sort_key)
        missing = [c for c in counts["chromosome"].unique() if c not in chrom_lengths]
        if missing:
            chroms.extend(sorted(missing, key=_natural_sort_key))
    else:
        chroms = sorted(counts["chromosome"].unique(), key=_natural_sort_key)

    color_cycle = ["#000000", "#7f7f7f"]

    xs = []
    ys = []
    tick_positions = []
    tick_labels = []
    offset = 0
    genome_length = 0

    for chrom_idx, chrom in enumerate(chroms):
        chrom_df = counts[counts["chromosome"] == chrom]
        chrom_length = None
        if chrom_lengths:
            chrom_length = chrom_lengths.get(chrom)

        if chrom_df.empty:
            length_for_spacing = chrom_length if chrom_length is not None else 0
            genome_length = offset + length_for_spacing
            tick_positions.append(offset + length_for_spacing / 2 if length_for_spacing > 0 else offset)
            tick_labels.append(chrom)
            offset += length_for_spacing + chrom_gap
            continue

        chrom_positions = chrom_df["position"].to_numpy()
        if chrom_length is None:
            chrom_length = chrom_positions.max()
        else:
            chrom_length = max(chrom_length, chrom_positions.max())

        chrom_positions = chrom_positions + offset
        xs.append(chrom_positions)
        ys.append(chrom_df["maf"].to_numpy())
        tick_positions.append(offset + chrom_length / 2)
        tick_labels.append(chrom)

        offset += chrom_length + chrom_gap
        genome_length = offset

    fig, ax = plt.subplots(figsize=(10, 3))
    for idx, (x_vals, y_vals) in enumerate(zip(xs, ys)):
        ax.scatter(x_vals, y_vals, s=6, color=color_cycle[idx % len(color_cycle)], alpha=0.8, linewidths=0)

    ax.set_ylim(-0.02, 0.52)
    ax.set_ylabel("MAF")
    ax.set_xlabel("Position")
    if tick_positions:
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels, rotation=90, fontsize=8)
    ax.set_title(f"{sample_id}")
    if genome_length > 0:
        ax.set_xlim(0, genome_length)
    ax.margins(x=0)
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    return out_path


def make_map_table(json_data):
    plt.style.use('ggplot')

    table_list = [["ID", "Mean DOC", "Median DOC", "BOC>=5x", "GG area", "Norm GG area", "Av. Qual", "Av. Insert Size", "Fws", "Hs", "MAF plot", "SNPs: total (unique)"]]

    for sample, report in json_data.items():
        mapping_stats = report["mapping"]["mapping_stats"]
        al_stats = report["mapping"]["allele_report"]

        het = report.get("heterozygosity", {}) or {}
        fws = het.get("fws")
        hs = het.get("hs")
        maf_plot = het.get("maf_plot") or het.get("maf_plot_path") or "NA"

        table_list.append([ sample,\
            np.round(mapping_stats["mean_depth_of_coverage"], 1), \
            np.round(mapping_stats["median_depth_of_coverage"], 1), \
            mapping_stats["coverage_breadth_hist"]["5"], \
            np.round(mapping_stats["GG_area"], 3), \
            np.round(mapping_stats["nGG_area"], 3), \
            np.round(mapping_stats["average_quality"], 1), \
            mapping_stats["insert_size_average"], \
            "NA" if fws is None else np.round(fws, 3), \
            "NA" if hs is None else np.round(hs, 3), \
            maf_plot, \
            f"{al_stats['total_snps']} ({len(al_stats['unique_snps'])})"])
    
    return table_list


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


def make_gg_plot(gg_dict, fig_height):
    gg_df = pd.DataFrame(gg_dict)

    ## Create a list of column names (excluding index)
    column_names = gg_df.columns
    w = list(int(i) for i in gg_df.index)
    traces = []

    ## Define hover template with column names and values
    for col in column_names:
        coldata = gg_df[col]
        dx=50
        norm_array = [x+1-np.max(coldata) for x in coldata]
        m_area = np.trapezoid(np.array(coldata), dx=dx)/np.max(w)
        m_areaN = np.trapezoid(norm_array, dx=dx)/np.max(w)
        hover_template = f"GG area={np.round(m_area, 3)}\nNormalised GG area={np.round(m_areaN, 3)}"
        traces.append(
            dict(
                x=gg_df.index,
                y=gg_df[col],
                name=col,
                hovertemplate=hover_template,
                mode="lines",
            )
        )

    # Create the figure
    fig = go.Figure(data=traces)
    fig.update_layout(height=fig_height)

    # Display the plot
    return fig
