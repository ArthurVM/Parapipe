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


def make_map_table(json_data):
    plt.style.use('ggplot')

    array_names = ["coverage_breadth_hist", "gg_array"]

    table_list = [["ID", "Mean DOC", "Median DOC", "BOC>=5x", "GG area", "Norm GG area", "Av. Qual", "Av. Insert Size", "Fws", "Het. k", "Bound. Prob", "SNPs: total (unique)"]]

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
            
        fws = report["heterozygosity"]["fws"]

        table_list.append([ sample,\
            np.round(mapping_stats["mean_depth_of_coverage"], 1), \
            np.round(mapping_stats["median_depth_of_coverage"], 1), \
            mapping_stats["coverage_breadth_hist"]["5"], \
            np.round(mapping_stats["GG_area"], 3), \
            np.round(mapping_stats["nGG_area"], 3), \
            np.round(mapping_stats["average_quality"], 1), \
            mapping_stats["insert_size_average"], \
            np.round(fws, 3), \
            het_k, \
            het_bp, \
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
        m_area = np.trapz(np.array(coldata), dx=dx)/np.max(w)
        m_areaN = np.trapz(norm_array, dx=dx)/np.max(w)
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