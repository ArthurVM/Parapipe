import os
import json
import pandas as pd
import numpy as np
import argparse
from io import StringIO
import plotly.express as px
import plotly.graph_objects as go

from plotting_funcs import *
from phylo_plotting_funcs import *

fig_height = 750


def write_html(reports, map_tab, gg_fig, qc_report, env, dist_matrix, run_samples):

    dist_prefix = path.basename(dist_matrix.split(".csv")[0])
    dist_df = pd.read_csv(dist_matrix, index_col=0)

    ## make logo embedding
    scriptdir = os.path.dirname(os.path.realpath(__file__))
    logo_path = f"{scriptdir}/../resources/./PP_logo.png"
    logo_prefix, logo_b64 = embed_image(logo_path)

    run_dt = env["date"] + " " + env["time"]
    report_dt = get_dt()

    ## get phylo plot embeddings
    phylo_plots_wrapper = make_WG_plots(dist_df, dist_prefix, dist_df.index)

    ## include parameters
    params_line = ""
    for key, value in env["params"].items():
        params_line += f"<li><b>{key}</b>: {value}</li>\n"

    ## make map table
    map_table = go.Table(
        header=dict(values=map_tab.columns),
        cells=dict(values=map_tab.values.T.tolist()))
    map_tab_fig = go.Figure(data=[map_table])
    map_tab_fig.update_layout(height=fig_height)

    ## make MOI plots
    moi_htmls = ""
    moi_html_links = []
    for prefix, report in reports.items():
        h_rep = report["heterozygosity"]
        moi_htmls += make_html_moi_plot(prefix, os.path.join(h_rep["het_dir"], h_rep["cluster_fig_path"]))
        moi_html_links.append(content_link(prefix, prefix, chunk_class="content-panelh2"))
    
    ## Generate content links
    content_links = "\n".join([
        content_link("Mapping Statistics", "map"),
        content_link("Gini-Granularity Curves", "gini"),
        content_link("Phylogenetic Analysis", "phylo"),
        content_link("SNP Report", "snp"),
        content_link("MOI Analysis", "moi"),
        "\n".join(moi_html_links)
    ])

    ## make stacked bar plot for SNPs
    snp_bar_dict = pd.DataFrame({ prefix : { "total_snps" : d["mapping"]["allele_report"]["total_snps"], "non_unique_snps" : d["mapping"]["allele_report"]["total_snps"] - len(d["mapping"]["allele_report"]["unique_snps"]), "unique_snps" : len(d["mapping"]["allele_report"]["unique_snps"]) } for prefix, d in reports.items() }).T
    stacked_bar = make_stacked_bar(snp_bar_dict, fig_height)

    ## make html
    report_html = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Parapipe Run Report</title>
        <style>
        body {{
            font-family: Helvetica, sans-serif;
        }}
        .content-panel_wrapper {{
            flex: 0 0 220px;
            padding: 10px;
            border-right: 3px solid darkgray;
            position: fixed;
            top: 0;
            left: 0;
            height: 100vh;
            background-color: #ededed;
            overflow-x: auto;
        }}
        .content-panel h1 {{
            font-size: 18px;
            text-align: center;
            margin: 0;
            border-bottom: 1px solid #ccc;
        }}
        .content-panel h1 a {{
            padding: 20px 0 15px;
        }}
        .content-panel h1 img {{
            height: 26px;
        }}
        .content-panel h1 small {{
            font-size: 12px;
        }}
        .content-panel .side-nav-title a {{
            color: #333;
            font-size: 16px;
            font-weight: normal;
            padding: 15px 0;
        }}
        .content-panel .mqc_loading_warning {{
            text-align: center;
            border-bottom: 1px solid #ccc;
            color: #ca424c;
        }}
        .content-panel p {{
            font-size: 0.8em;
            color: #999;
            padding: 10px;
        }}
        .content-panel a:hover,
        .content-panel a:active,
        .content-panel a:focus {{
            background-color: #dedede;
        }}
        .content-link {{
            display: block;
            margin-bottom: 5px;
        }}
        .content-panelh2 {{
            font-size: 12px;
            text-align: left;
            margin: 10px;
        }}
        .main-content {{
            flex: 1 0 0;
            padding: 10px;
            margin-left: 230px;
        }}
        .plot-container {{
            display: flex;
            flex-direction: row;
            margin-bottom: 10px;
        }}
        .full-width-image {{
            width: 100vw;
            height: 1000px;
        }}
        .img-column3x {{
            flex: 1 0 33.33%;
            margin-bottom: 10px;
        }}
        .plots-container {{
            display: flex;
            justify-content: space-between;
            width: 100%;
        }}
        .plot-wrapper {{
            flex: 0 0 45%;
            margin: 10px;
        }}
        .plot-prefix {{
            font-weight: bold;
            margin-right: 10px;
            cursor: pointer;
        }}
        .plot-image {{
            width: 400px;
        }}
        .disable-scrolling{{
            height: 100%;
            overflow: hidden;
        }}
        .hidden {{
            display: none;
        }}
        </style>
    </head>
    
    <div class="content-panel_wrapper">
        <div class="content-panel">
            <a href="https://github.com/ArthurVM/Parapipe"><img src="data:image/{prefix};base64,{logo_b64}" width=200></a>
            <hr style="width:80%;text-align:left;"> 
            <h2>Contents</h2>
            {content_links}
        </div>
    </div>

    <div class="main-content">
        <a href="https://github.com/ArthurVM/Parapipe"><img src="data:image/{prefix};base64,{logo_b64}" width=300></a>
        <p>Run initiated at {run_dt} BST</p>
        <p>Report generated at {report_dt} BST</p>

        <hr style="width:50%;text-align:left;margin-left:0"> 
        <h3>Input Parameters</h3>
        {params_line}
        <hr style="width:50%;text-align:left;margin-left:0"> 
        
        <div id="map">
        <h2>Mapping Statistics</h2>

        <p>Mapping statistics for samples within this Parapipe run.</p>
        <p>MultiQC report: <a href={qc_report[0]}>{qc_report[0]}</a></p>

        {map_tab_fig.to_html(full_html=False)}

        
        <div id="gini">
        <h2>Gini-Granularity Curves</h2>
        <p>GG curves characterise coverage inequality. Low GG-area indicates low coverage equality. Low normalised GG-area indicates reads aggregate in small\
            regions, whereas high normalised GG-area indicates reads aggregate in larger regions <a href="https://link.springer.com/chapter/10.1007/978-3-030-46970-2_11">[1]</a>.</p>
        {gg_fig.to_html(full_html=False)}
        

        <div id="phylo">
        <h2>Phylogenetic Analysis</h2>
        <p>All analysis carried out using absolute <a href="https://en.wikipedia.org/wiki/Hamming_distance">(Hamming)</a> SNP distance across the whole genome.</p>
        <p>Key:</p>
        <p>red = samples processed this run.</p>
        <p>grey = samples included in the input archive.</p>
        
        {phylo_plots_wrapper}

        <div id="snp">
        <h2>SNP Report</h2>
        {stacked_bar.to_html(full_html=False)}

        
        <div id="moi">
        <h1>MOI Analysis</h1>
        {moi_htmls}

    </div>

    </body>

    </html>
    """

    # Write HTML report to file
    with open('Parapipe_report.html', 'w') as f:
        f.write(report_html)


def parse_args(argv):

    parser = argparse.ArgumentParser()

    parser.add_argument('script_path', action='store', help=argparse.SUPPRESS)
    parser.add_argument('--env', action='store', nargs=1, help='JSON containing run environment info.')
    parser.add_argument('--report', action='store', nargs="*", help='JSON reports produced by Parapipe.')
    parser.add_argument('--moi_json', action='store', nargs="*", help='MOI JSON reports.')
    parser.add_argument('--multiqc', action='store', nargs="*", help='MultiQC report path.')
    parser.add_argument('--dist_matrix', action='store', help='Distance matrix constructed by Ardal.')

    args = parser.parse_args(argv)

    return args


def main(argv):

    args = parse_args(argv)

    env = parse_json(args.env[0])

    reports = {}
    run_samples = []
    for rep in args.report:
        prefix = os.path.basename(rep).split("_report.json")[0]
        reports[prefix] = parse_json(rep)
        run_samples.append(prefix)

    for rep in args.moi_json:
        prefix = os.path.basename(rep).split(".moi.json")[0]
        reports[prefix]["heterozygosity"] = parse_json(rep)
        reports[prefix]["heterozygosity"]["het_dir"] = "/".join(os.path.abspath(rep).split("/")[:-1])

    map_tab = make_map_table(reports)
    map_df = pd.DataFrame(map_tab[1:], columns=map_tab[0])

    map_df.to_csv("./run_results.csv")

    gg_dict = { k : {int(w) : g for w, g in v["mapping"]["mapping_stats"]["gg_array"].items()} for k, v in reports.items() }
    gg_fig = make_gg_plot(gg_dict, fig_height)

    ## construct html
    write_html(reports, map_df, gg_fig, args.multiqc, env, args.dist_matrix, run_samples)

main(sys.argv)