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


def summarise_gp60(reports):
    """
    Build summary and long-form tables indicating gp60 polyclonality per sample.
    """
    summary_rows = []
    stack_rows = []
    for sample, rep in reports.items():
        raw_poly = rep.get("gp60_polyclonality") or {}
        if isinstance(raw_poly, dict) and raw_poly:
            first_val = list(raw_poly.values())[0]
            poly = first_val if isinstance(first_val, dict) else {}
        else:
            poly = {}

        fam_counts = {}
        subtype_counts = []

        for fam, entries in poly.items():
            if not isinstance(entries, list):
                continue
            fam_counts[fam] = len(entries)
            for entry in entries:
                if not entry or len(entry) < 3:
                    continue
                _, subtype, count = entry
                try:
                    cnt_val = float(count)
                except Exception:
                    cnt_val = 0
                subtype_counts.append((fam, subtype, cnt_val))

        total_counts = sum(c for _, _, c in subtype_counts) if subtype_counts else 0
        summary_rows.append(
            {
                "sample_id": sample,
                "families": len(fam_counts),
                "allele_count": sum(fam_counts.values()),
                "polyclonal": sum(fam_counts.values()) > 1,
            }
        )

        for fam, subtype, cnt_val in subtype_counts:
            proportion = (cnt_val / total_counts) if total_counts else 0
            stack_rows.append(
                {
                    "sample_id": sample,
                    "family": fam,
                    "subtype": subtype,
                    "allele_count": cnt_val,
                    "proportion": proportion,
                }
            )

    return pd.DataFrame(summary_rows), pd.DataFrame(stack_rows)


def make_gp60_bar(df: pd.DataFrame):
    if df.empty:
        return None
    df = df.copy()

    families = df["family"].unique().tolist()
    base_colors = px.colors.qualitative.Set2

    subtype_colors = {}
    for row in df.itertuples():
        base = base_colors[families.index(row.family) % len(base_colors)] if families else "#777777"
        subtype_colors[row.subtype] = base

    fig = px.bar(
        df,
        x="sample_id",
        y="proportion",
        color="subtype",
        color_discrete_map=subtype_colors,
        title="gp60 polyclonality (subtype proportions by family)",
        labels={"proportion": "Allele proportion", "sample_id": "Sample ID", "subtype": "Subtype"},
        hover_data={"family": True, "allele_count": True},
    )
    fig.update_traces(marker_line_color="white", marker_line_width=1)
    fig.update_layout(height=fig_height, xaxis_tickangle=45, hovermode="x unified", barmode="stack", legend_title_text="Subtype (color by family)")
    return fig


def write_html(reports, map_tab, gg_fig, qc_report, env, dist_matrix, run_samples, gp60_summary=None, gp60_fig=None):

    dist_prefix = os.path.basename(dist_matrix.split(".csv")[0])
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
        h_rep = report.get("heterozygosity", {}) or {}
        maf_plot_path = h_rep.get("maf_plot_path") or h_rep.get("maf_plot")
        if maf_plot_path:
            moi_htmls += make_html_moi_plot(prefix, maf_plot_path)
        else:
            moi_htmls += f"<p>No MOI plot available for {prefix}</p>"
        moi_html_links.append(content_link(prefix, prefix, chunk_class="content-panelh2"))
    
    ## Generate content links
    content_links = "\n".join(
        [
            content_link("Run Summary", "summary"),
            content_link("Mapping Statistics", "map"),
            content_link("Gini-Granularity Curves", "gini"),
            content_link("Phylogenetic Analysis", "phylo"),
            content_link("SNP Report", "snp"),
            content_link("gp60 Subtyping", "gp60"),
            content_link("MOI Analysis", "moi"),
            "\n".join(moi_html_links),
        ]
    )

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
            font-family: "Helvetica Neue", Arial, sans-serif;
            background-color: #f6f7fb;
            color: #1f2d3d;
        }}
        a {{ color: #1e74c6; text-decoration: none; }}
        a:hover {{ text-decoration: underline; }}
        .content-panel_wrapper {{
            position: fixed;
            top: 0;
            left: 0;
            height: 100vh;
            width: 230px;
            padding: 14px;
            background: #e9ecf3;
            border-right: 1px solid #d1d7e0;
            overflow-y: auto;
        }}
        .content-panel h2 {{
            font-size: 16px;
            margin: 0 0 10px;
        }}
        .nav-list {{
            list-style: none;
            padding-left: 0;
            margin: 0;
        }}
        .nav-list li {{
            margin-bottom: 6px;
        }}
        .nav-list a {{
            display: block;
            padding: 6px 8px;
            border-radius: 6px;
        }}
        .nav-list a:hover {{
            background: #dbe4f3;
        }}
        .content-panelh2 {{
            font-size: 12px;
            margin-left: 8px;
        }}
        .main-content {{
            margin-left: 250px;
            padding: 16px 18px 40px;
        }}
        .card {{
            background: #fff;
            border: 1px solid #e1e4ea;
            border-radius: 10px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.04);
            margin-bottom: 18px;
            padding: 16px 18px;
        }}
        .card h2, .card h1 {{ margin-top: 0; }}
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(160px, 1fr));
            gap: 12px;
            margin: 12px 0;
        }}
        .metric-box {{
            background: #f4f7fb;
            border: 1px solid #dde3ed;
            border-radius: 8px;
            padding: 10px 12px;
            text-align: left;
        }}
        .metric-label {{
            font-size: 12px;
            color: #5b6b7a;
            margin: 0;
        }}
        .metric-value {{
            font-size: 20px;
            font-weight: 700;
            margin: 2px 0 0;
        }}
        table.dataframe {{
            border-collapse: collapse;
            width: 100%;
            font-size: 12px;
        }}
        table.dataframe td, table.dataframe th {{
            border: 1px solid #e1e4ea;
            padding: 6px 8px;
        }}
        table.dataframe tr:nth-child(even) {{ background: #f8f9fb; }}
        table.dataframe th {{
            background: #eef2f8;
            font-weight: 600;
        }}
        .plot-wrapper img {{
            width: 100%;
        }}
    </style>
    </head>
    
    <div class="content-panel_wrapper">
        <div class="content-panel">
            <a href="https://github.com/ArthurVM/Parapipe"><img src="data:image/{prefix};base64,{logo_b64}" width=200></a>
            <hr style="width:80%;text-align:left;"> 
            <h2>Contents</h2>
            <ul class="nav-list">
            {content_links}
            </ul>
        </div>
    </div>

    <div class="main-content">
        <a href="https://github.com/ArthurVM/Parapipe"><img src="data:image/{prefix};base64,{logo_b64}" width=300></a>
        <p>Run initiated at {run_dt} BST</p>
        <p>Report generated at {report_dt} BST</p>

        <div id="summary" class="card">
          <h2>Run Summary</h2>
          <div class="metrics-grid">
            <div class="metric-box">
              <p class="metric-label">Samples</p>
              <p class="metric-value">{len(reports)}</p>
            </div>
            <div class="metric-box">
              <p class="metric-label">Polyclonal gp60</p>
              <p class="metric-value">{0 if gp60_summary is None or gp60_summary.empty else int(gp60_summary['polyclonal'].sum())}</p>
            </div>
            <div class="metric-box">
              <p class="metric-label">MOI results</p>
              <p class="metric-value">{sum(1 for r in reports.values() if r.get('heterozygosity'))}</p>
            </div>
            <div class="metric-box">
              <p class="metric-label">Total SNPs</p>
              <p class="metric-value">{int(sum(d['mapping']['allele_report']['total_snps'] for d in reports.values()))}</p>
            </div>
          </div>
          <p>Input Parameters</p>
          <div>{params_line}</div>
          <p>MultiQC report: <a href={qc_report[0]}>{qc_report[0]}</a></p>
        </div>

        <div id="map" class="card">
        <h2>Mapping Statistics</h2>
        <p>Mapping statistics for samples within this Parapipe run.</p>
        {map_tab_fig.to_html(full_html=False)}
        </div>
        
        <div id="gini" class="card">
        <h2>Gini-Granularity Curves</h2>
        <p>GG curves characterise coverage inequality. Low GG-area indicates low coverage equality. Low normalised GG-area indicates reads aggregate in small\
            regions, whereas high normalised GG-area indicates reads aggregate in larger regions <a href="https://link.springer.com/chapter/10.1007/978-3-030-46970-2_11">[1]</a>.</p>
        {gg_fig.to_html(full_html=False)}
        </div>

        <div id="phylo" class="card">
        <h2>Phylogenetic Analysis</h2>
        <p>All analysis carried out using absolute <a href="https://en.wikipedia.org/wiki/Hamming_distance">(Hamming)</a> SNP distance across the whole genome.</p>
        <p>Key:</p>
        <p>red = samples processed this run.</p>
        <p>grey = samples included in the input archive.</p>
        
        {phylo_plots_wrapper}
        </div>

        <div id="snp" class="card">
        <h2>SNP Report</h2>
        {stacked_bar.to_html(full_html=False)}
        </div>

        <div id="gp60" class="card">
        <h2>gp60 Subtyping</h2>
        <p>Subtype proportions stacked per sample; colors indicate gp60 family. Polyclonal calls are flagged in the table.</p>
        {"<p>No gp60 results available.</p>" if gp60_summary is None or gp60_summary.empty else ""}
        {"<div>"+gp60_fig.to_html(full_html=False)+"</div>" if gp60_fig else ""}
        {"<div>"+gp60_summary.assign(polyclonal=gp60_summary['polyclonal'].map({True:'<span style=\"font-weight:bold;color:#b22222\">Yes</span>',False:'No'})).to_html(index=False, escape=False)+"</div>" if gp60_summary is not None and not gp60_summary.empty else ""}
        </div>

        <div id="moi" class="card">
        <h2>MOI Analysis</h2>
        {moi_htmls}
        </div>

    </div>
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
        
        if prefix not in reports:
            continue
        
        het_data = parse_json(rep)
        maf_plot = het_data.get("maf_plot")
        if maf_plot:
            het_data["maf_plot_path"] = os.path.join(os.path.abspath(os.path.dirname(rep)), maf_plot)
        reports[prefix]["heterozygosity"] = het_data

    map_tab = make_map_table(reports)
    map_df = pd.DataFrame(map_tab[1:], columns=map_tab[0])

    map_df.to_csv("./run_results.csv")

    gg_dict = { k : {int(w) : g for w, g in v["mapping"]["mapping_stats"]["gg_array"].items()} for k, v in reports.items() }
    gg_fig = make_gg_plot(gg_dict, fig_height)

    gp60_summary, gp60_stack = summarise_gp60(reports)
    gp60_fig = make_gp60_bar(gp60_stack)

    ## construct html
    write_html(reports, map_df, gg_fig, args.multiqc, env, args.dist_matrix, run_samples, gp60_summary, gp60_fig)

main(sys.argv)
