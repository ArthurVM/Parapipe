import json
import os, sys
from collections import defaultdict
from matplotlib import pyplot as plt


def readJSON(json_file):
    with open(json_file, 'r') as j:
        return json.loads(j.read())


def process_json(sample_id, mapstats_json, typing_json, al_data):

    report_json = defaultdict(lambda: defaultdict(dict))

    report_json["mapping"] = readJSON(mapstats_json)

    ## check the Molecular typing scheme was run
    if typing_json == "None":
        ## else initialise empty fields
        report_json["typing"] = None
    
    else:
        typing_data = readJSON(typing_json)

        for scheme, subdict in typing_data.items():
            tmp = subdict["locus_coverage"][sample_id]
            report_json["typing"][scheme]["locus_coverage"] = tmp

            report_json["typing"][scheme]["phylo"] = {}
            report_json["typing"][scheme]["phylo"]["contree"] = subdict["contree"]
            report_json["typing"][scheme]["phylo"]["mldist"] = subdict["mldist"]

    report_json["mapping"]["allele_report"] = al_data
    json.dump(report_json, indent=4, fp=sys.stdout)
    
    return report_json


def make_plots(cov_array, gg_array):
    plt.style.use('ggplot')

    data = [b for c, b in cov_array.items()]
    x_pos = [i for i, _ in enumerate(data)]
    x_labels = [c for c, b in cov_array.items()]

    plt.bar(x_pos, data, color='blue' )
    plt.xticks(x_pos, x_labels)
    plt.xlabel("Depth of coverage")
    plt.ylabel("% Genome covered")
    plt.title("Genome covered to depth")

    plt.tight_layout()

    plt.savefig("./covplot.jpeg")


def makehtml(map_stats, phylo_tree, sampleID):
    html = f"""
    <html>
        <head>
            <title> Parapipe Report </title>
        </head>
        <body>
            <p style="font-size:30px">{sampleID} Parapipe Report </p>
            <p>Reads mapped: {map_stats['reads_mapped']}</p>
            <p>GG area: {map_stats['GG_area']}</p>
            <p>nGG area: {map_stats['nGG_area']}</p>
            <img src='covplot.jpeg' width="400">
            <img src='{phylo_tree}' width="600">
        </body>
    </html>
    """

    with open("./report.html", "w") as htmlout:
        htmlout.write(html)


def main(sample_id, mapstats_json, typing_json, al_stats_json):

    prefix = os.path.basename(mapstats_json).split("_mapstats.json")[0]

    with open(al_stats_json, 'r') as fin:
        al_data = json.load(fin)[prefix]

    process_json(sample_id, mapstats_json, typing_json, al_data)
    # make_plots(map_stats['coverage_breadth'], map_stats['gg_array'])
    #
    # makehtml(map_stats, phylo_tree, os.path.splitext(os.path.basename(json_file))[0])


if __name__=="__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
