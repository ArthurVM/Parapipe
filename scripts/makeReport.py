import json
import os
from matplotlib import pyplot as plt

def read_json(json_file):
    with open(json_file, 'r') as j:
        report_json = json.loads(j.read())

    map_stats = report_json["bam_stats"]
    return map_stats

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

def main(json_file, phylo_tree):
    map_stats = read_json(json_file)
    make_plots(map_stats['coverage_breadth'], map_stats['gg_array'])

    makehtml(map_stats, phylo_tree, os.path.splitext(os.path.basename(json_file))[0])

if __name__=="__main__":
    # main(sys.argv[1])
    main("/home/amorris/BioInf/Parapipe/report_WD/UKH23_1-1_mapstats.json", "/home/amorris/BioInf/Parapipe/report_WD/UKH23_1-1.sorted.bam.png")
