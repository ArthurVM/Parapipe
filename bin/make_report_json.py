import json
import os, sys
import pandas as pd
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


def main(sample_id, mapstats_json, typing_json, al_stats_json):

    prefix = os.path.basename(mapstats_json).split("_mapstats.json")[0]

    with open(al_stats_json, 'r') as fin:
        al_data = json.load(fin)[prefix]

    process_json(sample_id, mapstats_json, typing_json, al_data)


if __name__=="__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
