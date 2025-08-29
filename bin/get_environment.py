#!/usr/bin/env python3

import re
import json
import os
import glob
import sys
import datetime

def get_versions(singpath):

    versions = {}

    paths = list(set(glob.glob(os.path.join(singpath, "singularity.*"))) - \
                    set(glob.glob(os.path.join(singpath, "*.img"))))

    for s_path in paths:
        container = s_path.split("/")[-1].split(".")[1]
        versions[container] = {}

        with open(s_path, 'r') as fin:
            for line in fin.readlines():
                sline = line.strip("\n").split(" ")[0].split("_")
                if len(sline)<2: continue
                if sline[1].startswith("version="):
                    tool = sline[0]
                    version = sline[1].split("=")[1]
                    versions[container][tool] = version

    return versions


def main(argv):
    env = {}
    
    input_dir, output_dir, ref, yaml, db, singpath, read_n_threshold, mincov, missing, maf, mac = argv

    current_datetime = datetime.datetime.now()
    formatted_date = current_datetime.strftime("%Y-%m-%d")
    formatted_time = current_datetime.strftime("%H:%M:%S")

    env["date"] = formatted_date
    env["time"] = formatted_time

    env["params"] = {}
    params = env["params"]
    params["input_dir"] = input_dir
    params["output_dir"] = output_dir
    params["ref"] = ref
    params["yaml"] = yaml
    params["db"] = db
    params["singpath"] = singpath
    params["read_n_threshold"] = read_n_threshold
    params["mincov"] = mincov
    params["missing"] = missing
    params["maf"] = maf
    params["mac"] = mac

    versions = get_versions(singpath)

    env["versions"] = versions

    with open("env.json", 'w') as fout:
        json.dump(env, fout, indent=4)


if __name__ == "__main__":
    main(sys.argv[1:])
