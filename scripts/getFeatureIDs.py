import os
import sys

def getFeatures(gff):
    """ Takes a gff file and yield features
    """
    with open(gff, "r") as fin:
        for line in fin.readlines():
            if line.startswith("#"):
                continue
            else:
                yield line.split("\t")[2]

def main(gff):
    """ main function
    """
    featureset = set([feature for feature in getFeatures(gff)])
    with open("featureset.txt", "w") as fout:
        fout.write("\n".join(featureset))

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("USAGE: getFeatureIDs.py <gff>")
        sys.exit(1)
    else:
        main(sys.argv[1])
