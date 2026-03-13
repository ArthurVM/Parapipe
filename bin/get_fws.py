import pandas as pd
import sys
import json, os, shutil, glob
from VCF import *

def main(vcf):
    vcf_obj = VCF(vcf)
    sample_ids = vcf_obj.getSampleIds()
    fws_values, hs_values = vcf_obj.fws()
    print(fws_values)

    df = pd.DataFrame({"fws": fws_values, "hs": hs_values})
    df.index.name = "sample_id"
    df.to_csv("fws.csv")

if __name__=="__main__":
    main(sys.argv[1])
