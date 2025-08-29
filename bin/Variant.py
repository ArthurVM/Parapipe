class Variant(object):
    """class for storing variants defined in a vcf file"""


    def __init__(self, ):
        super(Variant, self).__init__()


    def parseLine(self, line: str, sample_ids: list):
        sline = line.strip("\n").split("\t")
        if len(sline) < 9:
            raise ValueError("Malformed VCF line: fewer than 9 fields")

        chrom, pos, id, ref, alt, qual, filter, info, fmt = sline[:9]

        self.chrom = chrom
        self.pos = int(pos)
        self.id = id
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.filter = filter
        self.info = info
        self.format = fmt

        sample_allele_data = sline[9:]
        sample_data = zip(sample_ids, sample_allele_data)

        self.samples = {}

        for sample_id, allele_data in sample_data:
            ztmp = zip(fmt.split(":"), allele_data.split(":"))
            subdict = {}
            for f, d in ztmp:
                tmp_d = d.split(",")
                if len(tmp_d) == 1:
                    try:
                        d = int(tmp_d[0])
                    except ValueError:
                        d = tmp_d[0]
                else:
                    try:
                        d = [int(x) if x != "." else None for x in tmp_d]
                    except ValueError:
                        d = [x if x != "." else None for x in tmp_d]
                subdict[f] = d

            if "AD" not in subdict:
                raise ValueError("AD tag not found. VCF files must include Allelic Depth (AD) field.")

            self.samples[sample_id] = subdict
