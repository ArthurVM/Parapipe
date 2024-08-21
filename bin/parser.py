import sys
import argparse
from ._version import __version__

def run(args):
    """ Level 0 run function
    """
    None

base_parser = argparse.ArgumentParser(add_help=True, description=__doc__)

base_parser.add_argument('-v', '--version', action='version',
    version='%(prog)s {version}'.format(version=__version__))

subparsers = base_parser.add_subparsers(
    title="[sub-commands]", dest="command"
)

## run pipeline args parser
parser_runPipeline = subparsers.add_parser(
    "run_pipeline",
    help="Run the Parapipe pipeline.",
)

parser_runPipeline.add_argument('ID_file',
    type=isFile,
    action='store',
    help='List of line seperated IDs to download. By default, these will be assumed to be species names. If the -a flag is provided, they will be assumed as accession IDs.')

parser_runPipeline.add_argument('-a', '--accessions',
    action='store_true',
    default=False,
    help='Flag to specify that ID_file contains accession IDs, rather than species names.')

parser_runPipeline.add_argument('-n', '--num_assemblies',
    type=int,
    default=1,
    action='store',
    help='The number of assemblies for each species to download. If the -a flag is provided, this is ignored. Default=1.')

parser_runPipeline.add_argument('-o', '--output_prefix',
    type=str,
    default='assemblies',
    action='store',
    help='Output prefix for this run. Default=assemblies.')


## container setup
parser_runPipeline = subparsers.add_parser(
    "setup",
    help="Build singularity containers. NOTE: This will overwrite existing containers.",
)