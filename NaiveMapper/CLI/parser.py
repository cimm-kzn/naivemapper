from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from importlib.util import find_spec
from .main_cross_validation import cross_validation_core
from .main_predict import predict_core
from .main_train import train_core
from ..version import version


def _common(parser):
    parser.add_argument("--pairs", "-p", type=str, default="simple",
                        help="Type of union of atoms pairs:\n"
                             "0 = 'sim' - uniting atoms with the same name (in periodic table),\n"
                             "1 = 'eqv' - uniting same name with atoms symmetries refinement.")
    parser.add_argument("--duplicate", "-d", type=bool, default=True,
                        help="Accounted the atomic pairs information duplicates:\n"
                             "True - doesn't duplicate,\n"
                             "False - does all duplicate.")  # "2-does 'False' duplicate")
    parser.add_argument("--fragment_type", "-ft", type=int, default=2,
                        help="Method of fragmentation of a molecule:\n"
                             "0 = 'seq' - sequenced fragments,\n"
                             "1 = 'aug' - augmented fragments,\n"
                             "2 = 'seq+aug' - sequenced and augmented fragments.")
    parser.add_argument("--min", "-m", type=int, default=1,
                        help="The minimal sequenced fragments length")
    parser.add_argument("--max", "-M", type=int, default=8,
                        help="The maximal sequenced fragments length")
    parser.add_argument("--deep", "-deep", type=int, default=3,
                        help="The maximum number of levels of augmented fragments")
    parser.add_argument("--fragment_count", "-fc", type=bool, default=False,
                        help="Accounted for the number of fragments of each type:\n"
                             "False - to ignored,\n"
                             "True - to account.")
    parser.add_argument("--bitstring", "-b", type=int, default=0,
                        help="Type of union bit-strings the reagents (A) and the products (B):\n"
                             "0 = 'and' - intersection of information [A & B],\n"
                             "1 = [A+B+(A*B)],\n"
                             "2 = [(A!B)+(B!A)+(A*B)],\n"
                             "3 = [A+B],\n"
                             "4 = 'xor' - United 'symmetric difference' and 'intersection' [(A ^ B) + (A & B)],\n"
                             "5 = 'kron' - Tensor product of information.")
    parser.add_argument("--length", "-l", type=int, default=2048,
                        help="Length the bit-strings")
    parser.add_argument("--chunk", "-c", type=int, default=None,
                        help="Necessary partitioning of the process of creating bit strings, "
                             "if the LENGTH value exceeding 100,000")


def train(subparsers):
    parser = subparsers.add_parser('train', help='The stage of the mapping learning on the reaction sets',
                                   formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", "-i", default="input.rdf", type=FileType('r'),
                        help="RDF input file on which to learn to create mapping")
    parser.add_argument("--model", "-n", default="model.dat", type=FileType('wb'),
                        help="File with trained model")
    parser.add_argument("--debug", action='store_true', help="debug mod")
    _common(parser)
    parser.set_defaults(func=train_core)


def predict(subparsers):
    parser = subparsers.add_parser('predict', help='The stage of the mapping prediction on the new reaction  sets',
                                   formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", "-i", default="input.rdf", type=str,
                        help="RDF input filename, to which want to create the mapping")
    parser.add_argument("--model", "-n", default="model.dat", type=FileType('rb'),
                        help="File with trained model")
    parser.add_argument("--output", "-o", default="output.rdf", type=str,
                        help="RDF outputfile")
    parser.add_argument("--dfs", "-dfs", type=int, default=0,
                        help="Choice of the revision method (Depth-first search):\n"
                             "0 - by the symmetrically equivalent groups,\n"
                             "1 - by the values of probabilities.")
    '''
    parser.add_argument("--rank", "-ro", default="rank/rank.txt", type=FileType('w'),
                        help="The debug file with average values of the mapping probability a reaction atoms "
                             "at the mappings value True/False")
    '''
    parser.add_argument("--debug", action='store_true', help="debug mod")
    parser.set_defaults(func=predict_core)


def cross_validation(subparsers):
    parser = subparsers.add_parser('cross_validation', help='The stage of the process cross-validation',
                                   formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", "-i", default="input.rdf", type=str,
                        help="RDF input file")
    parser.add_argument("--fold", "-k", type=int, default=5,
                        help="Split the data into k consecutive folds.")
    parser.add_argument("--repeat", "-r", type=int, default=1,
                        help="The number of repetitions of the cross-validation procedure.")
    parser.add_argument("--dfs", "-dfs", type=int, default=0,
                        help="Choice of the revision method (Depth-first search):\n"
                             "0 - by the symmetrically equivalent groups,\n"
                             "1 - by the values of probabilities.")
    parser.add_argument("--debug", action='store_true', help="debug mod")
    _common(parser)
    parser.set_defaults(func=cross_validation_core)


def argparser():
    parser = ArgumentParser(description="CGRtools", epilog="(c) Dr. Ramil Nugmanov", prog='cgrtools')
    parser.add_argument("--version", "-v", action="version", version=version(), default=False)
    subparsers = parser.add_subparsers(title='subcommands', description='available utilities')

    train(subparsers)
    predict(subparsers)
    cross_validation(subparsers)

    if find_spec('argcomplete'):
        from argcomplete import autocomplete
        autocomplete(parser)

    return parser
