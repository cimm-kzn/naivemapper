import argparse


def parse_args():
    parser = argparse.ArgumentParser(description="", epilog="", prog='naivemapper')
    parser.add_argument("--version", "-v", action="version", version="0.1", default=False)
    parser.add_argument("--input", "-i", default="input.rdf", type=str, help="RDF inputfile")
    parser.add_argument("--output", "-o", default="output.rdf", type=str, help="RDF outputfile")
    parser.add_argument("--output_rank", "-o2", default="rank/out_rank.txt", type=str,
                        help="the .txt-file with probability and mapping values ")
    parser.add_argument("--model", "-n", default="model.dat", type=str, help="Model file")
    parser.add_argument("--stage", "-s", type=int, default=2,
                        help="type of stage of the process: 0-learning, 1-prediction, 2-cross_validation")
    parser.add_argument("--f_type", "-ft", type=int, default=0, help="Method of fragmentation of a molecule")
    parser.add_argument("--min", "-m", type=int, default=1, help="minimal fragments length")
    parser.add_argument("--max", "-M", type=int, default=8, help="maximal fragments length")
    parser.add_argument("--deep", "-deep", type=int, default=3, help="maximal augmented-fragments length")
    parser.add_argument("--f_count", "-fc", type=int, default=0,
                        help="The way of accounting information on the number of fragments of each type in the set")
    parser.add_argument("--pairs", "-p", type=int, default=1,
                        help="type of united respective pairs: 0 - simple, 1 - equivalent")
    parser.add_argument("--duplicate", "-d", type=int, default=0,
                        help="duplicate availability:0-does all duplicate,1-has all duplicate,2-does 'False' duplicate")
    parser.add_argument("--bitstring", "-b", type=int, default=2,
                        help="type of united bitstring for atomic bitstrings A (substrats) and B (products): "
                             "0 - [A*B], 1 - [A+B+(A*B)], 2 - [(A!B)+(B!A)+(A*B)], 3 - [A+B], 4 - [(AxorB)+(A*B)], "
                             "5-Based on the pairwise (name_frag_in_sub_set, name_frag_in_prod_set)")
    parser.add_argument("--b_length", "-l", type=int, default=2048, help="lenght of bitstrings")
    parser.add_argument("--folds", "-N", type=int, default=5,
                        help="the number of partitions of the data: sets to training and test")
    parser.add_argument("--repeat", "-r", type=int, default=1,
                        help="the number of repetitions of the cross-validation procedure")
    parser.add_argument("--dfs", "-dfs", type=int, default=0,
                        help="Choice of method for reconsidering the assignment: "
                             "0-for symmetrically equivalent groups, 1-for the values of probabilities")
    parser.add_argument("--debug", action='store_true', help="debug mod")
    return parser
