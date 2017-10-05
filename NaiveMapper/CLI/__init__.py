from .parser import parse_args
from ..core import mapper_core


def main():
    parser = parse_args()
    args = parser.parse_args()
    mapper_core(**vars(args))
