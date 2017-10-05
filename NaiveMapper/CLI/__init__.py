from .parser import argparser


def launcher():
    _argparser = argparser()
    args = _argparser.parse_args()
    if 'func' in args:
        args.func(**vars(args))
    else:
        _argparser.print_help()
