#!/usr/bin/env python3
"""
Author : oliv123 <oliv123@localhost>
Date   : 2024-03-06
Purpose: a script to return cluster information from a given cluster ID
"""

import argparse
from typing import NamedTuple, TextIO
import pandas as pd
from pathlib import Path
import warnings
from collections import Counter
import ast

class Args(NamedTuple):
    """ Command-line arguments """
    positional: str
    string_arg: str
    int_arg: int
    file: TextIO
    on: bool


# --------------------------------------------------
def get_args() -> Args:
    """ Get command-line arguments """

    parser = argparse.ArgumentParser(
        description='Script to return cluster information for a given cluster ID or IDs from a .txt file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('cluster_input',
                        help='Cluster ID or a path to a .txt file with one ID on each line',
                        type=str)

    return parser.parse_args


# --------------------------------------------------
def main() -> None:
    """ Make a jazz noise here """

    args = get_args()
    str_arg = args.string_arg
    int_arg = args.int_arg
    file_arg = args.file
    flag_arg = args.on
    pos_arg = args.positional

    print(f'str_arg = "{str_arg}"')
    print(f'int_arg = "{int_arg}"')
    print('file_arg = "{}"'.format(file_arg.name if file_arg else ''))
    print(f'flag_arg = "{flag_arg}"')
    print(f'positional = "{pos_arg}"')


# --------------------------------------------------
if __name__ == '__main__':
    main()