#!/usr/bin/env python
from __future__ import print_function
import os
from os.path import abspath, dirname, join as path_join
import sys
from argparse import ArgumentParser
import re

mydir = dirname(abspath(__file__))
rootdir = dirname(mydir)
default_file_pattern = ".*\.(cc|h)$"

p = ArgumentParser(
        description='Lists all known described classes in the MPQC source.'
    )
p.add_argument(
    'root',
    metavar="MPQC_ROOT",
    nargs='?',
    default=rootdir,
    help="Root directory of the MPQC source, e.g. the directory containing the src directory (default: {0})".format(rootdir)
)

p.add_argument(
    '-p', '--pattern',
    metavar="PAT",
    nargs=1,
    default=[default_file_pattern],
    type=str,
    help="Regexp pattern for source files to search for described classes (default: {0})".format(default_file_pattern)
)

p.add_argument(
    '-x', '--exclude',
    metavar="PAT",
    nargs=1,
    default=[None],
    type=str,
    help="Regexp pattern for source files to exclude when looking for described classes (default: None)"
)

p.add_argument(
    '-d', '--include-duplicates',
    dest="dupes",
    action='store_true',
    default=False,
    help='Include duplicate class names in output (default: False)',
)

p.add_argument(
    '-v', '--verbose',
    action='store_true',
    default=False,
    help='Be verbose',
)
args = p.parse_args()

classes = []

pat = re.compile(r"""
        ClassDesc\s+
        [\w:]+\s*\(\s*
        typeid\s*
        \(\s*\w+\s*\)\s*,\s*"(\w+)"\s*,
    """,
    re.VERBOSE|re.M
)

if args.verbose:
    print("Searching directory tree {0}...".format(args.root))

for root, dirnames, filenames in os.walk(path_join(args.root, "src")):
    file_iter = (
        f for f in filenames if re.search(args.pattern[0], path_join(root, f)) 
                and (args.exclude[0] is None or not re.search(args.exclude[0], path_join(root, f)))
    )
    for filename in file_iter:
        full_file = abspath(path_join(root, filename))
        with open(full_file) as f:
            if args.verbose:
                print("Looking in {0}...".format(full_file))
            data = f.read()
            for m in pat.finditer(data):
                if args.verbose:
                    print('  Found class "{0}"'.format(m.group(1)))
                classes.append(m.group(1))

if args.verbose:
    print("Final sorted list of all classes:")

if args.dupes:
    print(" ".join(sorted(classes)))
else:
    print(" ".join(sorted(list(set(classes)))))
