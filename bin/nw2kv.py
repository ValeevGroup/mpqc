#!/usr/bin/env python
from __future__ import print_function, division
from os.path import basename, splitext
import re
import sys

infile = sys.argv[1]
opt_gen_contraction = True

class Shell(object):
    def __init__(self, am):
        self.am = am
        self.coefs = []
        self.expons = []

with open(infile) as f:
    data = f.read()

    bsname = data.splitlines()[0].rstrip()
    mpqc_filename = sys.argv[2] if len(sys.argv) > 2 else bsname.lower() + ".kv"

    print('Writing basis "{}" from {} to file {}'.format(bsname, infile, mpqc_filename))

    with open(mpqc_filename, "w+") as o:

        o.write("basis:(\n")

        for m in re.finditer(
                r"""(?:
                    ^([A-Z]+)         # Element name
                    \s+
                    (?:
                        [SPDFGHIKLMNORSTUVWXYZ]\s+\d+  # angular momentum
                        (?:\s+\d+\s+-?\d+\.\d+\s+-?\d+\.\d+)+\s+ # data
                    )+)
                    |^!(.*)$
                    """,
                data,
                flags=re.MULTILINE|re.VERBOSE
        ):
            if m.group(2) is not None:
                o.write("%{}\n".format(m.group(2)))
            else:
                element_name = m.group(1)
                shells = []
                for mm in re.finditer(r"""
                    ([SPDFGHIKLMNORSTUVWXYZ])\s+(\d+)  # angular momentum
                    ((?:\s+\d+\s+-?\d+\.\d+\s+-?\d+\.\d+e?-?\d*)+) # data
                """, m.group(0), flags=re.MULTILINE|re.VERBOSE):
                    am = mm.group(1)
                    coefexpon = mm.group(3)
                    shell = Shell(am)
                    for line in coefexpon.strip().splitlines():
                        _, expon, coef = line.split()
                        shell.expons.append(expon)
                        shell.coefs.append(coef)
                    if len(shell.expons) != int(mm.group(2)):
                        raise ValueError("Invalid NWChem file format.  Expected {} coefficients and exponents, got {}")
                    shells.append(shell)
                grouped_shells = [[shells[0]]]
                spot = 1
                while spot < len(shells):
                    if shells[spot].expons == grouped_shells[-1][0].expons:
                        grouped_shells[-1].append(shells[spot])
                    else:
                        grouped_shells.append([shells[spot]])
                    spot += 1
                o.write('{}: "{}": [\n'.format(element_name.lower(), bsname))
                for shgrp in grouped_shells:
                    o.write('  (type: [')
                    for sh in shgrp:
                        if sh.am in "SP":
                            o.write(" am = {}".format(sh.am.lower()))
                        else:
                            o.write(" (am = {} puream = 1)".format(sh.am.lower()))
                    o.write(" ]\n    {exp")
                    for i, _ in enumerate(shgrp):
                        o.write(" coef:{}".format(i))
                    o.write("} = {\n      ")
                    for i in range(len(shgrp[0].expons)):
                        o.write(shgrp[0].expons[i])
                        for sh in shgrp:
                            o.write("   " + sh.coefs[i])
                        o.write("\n      ")
                    o.write("    })\n")
                o.write("]\n")

        o.write(")")














