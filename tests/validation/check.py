#
# compares MPQC4 outputs
# usage: check.py <output> <reference_output>
#

##########################################################
# util
##########################################################

# should work with python 2 and 3
from __future__ import absolute_import, division, print_function, unicode_literals
import sys, re, math
import json

def validate(label, data, refdata, tolerance):
    if not isinstance(data,list):
        data = [data]
    if not isinstance(refdata,list):
        refdata = [refdata]
    ok = True
    ndata = len(refdata)
    for i in range(ndata):
        datum = float(data[i])
        refdatum = float(refdata[i])
        if (math.fabs(refdatum - datum) > tolerance):
            ok = False
            break
    return ok

def pat_numbers(n):
    result = ''
    for i in range(n):
        result += '\s*([+-e\d.]+)'
    return result

def parse_json(file_name):
    match = False
    json_lines = ""
    with open(file_name, 'r') as file:
        for line in file:
            if match:
                json_lines += line
            if not match:
                line = line.decode('utf-8')
                match = line.find("Output KeyVal (format=JSON):") != -1

    result = json.loads(json_lines)
    return result


def total_energy(file_name):
    json = parse_json(file_name)
    return json["property"]["value"]["value"]
    # with open(file_name, 'r') as file:
    #     ifile = iter(file)
    #     for line in ifile:
    #         match1 = re.match('\A\s*Property "Energy" computed with Wavefunction ".*":', line)
    #         if match1:
    #             line = next(ifile)
    #             match1 = re.match('\s*value:\s*' + pat_numbers(1), line)
    #             return match1.groups()

##########################################################
# main
##########################################################
file_name = sys.argv[1]
energy = total_energy(file_name)
ref_file_name = sys.argv[2]
ref_energy = total_energy(ref_file_name)

eok = validate("total energy", energy, ref_energy, 1e-10)

ok = eok
if not ok: sys.exit(1)
