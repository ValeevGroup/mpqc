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

# reload(sys)
# sys.setdefaultencoding("utf-8")

default_precision = {
    "Energy" : 1.0e-9,
    "GFRealPole" : 1.0e-4,
    "ExcitationEnergy" : 1.0e-6
}

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
            print(refdatum)
            print(datum)
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
                match = re.match('^\s*Output KeyVal \(format=JSON\):',line)

    result = json.loads(json_lines)
    return result


def total_energy(json):
    return json["property"]["value"]["value"]


def get_precision(json):
    if "precision" not in json["property"]:
        property = json["property"]["type"]
        return default_precision[property]
    else:
        return float(json["property"]["precision"])
##########################################################
# main
##########################################################
file_name = sys.argv[1]
output_json = parse_json(file_name)
value = total_energy(output_json)

ref_file_name = sys.argv[2]
ref_json = parse_json(ref_file_name)
ref_value = total_energy(ref_json)

# how precise should we expect the results to agree? Depends on precision of both results
precision = max( get_precision(output_json), get_precision(ref_json) )


eok = validate("value", value, ref_value, precision)
ok = eok
if not ok: sys.exit(1)
