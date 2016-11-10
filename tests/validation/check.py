#
# compares MPQC4 outputs
# usage: check.py <output> <reference_output>
#

##########################################################
# util
##########################################################

import sys, re, math

def validate(label, data, refdata, tolerance):
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

def total_energy(file_name):
    file = open(file_name, 'r')
    for line in file:
        match1 = re.match('Wfn energy is: ' + pat_numbers(1), line)
        if match1:
            return match1.groups()

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
