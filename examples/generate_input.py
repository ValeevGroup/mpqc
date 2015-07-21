#!/usr/bin/python

import json, sys, getopt
from pprint import pprint

def main(argv):

    try:
        opts, args = getopt.getopt(argv,"h",["json=","xyz=","nobs=","ndfbs=","nocc=","obs=","dfbs=","sthresh=","lthresh="])
    except getopt.GetoptError:
        print 'random_scan.py'
        print '--json <json template file>'
        print '--xyz <molecule>'
        print '--nobs <num obs clusters>'
        print '--ndfbs <num dfbs clusters>'
        print '--nocc <num occupied clusters>'
        print '--obs <basis>'
        print '--dfbs <df basis>'
        print '--sthresh <block sparse threshold>'
        print '--lthresh <low rank threshold>'
        sys.exit(2)


    json_file=''

    for opt, arg in opts:
        if opt == '-h':
            print 'random_scan.py'
            print '--json <json template file>'
            print '--xyz <molecule>'
            print '--nobs <num obs clusters>'
            print '--ndfbs <num dfbs clusters>'
            print '--nocc <num occupied clusters>'
            print '--obs <basis>'
            print '--dfbs <df basis>'
            print '--sthresh <block sparse threshold>'
            print '--lthresh <low rank threshold>'
            sys.exit()
        elif opt == '--json':
            json_file=arg

    with open(json_file) as input_json:
        json_data = json.load(input_json)


    for opt, arg in opts:
        if opt == "--xyz":
            json_data['xyz file']=str(arg)
        elif opt == "--nobs":
            json_data["number of bs clusters"]=int(arg)
        elif opt == "--ndfbs":
            json_data["number of dfbs clusters"]=int(arg)
        elif opt == "--nocc":
            json_data["number of occupied clusters"]=int(arg)
        elif opt == "--obs":
            json_data["basis"]=str(arg)
        elif opt == "--dfbs":
            json_data["df basis"]=str(arg)
        elif opt == "--sthresh":
            json_data["block sparse threshold"]=float(arg)
        elif opt == "--lthresh":
            json_data["low rank threshold"]=float(arg)
    
    print json.dumps(json_data)


if __name__ == "__main__":
    main(sys.argv[1:])
