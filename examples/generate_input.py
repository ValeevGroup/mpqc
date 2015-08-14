#!/usr/bin/python

import json, sys, getopt
from pprint import pprint

def main(argv):

    try:
        opts, args = getopt.getopt(argv,"h",["xyz=","nobs=","ndfbs=","nocc=","obs=","dfbs=","sthresh=","lthresh=","print_clusters=","mp2=","cholesky_vectors=","cluster_occ="])
    except getopt.GetoptError:
        print 'Input Generator'
        print '--xyz <molecule>'
        print '--nobs <num obs clusters>'
        print '--ndfbs <num dfbs clusters>'
        print '--nocc <num occupied clusters>'
        print '--obs <basis>'
        print '--dfbs <df basis>'
        print '--sthresh <block sparse threshold>'
        print '--lthresh <low rank threshold>'
        print '--print_clusters <Prefix for cluster files>'
        print '--mp2 <Do mp2>'
        print '--cholesky_vectors <Use Cholesky Vectors>'
        print '--cluster_occ <Cluster The occupied vectors>'
        sys.exit(2)



    for opt, arg in opts:
        if opt == '-h':
            print 'Input Generator'
            print '--xyz <molecule>'
            print '--nobs <num obs clusters>'
            print '--ndfbs <num dfbs clusters>'
            print '--nocc <num occupied clusters>'
            print '--obs <basis>'
            print '--dfbs <df basis>'
            print '--sthresh <block sparse threshold>'
            print '--lthresh <low rank threshold>'
            print '--print_clusters <Prefix for cluster files>'
            print '--mp2 <T/F>'
            print '--cholesky_vectors <T/F>'
            print '--cluster_occ <T/F>'
            sys.exit()

    json_data = json.loads('{}')
    
    # Define the defaults 
    json_data['xyz file']='file.xyz'
    json_data["number of bs clusters"]=5
    json_data["number of dfbs clusters"]=2
    json_data["number of occupied clusters"]=1
    json_data["basis"]="cc-pvdz"
    json_data["df basis"]="cc-pvdz-ri"
    json_data["block sparse threshold"]=1e-13
    json_data["low rank threshold"]=1e-8
    json_data["print clusters"]=False
    json_data["debug break"]=0
    json_data["do mp2"]=False
    json_data["use cholesky vectors"]=False
    json_data["cluster orbitals"]=False

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
        elif opt == "--print_clusters":
            json_data["print clusters"]=True
            json_data["basis clusters file"]=str(arg)+"_obs.xyz"
            json_data["df basis clusters file"]=str(arg)+"_dfbs.xyz"
        elif opt == "--mp2":
            json_data["do mp2"]=bool(arg)
        elif opt == "--cholesky_vectors":
            json_data["use cholesky vectors"]=bool(arg)
        elif opt == "--cluster_occ":
            json_data["cluster orbitals"]=bool(arg)


    
    print json.dumps(json_data, indent=4, 
            separators=(',', ': '))

if __name__ == "__main__":
    main(sys.argv[1:])
