#!/opt/anaconda3/bin/python3.4

import numpy as np
import pandas as pd
import sys
import os
import csv
import re
import math
from ete3 import Tree

if (len(sys.argv) < 4) :
    print("<inputfile><iterationToRecover><outputfile><param>")
    print("available paramameter to recover")
    print("iter : ", "whole iteration ")
    print("tree : ", "phylogenetic tree ")
    print("tl : "  , "tree length ")
    sys.exit(0)

model = sys.argv[1]
inputfile  = sys.argv[2]
index  = int(sys.argv[3])
outputfile = sys.argv[4]
param = sys.argv[5]
Nparam = len(sys.argv)

def readChainMSAA(pathtochain,index):

    lines = linesParser(pathtochain,index)

    ssaaprofiles = None #np.matrix
    allocations = None #np.array

    tree = lines[0]
    tl = extrac_treelenght(tree)
    L1  = lines[1]
    L2  = lines[2]
    nucp = np.array(lines[3].strip().split("\t"),dtype=float)
    L4  = lines[4]
    nucrr = np.array(lines[5].strip().split("\t"),dtype=float)
    L6  = lines[6]
    codonprofile = np.array(lines[7].strip().split("\t"),dtype=float)
    L8  = lines[8]
    omega = float(lines[9])
    L10 = lines[10]
    Nsite_codon = int(lines[11])
    L12 = np.array(lines[12].strip().split("\t"),dtype=float)
    L13 = lines[13]
    ssaaprofiles_list = [line.strip().split("\t") for line in lines[14:Nsite_codon+14]]
    ssaaprofiles = np.array(ssaaprofiles_list,dtype=float)

    allocations =  np.array(lines[Nsite_codon+14].strip().split("\t"),dtype=float)


    dic = {}
    dic["iter"] = lines
    dic["tl"] = tl
    param = ["tree","tl","nucp","nucrr","codonprofile","omega","Nsite_nuc","ssaaprofiles","allocations","iter"]
    dic["tree"] = tree
    dic["L1"] = L1
    dic["L2"] = L2
    dic["nucp"] = nucp
    dic["L4"] = L4
    dic["nucrr"] = nucrr
    dic["L6"] = L6
    dic["codonprofile"] = codonprofile
    dic["L8"] = L8
    dic["omega"] = omega
    dic["L10"] = L10
    dic["Nsite_codon"] = Nsite_codon
    dic["L12"] = L12
    dic["L13"] = L13
    dic["ssaaprofiles"] = ssaaprofiles
    dic["allocations"] = allocations
    return dic
def  linesParser(pathtochain,index):
    chaintoread = open(pathtochain, "r")
    lines = chaintoread.readlines()
    chaintoread.close()

    print("path to chain " , pathtochain, " iteration to recover ", index, "number of lines in chain file ",  len(lines))
    pt= 1
    k = 1
    brack = 0
    start = 1
    end = 1
    a = True
    #b = 0
    new_start = 1
    old_start = 1
    step =  0

    for line in lines:
        #print(pt)
        if (line.startswith("(") and brack == 0 and pt == index):
            start = k
            brack = 1

        #if (line.startswith("(") and brack == 1 and pt == index+1):
        #    end = k
        if (line.startswith("(") and k > 1 and a == True):
            step = k-1
            a = False
            print("expected step size (base on first iteration) : ", step)
            print("only inconsitent iterations will be shown")

        #if end != 1 :
        #    break

        if (line.startswith("(")):
            pt+= 1

        if (line.startswith("(") and pt > 1):
            old_start = new_start
            new_start = k
            if (step != (new_start-old_start)):
                print("inconsitent chain ","pt ",pt," expected step size ", step," new step size " ,new_start-old_start, " start ", old_start, " end ",new_start)


        k+=1

    if(pt < index):
        print(pt, index)
        print("faulty chain")
        sys.exit(0)

    return lines[start-1:start-1+step]


def readChainMS(pathtochain,index):

    lines = linesParser(pathtochain,index)

    ssaaprofiles = None #np.matrix
    allocations = None #np.array

    tree = lines[0]
    tl = extrac_treelenght(tree)
    L1  = lines[1]
    L2  = lines[2]
    nucp = np.array(lines[3].strip().split("\t"),dtype=float)
    L4  = lines[4]
    nucrr = np.array(lines[5].strip().split("\t"),dtype=float)
    L6  = lines[6]
    codonprofile = np.array(lines[7].strip().split("\t"),dtype=float)
    L8  = lines[8]
    omega = float(lines[9])
    #L10 = lines[10]
    Nsite_codon = int(lines[10])
    L11 = np.array(lines[11].strip().split("\t"),dtype=float)
    L12 = lines[12]
    ssaaprofiles_list = np.array(lines[13].strip().split("\t"),dtype=float)
    allocations =  np.array(lines[14].strip().split("\t"),dtype=float)
    param = ["tree","tl","nucp","nucrr","codonprofile","omega","Nsite_nuc","ssaaprofiles","allocations","iter"]

    dic = {}
    dic["iter"] = lines
    dic["tl"] = tl
    dic["tree"] = tree
    dic["L1"] = L1
    dic["L2"] = L2
    dic["nucp"] = nucp
    dic["L4"] = L4
    dic["nucrr"] = nucrr
    dic["L6"] = L6
    dic["codonprofile"] = codonprofile
    dic["L8"] = L8
    dic["omega"] = omega
    #dic["L10"] = L10
    dic["Nsite_codon"] = Nsite_codon
    dic["L11"] = L11
    dic["L12"] = L12
    dic["ssaaprofiles"] = ssaaprofiles
    dic["allocations"] = allocations
    return dic



def extrac_treelenght(tree):
        t = Tree(tree)
        s = np.array([node.dist for node in t.traverse("postorder") ])
        s = s.sum()
        return s

if __name__ == '__main__':


    print("model ",model)
    print("inputfile ",inputfile)
    print("index ",index)
    print("outputfile ",outputfile)
    print("param ",param)
    print("Nparam ",Nparam)


    pathtochain = inputfile + ".chain"
    print(model,pathtochain,index)
    dic = {}
    if (model == "MSAA") :
        dic = readChainMSAA(pathtochain,index)
    elif (model == "MS") :
        dic = readChainMS(pathtochain,index)

    outfile =  open(outputfile,"a")

    if (param in dic):
        if (param == "tree" or param == "tl" or param == "codonprofile" ):
            outfile.writelines(str(dic[param]))
    elif (param == "extrac"):
        outfile.writelines(dic["iter"])

    else :
        print("unknown param : ", param)
