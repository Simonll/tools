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
    
inputfile  = sys.argv[1]
index  = int(sys.argv[2])
outputfile = sys.argv[3]
param = sys.argv[4]
Nparam = len(sys.argv)

def readChainMSAA(pathtochain,index):
   
    chaintoread = open(pathtochain, "r")
    lines = chaintoread.readlines()
    chaintoread.close()
    print("path to chain " , pathtochain, " iteration to recover ", index, "number of lines in chain file ",  len(lines))
    
    pt= 0 
    k = 0
    brack = 0
    start = 0 
    end = 0
    a = True
    b = 0
    new_start = 0
    old_start = 0 

    for line in lines:
        if (line.startswith("(") and brack == 0 and pt == index):
            start = k
            brack = 1
            
        if (line.startswith("(") and brack == 1 and pt == index+1):
            end = k
            
        if end != 0 : 
            break
            
        if (line.startswith("(")):
            pt+= 1
            
        if (line.startswith("(") and k != 0 and a == True):
            step = k
            a = False 
            print("expected step size (base on first iteration) : ", step)
            print("only inconsitent iterations will be shown")
        
        if (line.startswith("(") and pt > 1):
            old_start = new_start
            new_start = k
            if (step != (new_start-old_start)):
                print("inconsitent chain ","pt ",pt," expected step size ", step," new step size " ,new_start-old_start, " start ", old_start, " end ",new_start)
                
            
        k+=1
    if(pt < index):
        print("faulty chain")
        sys.exit(0)
    
    lines = lines[start:end]
    
#     ssaaprofiles = None #np.matrix
#     allocations = None #np.array
    
    tree = lines[0]
    tl = extrac_treelenght(tree)
    
#     nucp = np.array(lines[3].strip().split("\t"),dtype=float)
    
#     nucrr = np.array(lines[5].strip().split("\t"),dtype=float)
    
#     codonprofile = np.array(lines[7].strip().split("\t"),dtype=float)
    
#     omega = float(lines[9])
    
#     Nsite_codon = int(lines[11])
    
    
#     ssaaprofiles_list = [line.strip().split("\t") for line in lines[14:Nsite_codon+14]]
#     ssaaprofiles = np.array(ssaaprofiles_list,dtype=float)
    
#     allocations =  np.array(lines[Nsite_codon+14].strip().split("\t"),dtype=float)
    
    
    dic = {} 
    param = ["tree","tl","nucp","nucrr","codonprofile","omega","Nsite_nuc","ssaaprofiles","allocations","iter"]
    dic["tree"] = tree
    dic["tl"] = tl
    dic["iter"] = lines
#     dic["nucp"] = nucp
#     dic["nucrr"] = nucrr
#     dic["codonprofile"] = codonprofile
#     dic["omega"] = omega
#     dic["Nsite_codon"] = Nsite_codon
#     dic["ssaaprofiles"] = ssaaprofiles
#     dic["allocations"] = allocations
    return dic

def extrac_treelenght(tree):
        t = Tree(tree)
        s = np.array([node.dist for node in t.traverse("postorder") ])
        s = s.sum()
        return s 
        
if __name__ == '__main__':
    
    
    pathtochain = inputfile + ".chain"
    print(pathtochain,index)                       
    dic = readChainMSAA(pathtochain,index)
    
    outfile =  open(outputfile,"w")
    if (param in dic):
        if (param == "tree" or param == "tl" ): 
            outfile.writelines(str(dic[param])+"\n")
        if (param == "iter"):
            outfile.writelines(dic[param])
    else : 
        print("unknown param : ", param)
    
    
    
    