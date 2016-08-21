#!/opt/anaconda3/bin/python3.4
import sys
import re

if (len(sys.argv) != 2):
    print("exec <phylip>")
    sys.exit(0)

file  = sys.argv[1]

def compute_CpG(puz):
    f = open(puz,"r")
    lines = f.readlines()
    f.close()
    lines_clean = []
    Nsite_nuc = 0
    for line,index in zip(lines,[i for i in range(0,len(lines))]):
        if(index != 0):
            lines_clean.append(line.split(" ")[-1])
        else:
            Nsite_nuc = line.split(" ")[-1].strip()

    Nsp = cpg12 = cpg23 = cpg31 = 0
    for line_clean in lines_clean:
        line_clean = line_clean.strip()
        Nsp += 1

        for i,j in zip([line_clean[i] for i in range(0,len(line_clean),3)],\
                       [line_clean[i] for i in range(1,len(line_clean),3)]):
            if ((i=="C" or i=="c") and (j=="G" or j == "g")):
                cpg12 +=1
        for i,j in zip([line_clean[i] for i in range(1,len(line_clean),3)],\
                       [line_clean[i] for i in range(2,len(line_clean),3)]):
            if ((i=="C" or i == "c") and (j=="G" or j=="g")):
                cpg23 +=1
        for i,j in zip([line_clean[i] for i in range(2,len(line_clean),3)],\
                       [line_clean[i] for i in range(3,len(line_clean),3)]):
            if ((i=="C" or i == "c") and (j=="G" or j=="g")):
                cpg31 +=1

    return (cpg12+cpg23+cpg31,cpg12,cpg23,cpg31,Nsp,Nsite_nuc)

tup = compute_CpG(file)
print("Ntaxa ", tup[4], " Nsite_nuc ", tup[5])
print("CpG ", round(tup[0],3))
print("CpG1 ", round(tup[1],3))
print("CpG2 ", round(tup[2],3))
print("CpG3 ", round(tup[3],3))
