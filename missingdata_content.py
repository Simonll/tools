#!/opt/anaconda3/bin/python3.4
import sys
import re

if (len(sys.argv) != 2):
    print("GC_content <phylip>")
    sys.exit(0)

file  = sys.argv[1]

def compute_missingdata(phylip):
    f = open(phylip,"r")
    lines = f.readlines()
    f.close()
    lines_clean = []
    Nsite_nuc = 0
    for line,index in zip(lines,[i for i in range(0,len(lines))]):
        if(index != 0):
            lines_clean.append(line.split(" ")[-1])
        else:
            Nsite_nuc = int(line.split(" ")[-1].strip())

    missing = Nsp = Nchar = 0
    for line_clean in lines_clean:
        line_clean = line_clean.strip()
        Nsp += 1
        lenth = len(line_clean)
        missing += len(re.findall("[Xx\-\?]",str([line_clean[i] for i in range(0,lenth)])))
        Nchar += len(re.findall("[ACGTacgtXx\-\?]",str([line_clean[i] for i in range(0,lenth)])))
    return (missing/(Nsp*Nsite_nuc),Nsp,Nsite_nuc,Nchar)

tup = compute_missingdata(file)
print("Ntaxa ", tup[1], " Nsite_nuc ", tup[2], " Nnuc[ACGTacgt-?] ",tup[3])
print("percent missing data ", round(tup[0],3))