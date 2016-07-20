#!/opt/anaconda3/bin/python3.4
import sys
import re

if (len(sys.argv) != 2):
    print("GC_content <phylip>")
    sys.exit(0)

file  = sys.argv[1]

def compute_GC(phylip):
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

    gc = gc1 = gc2 = gc3 = acgt = acgt1 = acgt2 = acgt3 = Nsp = 0
    for line_clean in lines_clean:
        line_clean = line_clean.strip()
        Nsp += 1
        lenth = len(line_clean)
        acgt += len(re.findall("[ACGT]",str([line_clean[i] for i in range(0,lenth)])))
        acgt1 += len(re.findall("[ACGT]",str([line_clean[i] for i in range(0,lenth,3)])))
        acgt2 += len(re.findall("[ACGT]",str([line_clean[i] for i in range(1,lenth,3)])))
        acgt3 += len(re.findall("[ACGT]",str([line_clean[i] for i in range(2,lenth,3)])))
        gc += len(re.findall("[GC]",str([line_clean[i] for i in range(0,lenth)])))
        gc1 += len(re.findall("[GC]",str([line_clean[i] for i in range(0,lenth,3)])))
        gc2 += len(re.findall("[GC]",str([line_clean[i] for i in range(1,lenth,3)])))
        gc3 += len(re.findall("[GC]",str([line_clean[i] for i in range(2,lenth,3)])))

    return (gc/acgt,gc1/acgt1,gc2/acgt2,gc3/acgt3,acgt,Nsp,Nsite_nuc)

tup = compute_GC(file)
print("Ntaxa ", tup[5], " Nsite_nuc ", tup[6], " Nnuc[ACGT] ",tup[4])
print("GC ", round(tup[0],3))
print("GC1 ", round(tup[1],3))
print("GC2 ", round(tup[2],3))
print("GC3 ", round(tup[3],3))
