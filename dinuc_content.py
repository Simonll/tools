#!/opt/anaconda3/bin/python3.4
import sys
import re
import numpy as np

if (len(sys.argv) != 2):
    print("exec <phylip>")
    sys.exit(0)

file  = sys.argv[1]

def compute_dinuc(puz):
    
    nuc = {"A" : 0, "C":1, "G":2, "T":3}
    
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

    Nsp = 0
    dinuc = np.zeros((4,4))
    dinuc12 = np.zeros((4,4))
    dinuc23 = np.zeros((4,4))
    dinuc31 = np.zeros((4,4))
    
    for line_clean in lines_clean:
        line_clean = line_clean.strip()
        Nsp += 1

        for i,j in zip([line_clean[i] for i in range(0,len(line_clean),3)],[line_clean[i] for i in range(1,len(line_clean),3)]):
            if (i.upper() in nuc and j.upper() in nuc ):
                dinuc12[nuc[i.upper()],nuc[j.upper()]]+=1 
        
        for i,j in zip([line_clean[i] for i in range(1,len(line_clean),3)],[line_clean[i] for i in range(2,len(line_clean),3)]):
            if (i.upper() in nuc and j.upper() in nuc ):
                dinuc23[nuc[i.upper()],nuc[j.upper()]]+=1 
        
        for i,j in zip([line_clean[i] for i in range(2,len(line_clean),3)],[line_clean[i] for i in range(3,len(line_clean),3)]):
            if (i.upper() in nuc and j.upper() in nuc ):
                dinuc31[nuc[i.upper()],nuc[j.upper()]]+=1 
       
        dinuc = dinuc12 + dinuc23 + dinuc31
        dinuc /= dinuc.sum()
        dinuc12 /= dinuc12.sum()
        dinuc23 /= dinuc23.sum()
        dinuc31 /= dinuc31.sum()
        
    return (dinuc,dinuc12,dinuc23,dinuc31,Nsp,Nsite_nuc)


def printdinuc(dinuc):
    nuc = ["A","C","G","T"]
    for i in range(0,len(nuc)):
        for j in range(0,len(nuc)):
            print(nuc[i],nuc[j],round(dinuc[i,j],3))
        

tup = compute_dinuc(file)
print("Ntaxa ", tup[4], " Nsite_nuc ", tup[5])
print("dinuc")
printdinuc(tup[0])
print("dinuc12")
printdinuc(tup[1])
print("dinuc21")
printdinuc(tup[2])
print("dinuc31")
printdinuc(tup[3])