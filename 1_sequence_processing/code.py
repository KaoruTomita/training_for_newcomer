from typing import List, Union
import numpy.typing as npt
import numpy as np

def base_count(fastafile: str) -> List[int]:
    # 課題 1-1
    ans=[0,0,0,0]
    with open(fastafile, "r") as f:
        name = f.readline()
        while True:
            line=f.readline()
            for i in range(len(line)):
                if line[i]=="A":
                    ans[0]+=1
                elif line[i]=="T":
                    ans[1]+=1
                elif line[i]=="G":
                    ans[2]+=1
                elif line[i]=="C":
                    ans[3]+=1
            if not line:
                break
    
    return ans # A, T, G, C

def gen_rev_comp_seq(fastafile: str) -> str:
    # 課題 1-2
    ans=""
    with open(fastafile, "r") as f:
        name = f.readline()
        while True:
            line=f.readline()
            for i in range(len(line)):
                if line[i]=="A":
                    ans+="T"
                elif line[i]=="T":
                    ans+="A"
                elif line[i]=="G":
                    ans+="C"
                elif line[i]=="C":
                    ans+="G"
            if not line:
                break
    ans=ans[::-1]
    return ans

def calc_gc_content(fastafile: str, window: int=1000, step: int=300) -> Union[npt.NDArray[np.float_], List[float]]:
    # 課題 1-3
    # 値を出力するところまで。matplotlibを使う部分は別途実装してください。
    ans=[]
    dna=[]
    with open(fastafile, "r") as f:
        name = f.readline()
        lines=[line.strip() for line in f.readlines()]
        dna=list("".join(lines))
    for i in range(0,len(dna)-window+1,step):
        gc=0
        for j in range(window):
            if dna[i+j]=="G" or dna[i+j]=="C":
                gc+=1
        ans.append(gc/window*100)
    return ans

def search_motif(fastafile: str, motif: str) -> List[str]:
    # 課題 1-4
    with open(fastafile, "r") as f:
        name = f.readline()
        lines=[line.strip() for line in f.readlines()]
        dna=list("".join(lines))
    rev=list(gen_rev_comp_seq(fastafile))
    ans=[]
    for i in range(len(dna)-len(motif)+1):
        if dna[i:i+len(motif)]==list(motif):
            ans.append("F"+str(i+1))
    for i in range(len(rev)-len(motif)+1):
        if rev[i:i+len(motif)]==list(motif):
            ans.append("R"+str(len(dna)-i))
    return ans

def translate(fastafile: str) -> List[str]:
    # 課題 1-5
    with open(fastafile, "r") as f:
        name = f.readline()
        lines=[line.strip() for line in f.readlines()]
        dna=list("".join(lines))
    rev=list(gen_rev_comp_seq(fastafile))
    ans=[]
    codon_table = {
        "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
        "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
        "UAU": "Y", "UAC": "Y", "UAA": "_", "UAG": "_",
        "UGU": "C", "UGC": "C", "UGA": "_", "UGG": "W",
        "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
        "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
        "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
        "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }
    tmp=""
    i=0
    while i<len(dna)-2:
        codon="".join(dna[i:i+3]).replace("T","U")
        if codon_table[codon]=="M":
            tmp+="M"
            i+=3
        elif len(tmp)>0:
            tmp+=codon_table[codon]
            if codon_table[codon]=="_":
                ans.append(tmp)
                tmp=""
            i+=3
        else :
            i+=1
    if len(tmp)>0:
        ans.append(tmp)
    tmp=""
    i=0
    while i<len(rev)-2:
        codon="".join(rev[i:i+3]).replace("T","U")
        if codon_table[codon]=="M":
            tmp+="M"
            i+=3
        elif len(tmp)>0:
            tmp+=codon_table[codon]
            if codon_table[codon]=="_":
                ans.append(tmp)
                tmp=""
            i+=3
        else :
            i+=1
    if len(tmp)>0:
            ans.append(tmp)

    return ans

if __name__ == "__main__":
    filepath = "data/NT_113952.1.fasta"
    # 課題 1-1
    print(base_count(filepath))
    # 課題 1-2
    print(gen_rev_comp_seq(filepath))
    # 課題 1-3
    print(calc_gc_content(filepath))
    # 課題 1-4
    print(search_motif(filepath, "ATG"))
    # 課題 1-5
    print(translate(filepath))
