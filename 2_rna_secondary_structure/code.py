from typing import List, Tuple, Union
import numpy.typing as npt
import numpy as np

def enumerate_pairs(fastafile: str) -> List[Tuple[int, int]]:
    # 課題 2-1
    ans=[]
    with open(fastafile, 'r') as f:
        name= f.readline()
        lines=[line.strip() for line in f.readlines()]
        dna=list("".join(lines))
    A=[]
    U=[]
    G=[]
    C=[]
    for i in range(len(dna)):
        if dna[i]=="A":
            A.append(i+1)
        elif dna[i]=="T":
            U.append(i+1)
        elif dna[i]=="G":
            G.append(i+1)
        elif dna[i]=="C":
            C.append(i+1)
    for a in A:
        for u in U:
            if a<u:
                ans.append((a,u))
            else :
                ans.append((u,a))
    for g in G:
        for c in C:
            if g<c:
                ans.append((g,c))
            else :
                ans.append((c,g))
    return ans

def enumerate_possible_pairs(fastafile: str, min_distance: int=4) -> List[Tuple[int, int]]:
    # 課題 2-2
    ans=enumerate_pairs(fastafile)
    ans=[pair for pair in ans if pair[1] - pair[0] >= min_distance]
    return ans

def enumerate_continuous_pairs(fastafile: str, min_distance: int=4, min_length: int=2) -> List[Tuple[int, int, int]]:
    # 課題 2-3
    possible_parirs=enumerate_possible_pairs(fastafile, min_distance)
    possible_parirs=sorted(possible_parirs, key=lambda x: x[0],reverse=True)
    ans=[]
    with open(fastafile, 'r') as f:
        name= f.readline()
        lines=[line.strip() for line in f.readlines()]
        dna=list("".join(lines))
    dict={"A":"T","T":"A","G":"C","C":"G"}
    for l,r in possible_parirs:
        count=0
        while l>=1 and r<=len(dna):
            if dict[dna[l-1]]!=dna[r-1]:
                break
            l-=1
            r+=1
            count+=1
        l+=1
        r-=1
        if count>= min_length:
            ans.append((l,r,count))
    return ans

def create_dotbracket_notation(fastafile: str, min_distance: int=4, min_length: int=2) -> str:
    # 課題 2-4
    used=set()
    continuous_pairs=enumerate_continuous_pairs(fastafile, min_distance, min_length)
    with open(fastafile, 'r') as f:
        name= f.readline()
        lines=[line.strip() for line in f.readlines()]
        dna=list("".join(lines))
    ans=["."]*len(dna)
    for l,r,count in continuous_pairs:
        flg=True
        for i in range(count):
            if l+i-1 in used or r-i-1 in used:
                flg=False
                break
        if not flg:
            continue
        for i in range(count):
            ans[l+i-1]="("
            ans[r-i-1]=")"
            used.add(l+i-1)
            used.add(r-i-1)
    return "".join(ans)

if __name__ == "__main__":
    filepath = "data/AUCGCCAU.fasta"
    # 課題 2-1
    print(enumerate_pairs(filepath))
    # 課題 2-2
    print(enumerate_possible_pairs(filepath))
    # 課題 2-3
    print(enumerate_continuous_pairs(filepath, 2))
    # 課題 2-4
    print(create_dotbracket_notation(filepath, 2))

    #2-5
    print(create_dotbracket_notation("data/NM_014495.4.fasta", 2))


