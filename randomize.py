#!/usr/bin/env python
#=========================================================================
# This is OPEN SOURCE SOFTWARE governed by the GPL version 3.
# Copyright (C)William H. Majoros (bmajoros@alumni.duke.edu)
#=========================================================================
from __future__ import (absolute_import, division, print_function, 
   unicode_literals, generators, nested_scopes, with_statement)
from builtins import (bytes, dict, int, list, object, range, str, ascii,
   chr, hex, input, next, oct, open, pow, round, super, filter, map, zip)
# The above imports should allow this program to run in both Python 2 and
# Python 3.  You might need to update your version of module "future".
import sys
import ProgramName
from Shuffler import Shuffler

NUM_SAMPLES=1000
MAX_DIST=1000000

class Locus:
    def __init__(self,ID,pos):
        self.ID=ID
        self.pos=pos
    def addNeighbor(self,neighbor):
        pass
    def distanceTo(self,other):
        return abs(self.pos-other.pos)

class Gene(Locus):
    def __init__(self,ID,pos):
        Locus.__init__(self,ID,pos)
        self.nearbyEnhancers=[]
        self.regulatedBy=[]
    def addNeighbor(self,neighbor):
        if(isinstance(neighbor,Enhancer)):
            self.nearbyEnhancers.append(neighbor)
    def inDegree(self):
        return len(self.regulatedBy)
        
class Enhancer(Locus):
    def __init__(self,ID,pos):
        Locus.__init__(self,ID,pos)
        self.regulates=[]
    def outDegree(self):
        return len(self.regulates)

def load(filename):
    records=[]
    with open(filename,"rt") as IN:
        IN.readline() # discard header
        for line in IN:
            fields=line.rstrip().split()
            records.append(fields)
    return records

def loadGenes(filename):
    genes=[]
    recs=load(filename)
    for rec in recs:
        if(len(rec)!=3): raise Exception("bad line in file"+filename)
        (Chr,tss,ID)=rec
        tss=int(float(tss))
        gene=Gene(ID,tss)
        genes.append(gene)
    return genes

def loadEnhancers(filename):
    enhancers=[]
    recs=load(filename)
    for rec in recs:
        if(len(rec)!=4): raise Exception("bad line in file"+filename)
        (Chr,begin,end,ID)=rec
        pos=int((int(begin)+int(end))/2)
        enhancer=Enhancer(ID,pos)
        enhancers.append(enhancer)
    return enhancers

def makeHash(recs):
    h={}
    for rec in recs:
        h[rec.ID]=rec
    return h

def loadPairs(filename,geneHash,enhancerHash):
    pairs=load(filename)
    for rec in pairs:
        if(len(rec)!=2): raise Exception("bad line in file"+filename)
        (enhID,geneID)=rec
        enhancer=enhancerHash[enhID]
        gene=geneHash[geneID]
        if(gene.distanceTo(enhancer)<MAX_DIST):
            gene.regulatedBy.append(enhancer)
            enhancer.regulates.append(gene)

def getNeighborhoods(genes,enhancers,MAX_DIST):
    loci=[x for x in genes]
    loci.extend(enhancers)
    loci.sort(key=lambda x: x.pos)
    n=len(loci)
    for i in range(n):
        locus1=loci[i]
        j=i+1
        while(j<n):
            locus2=loci[j]
            if(locus1.distanceTo(locus2)<MAX_DIST):
                locus1.addNeighbor(locus2)
                locus2.addNeighbor(locus1)
                j+=1
            else: break

def getStats(enhancers,OUT):
    for enh in enhancers:
        print(enh.outDegree(),file=OUT)

def sample(genes,enhancers):
    for enh in enhancers: enh.regulates=[]
    for gene in genes:
        inDegree=gene.inDegree()
        gene.regulatedBy=[]
        Shuffler.shuffleArray(gene.nearbyEnhancers)
        for i in range(inDegree):
            enh=gene.nearbyEnhancers[i]
            gene.regulatedBy.append(enh)
            enh.regulates.append(gene)
    #i=50
    #print("gene i regulated by #1:",genes[i].regulatedBy[0].ID)
    #print("gene i regulated by",genes[i].inDegree(),"neighborhood=",
    #      len(genes[i].nearbyEnhancers))
            
#=========================================================================
# main()
#=========================================================================
if(len(sys.argv)!=6):
    exit(ProgramName.get()+" <genes> <enhancers> <enhancer-gene-pairs> <out1> <out2>\n")
(geneFile,enhancerFile,pairFile,outFile1,outFile2)=sys.argv[1:]

genes=loadGenes(geneFile)
enhancers=loadEnhancers(enhancerFile)
print(len(genes),"genes,",len(enhancers),"enhancers")
geneHash=makeHash(genes)
enhancerHash=makeHash(enhancers)
loadPairs(pairFile,geneHash,enhancerHash)
getNeighborhoods(genes,enhancers,MAX_DIST)
#for gene in genes: print(len(gene.regulatedBy),"enhancers")
OUT1=open(outFile1,"wt")
getStats(enhancers,OUT1)
OUT1.close()
OUT2=open(outFile2,"wt")
for i in range(NUM_SAMPLES):
    sample(genes,enhancers)
    getStats(enhancers,OUT2)
#    k=50
#    print("enhancer",k,"regulates:")
#    for gene in enhancers[k].regulates:
#        print("\t",gene.ID,sep="")
OUT2.close()
