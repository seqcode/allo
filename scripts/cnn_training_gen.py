####Custom training the CNN in for allo input##
#Arg 1 is samfile
#Arg 2 is positive peak set
#Arg 3 is negative peak set
#Arg 4 is output name


import os
import sys
import shutil
import math

#Used to get counts in regions
def getArray(chr, start, stop, genLand):
    array = []
    for k in range (int(start),int(stop+1)):
        key = chr + ";" + str(k)
        if key in genLand:
            array.append(genLand[key])
        else:
            array.append(0)
    return array
    

#Parsing reads and putting UMRs into dictionary
def parseUniq(sam, AS):
    genLandCur = {}
    rBlock = []
    numR = 0
    with open(sam) as f:
        for line in f:
            if not line.strip():
                break
            
            #Ignore header lines
            if line.startswith("@"):
                continue
                
            #Splitting columns
            r = line.split('\t')
                        

            #Appending if the block is empty and going to next line, happens first line
            if len(rBlock) == 0:
                rBlock.append(r)
                numR = numR + 1
                continue
                
            #If Bowtie1 -m 1 -k x was used, reads should already be sorted for truly MM reads so
            #dont need to compare alignment scores
            if AS==0:
                #Adding to block
                leadr = rBlock[0]
                if r[0] == leadr[0]:
                    rBlock.append(r)
                    continue
                #New block
                if not r[0] == leadr[0]:
                    numR = numR + 1
                    #Put uniquely mapped reads into a file
                    if len(rBlock) == 1:
                        genLandCur = readAssign(rBlock, genLandCur)
                    #Creating a new read block for the next read
                    rBlock = []
                    rBlock.append(r)
                
            #If Bowtie2 or BWA is used have to actually look at the alignment scores and compare
            else:
                leadr = rBlock[0]
                try:
                    rScore = int(r[AS].split(':')[2])
                    leadScore = int(leadr[AS].split(':')[2])
                except:
                    ASt = 0
                    if r[2] == "*":
                        rScore = 0
                    else:
                        ASt = 0
                        for t in range(0,len(rBlock[-1])):
                            if rBlock[-1][t].startswith("AS:"):
                                ASt = t
                        rScore = int(rBlock[-1][ASt].split(':')[2])
                    
                    if leadr[2] == "*":
                        leadScore = 0
                    else:
                        ASt = 0
                        for t in range(0,len(leadr)):
                            if leadr[t].startswith("AS:"):
                                ASt = t
                        leadScore = int(leadr[ASt].split(':')[2])
                    
                
                #Adding to block
                if r[0] == leadr[0] and rScore == leadScore:
                    rBlock.append(r)
                    continue
                
                #Deleting old block if read found with better score
                if r[0] == leadr[0] and rScore > leadScore:
                    rBlock = []
                    rBlock.append(r)
                    continue
                
                #New block
                if not r[0] == leadr[0]:
                    numR = numR + 1
                    #Put uniquely mapped reads into a file and into dictionary
                    if len(rBlock) == 1:
                        genLandCur = readAssign(rBlock, genLandCur)
                    #Creating a new read block for the next read
                    rBlock = []
                    rBlock.append(r)
    
    return genLandCur



#Assign reads (straight to dictionary for uniq and actual assign for multi-mapped)
def readAssign(rBlock, genLand):
    
    ##Uniquely mapped reads##
    if len(rBlock) == 1:
        #Adding to genetic landscape
        key = rBlock[0][2] + ";" + str(int(rBlock[0][3]))
        if key in genLand:
            genLand[key] = genLand[key] + 1
        else:
            genLand[key] = 1
        return genLand

    
    return genLand


##Main Method##          
if __name__ == '__main__':

    #Samfile
    samfile = sys.argv[1]   
    pos = sys.argv[2]
    neg = sys.argv[3]
    counts = open(sys.argv[4], "a") #file to write to
    
    ####Creating dictionary with uniquely mapped reads###
    #Getting column where alignment score is stored in sam (different for every aligner). If none listed assumes bowtie1 -m 1 -k x
    with open(samfile) as f:
        for line in f:
            if not line.startswith("@") and not line[2] == "*":
                l = line
                break
    AS = 0
    l = l.split('\t')
    for i in range(0,len(l)):
        if l[i].startswith("AS:"):
            AS = i
    #Parsing unique reads and putting into dictionary
    genLand = parseUniq(samfile, AS)
    
    
    
    ###Getting counts at given peak regions###
    n = 0
    with open(pos) as f:
        for line in f:
            r = line.split('\t')
            chr = r[0]
            start = int(r[1])
            stop = int(r[2])
            mid = (start + stop) / 2
            arr = getArray(chr,mid-250,mid+250,genLand)
            counts.write(','.join(map(str,arr)) + '\t1\n')
            n = n + 1 #Used to get equal amounts of negative regions for a balanced dataset
            
            
    #Limiting number of zero regions to 30% in negative set to allow CNN to learn more patterns
    n_zero = int(n*0.3)
    z = 0
    c = 0
    with open(neg) as f:
        for line in f:
            if c >= n:
                sys.exit(0)
            r = line.split('\t')
            chr = r[0]
            start = int(r[1])
            stop = int(r[2])
            mid = (start + stop) / 2
            arr = getArray(chr,mid-250,mid+250,genLand)
            tot = sum(arr)
            if tot == 0:
                z = z + 1
                if z > n_zero:
                    continue
                else:
                    counts.write(','.join(map(str,arr)) + '\t0\n')
                    c = c + 1
                    continue
            else:
                counts.write(','.join(map(str,arr)) + '\t0\n')
                c = c + 1
