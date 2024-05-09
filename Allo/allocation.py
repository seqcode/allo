#!/usr/bin/env python
#Lexi Morrissey, Mahony Lab @ Pennsylvania State University
#Contains methods for read allocation procedure of Allo.

from Allo import predictPeak
import math
import random
import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.models import load_model, Model
import pickle
import time 
import numpy as np
import sys
import multiprocessing
import re
import absl.logging
absl.logging.set_verbosity(absl.logging.ERROR)
import contextlib, io

#Add reads to UMR dictionary
def addToDict(tempFile, genLand, spliceD, seq):
    count = 0
    with open(tempFile) as f:
        for line in f:
            l = line.strip().split('\t')
            if l[2] == "*":
                continue
            elif "N" in l[9]:
                continue
            count = count + 1
            if seq == 0 and not (count % 2) == 0:
                continue
            l[3] = int(l[3])
            #Add to dictionary
            key = l[2] + ";" + str(l[3])
            if key in genLand:
                genLand[key] = genLand[key] + 1
            else:
                genLand[key] = 1
            if not spliceD == None:
                chars = re.findall('[A-Za-z]+', l[5])
                gaps = [i for i, x in enumerate(chars) if x == "N"]
                if len(gaps) > 0:
                    num = re.findall('\d+\.\d+|\d+', l[5])
                    curLoc = int(l[3])
                    for i in range(0,len(chars)):
                        if chars[i]=="M" or chars[i]=="D" or chars[i]=="=" or chars[i]=="X":
                            curLoc = curLoc + int(num[i])
                        elif chars[i]=="N":
                            if l[2]+";"+str(curLoc) in spliceD:
                                if spliceD[l[2]+";"+str(curLoc)] > curLoc + int(num[i]):
                                    pos = int(curLoc)
                                    size = int(num[i])
                                    spliceD[l[2]+";"+str(pos)] = size
                                    spliceD[l[2]+";"+str(pos+size)] = size
                            else:
                                pos = int(curLoc)
                                size = int(num[i])
                                spliceD[l[2]+";"+str(pos)] = size
                                #Adding splice in reverse direction
                                spliceD[l[2]+";"+str(pos+size)] = -1*(size)
                            curLoc = curLoc + int(num[i])

#Used to get counts in regions around multimapped reads
def getArray(read, winSize, genLand, spliceD):
    #Spliced array using junction info
    if spliceD:
        #Parsing cigar string
        chars = re.findall('[A-Za-z]+', read[5])
        num = re.findall('\d+\.\d+|\d+', read[5])
        gap_loc = {}
        start = int(read[3])
        chr = read[2]
        r_end = start #Get information about end of read to use for splice variants
        #print("cigar: " + read[5], flush=True)
        for i in range(0,len(chars)):
            if chars[i]=="M" or chars[i]=="D" or chars[i]=="=" or chars[i]=="X":
                r_end = r_end + int(num[i]) - 1
            elif chars[i]=="N":
                gap_loc[r_end+1]=int(num[i])
                r_end = r_end + int(num[i])
        array = []
        k = start
        l = 0
        #Upstream counts
        while l <= math.floor(winSize/2):
            #print("up: " + str(k), flush=True)
            key = chr + ";" + str(k)
            '''if key in spliceD and spliceD[key] < 0:
                #print("splice: " + str(splice[key]), flush=True)
                k = k + spliceD[key]
                key = chr + ";" + str(k)'''
            if key in genLand:
                array.insert(0,genLand[key])
            else:
                array.insert(0,0)
            k -= 1
            l += 1
        start_pos = k
        #Downstream counts 
        k = start
        l = 0
        while l <= math.floor(winSize/2):
            #Taking splice information from read if present
            while k < r_end:
                key = chr + ";" + str(k)
                if k in gap_loc:
                    k = k + gap_loc[k]
                    continue
                elif key in genLand:
                    array.append(genLand[key])
                else:
                    array.append(0)
                k += 1
                l += 1
            #Else use previously found splice info
            key = chr + ";" + str(k)
            '''if key in spliceD and spliceD[key] > 0:
                k = k + spliceD[key]
                key = chr + ";" + str(k)'''
            if key in genLand:
                array.append(genLand[key])
            else:
                array.append(0)
            k += 1
            l += 1
    #Non-spliced array
    else:
        array = []
        pos = round(int(read[3])/100)*100
        for k in range (pos-math.floor(winSize/2),int(pos)+math.floor(winSize/2)):
            key = read[2] + ";" + str(k)
            #Seeing if current pos in genetic landscape
            if key in genLand:
                array.append(genLand[key])
            else:
                array.append(0)
    return array


#Assign reads (straight to dictionary for uniq and actual assign for multi-mapped)
def readAssign(rBlock, samOut, winSize, genLand, model, cnn_scores, rc, rmz, spliceD):
    random.seed(7)      #To make results reproducible

    ##Multi-mapped reads##
    ###CNN###
    scores_rc = []
    scores_nn = []
    allZ = True #seeing if all zero regions
    for i in rBlock:
        #Find closest 100 window, use that score instead if it's already been assigned, saves time
        if spliceD:
            pos = i[2]+";"+str(i[3])
        else:
            pos = i[2]+";"+str(round(int(i[3])/100)*100)
        if pos in cnn_scores:
            scores_nn.append(cnn_scores[pos])
            allZ = False
        else:
            countArray = getArray(i, winSize, genLand, spliceD)
            s = sum(countArray)
            if s > 0:
                allZ = False
            #Allocation options
            if rc == 1:
                if s == 0:
                    scores_rc.append(1)
                else:
                    scores_rc.append(s+1)
                continue
            if rc == 2:
                scores_rc.append(1)
                continue
            #Use no read score if zero region
            if s == 0:
                scores_nn.append(0.0012*(s+1))
            else:
                nn = predictPeak.predictNN(countArray, winSize, model)
                scores_nn.append(nn*(s+1))
                cnn_scores[pos] = (nn*(s+1))
    
    #Removing reads that mapped to all zero regions
    if allZ and rmz == 1:
        return
        
    ##Choosing which read to keep based on probabilities from nearby read counts
    #List clean up
    if len(scores_nn)>0:
        try:
            percs = scores_nn / np.sum(scores_nn)
        except:
            print("Issue with array")
            print(scores_nn)
            print(scores_rc, flush=True)
            sys.exit(0)
    else:
        percs = scores_rc / np.sum(scores_rc)

    #Picking read to allocate based on percentages
    np.random.seed(7)
    choice = np.random.choice(range(0,len(rBlock)), p=percs, replace=True, size=1)[0]
    
    #Adding to file
    rBlock[choice][-1] = rBlock[choice][-1].strip()
    if not allZ:
        rBlock[choice].append("ZA:Z:" + str(len(rBlock)) + "\n")
    else:
        rBlock[choice].append("ZZ:Z:" + str(len(rBlock)) + "\n")
    samOut.write('\t'.join(rBlock[choice]))
    return


#Getting uniquely mapped reads into temp files and putting them in the dictionary
def parseUniq(tempFile, winSize, cnn_scores, AS, rc, keep):
    #Variables for stats
    cu = 0  #UMRs
    cf = 0  #Filtered
    
    model = None
    rmz = None
    rBlock = []
    if "border" in tempFile:
        b = 1
    else:
        b = 0
    #File with uniquely mapped reads
    UM = open(tempFile + "UM","w+")
    #File with unallocated multimapped reads
    MM = open(tempFile + "MM","w+")
    if b == 0:
        B = open(tempFile + "B","w+")
    
    numR = 0 #Need to get the reads bordering on the thread cuts
    with open(tempFile) as f:
        for line in f:
            if not line.strip():
                break
            
            #First file will have sam header in it
            if line.startswith("@"):
                continue
                
            #Splitting columns
            r = line.split('\t')
            
            #Keep unmapped reads if user chooses
            if keep == 0:
                if r[2] == "*":
                    cf += 1
                    continue
                if "N" in r[9]:
                    cf += 1
                    continue
              

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
                    if numR == 1 and b == 0:
                        for i in range(0,len(rBlock)):
                            B.write('\t'.join(rBlock[i]))
                        rBlock = []
                        rBlock.append(r)
                        numR = numR + 1
                        continue
                    numR = numR + 1
                    #Put uniquely mapped reads into a file
                    if len(rBlock) == 1:
                        UM.write('\t'.join(rBlock[0]))
                        cu += 1
                    #Put multimapped reads into a file
                    if len(rBlock) > 1:
                        for i in range(0,len(rBlock)):
                            MM.write('\t'.join(rBlock[i]))
                            
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
                    if numR == 1 and b == 0:
                        for i in range(0,len(rBlock)):
                            B.write('\t'.join(rBlock[i]))
                        rBlock = []
                        rBlock.append(r)
                        numR = numR + 1
                        continue
                    numR = numR + 1
                    #Put uniquely mapped reads into a file and into dictionary
                    if len(rBlock) == 1:
                        UM.write('\t'.join(rBlock[0]))
                        cu += 1
                    #Put multimapped reads into a file only, no dictionary assign
                    if len(rBlock) > 1:
                        for i in range(0,len(rBlock)):
                            MM.write('\t'.join(rBlock[i]))
                            
                    #Creating a new read block for the next read
                    rBlock = []
                    rBlock.append(r)
    
    #Running function one more time for last read
    if b == 0:
        for i in range(0,len(rBlock)):
            B.write('\t'.join(rBlock[i]))
    if b == 1:
        if len(rBlock) == 1:
            UM.write('\t'.join(rBlock[0]))
            cu += 1
    #Put multimapped reads into a file
        if len(rBlock) > 1:
            for i in range(0,len(rBlock)):
                MM.write('\t'.join(rBlock[i]))
                

    UM.close()
    MM.close()
    if b == 0:
        B.close()
    os.remove(tempFile) #Removing old temp
    return [cu, cf]


def parseMulti(tempFile, winSize, genLand, modelName, cnn_scores, rc, keep, rmz, maxa, spliceD):   
    numLoc = [0,0] #Keep info on average number of places read maps to
    #Getting trained CNN
    if rc == 0:
        try:
            json_file = open(modelName+'.json', 'r')
            loaded_model_json = json_file.read()
            json_file.close()
            model = tf.keras.models.model_from_json(loaded_model_json)
            model.load_weights(modelName+'.h5')
            model = LiteModel.from_keras_model(model)
        except:
            print("Could not load Tensorflow model :( Allo was written with Tensorflow version 2.11")
            sys.exit(0)
    else:
        model = None
    #Exception that causes errors
    if os.stat(tempFile+"MM").st_size == 0:
        return numLoc
    
    #File with MM allocations
    AL = open(tempFile + "AL","w+")
    tempFile = tempFile + "MM"
    rBlock = []
    with open(tempFile) as f:
        for line in f:
            if not line.strip():
                break
            
            #Splitting columns
            r = line.split('\t')
            
            #Keep unmapped reads if user chooses
            if keep == 0:
                if r[2] == "*":
                    continue
                if "N" in r[9]:
                    continue
                
            #Appending if the block is empty and going to next line, happens first line
            if len(rBlock) == 0:
                rBlock.append(r)
                continue
            
            #Allocating multimapped reads
            else:
                leadr = rBlock[0]

                #Adding to block
                if r[0] == leadr[0]:
                    rBlock.append(r)
                    continue
                
                #New block, have to allocate reads
                if not r[0] == leadr[0]:
                    if maxa is not None and len(rBlock) > maxa:
                        rBlock = []
                        rBlock.append(r)
                        continue
                    readAssign(rBlock, AL, winSize, genLand, model, cnn_scores, rc, rmz, spliceD)
                    #Getting average number of locations mapped to
                    numLoc[0] = (numLoc[0]*numLoc[1] + len(rBlock)) / (numLoc[1]+1)
                    numLoc[1] = numLoc[1] + 1
                    #Creating a new read block for the next read
                    rBlock = [] 
                    rBlock.append(r)

    #For last read
    if maxa is None or len(rBlock) <= maxa:
        readAssign(rBlock, AL, winSize, genLand, model, cnn_scores, rc, rmz, spliceD)
        numLoc[0] = (numLoc[0]*numLoc[1] + len(rBlock)) / (numLoc[1]+1)
        numLoc[1] = numLoc[1] + 1

    os.remove(tempFile) #Removing old temp
    AL.close()
    return numLoc
    

#Getting uniquely mapped reads into temp files and putting them in the dictionary
def parseUniqPE(tempFile, winSize, cnn_scores, AS, rc, keep, r2):
    #Keeping information on read counts
    cu = 0
    cf = 0
    
    model = None
    rmz = None
    if "border" in tempFile:
        b = 1
    else:
        b = 0
    rBlock = []
    rBlock2 = []
    #File with uniquely mapped reads
    UM = open(tempFile + "UM","w+")
    #File with unallocated multimapped reads
    MM = open(tempFile + "MM","w+")
    #File with border reads due to thread split
    if b == 0:
        B = open(tempFile + "B","w+")
    
    numR = 0 #Need to get the reads bordering on the thread cuts
    curRead = '' #Keep track of current read name
    pairs = {} #Dictionary used to match pairs
    with open(tempFile) as f:
        for line in f:
            #Exception
            if not line.strip():
                break
            
            #First file will have sam header in it
            if line.startswith("@"):
                continue
        
            #Splitting columns
            r = line.split('\t')
            #For first read which is a border read
            if numR == 0:
                curRead = r[0]
                if b == 0:
                    B.write('\t'.join(r))
                else:
                    #Getting read info
                    r1chr = str(r[2])
                    r1pos = str(r[3])
                    key = r1chr + ";" + r1pos
                    #Getting mate info
                    r2chr = str(r[6])
                    if r2chr == "=":
                        r2chr = str(r[2])
                    if r2chr == "*":
                        r2chr = "*"
                        r2pos = "0"
                    else:
                        r2pos = str(r[7])
                    key_mate = r2chr + ";" + r2pos
                    #Restarting dictionary
                    pairs[key_mate + ":" + key] = r
                numR = 1
                continue
            
            #Just adding to arrays if it's the same read
            if curRead == r[0]:
                if numR == 1 and b == 0:
                    B.write('\t'.join(r)) 
                    continue
                else:
                    #Adding to dictionary of pairs or rBlock if a proper mate is found
                    #Gettin read info
                    r1chr = str(r[2])
                    r1pos = str(r[3])
                    key = r1chr + ";" + r1pos
                    #Getting mate info
                    r2chr = str(r[6])
                    if r2chr == "=":
                        r2chr = str(r[2])
                    if r2chr == "*":
                        r2chr = "*"
                        r2pos = "0"
                    else:
                        r2pos = str(r[7])
                    key_mate = r2chr + ";" + r2pos
                    #Getting pair key
                    currPair = key + ":" + key_mate
    
                    if currPair in pairs:
                        #Using R2 to make allocations if user specifies
                        if r2 == 0:
                            if bin(int(pairs[currPair][1]))[-7] == 1:    #Getting r1 based on sam flags
                                rBlock.append(pairs[currPair])
                                rBlock2.append(r)
                            else:
                                rBlock.append(r)
                                rBlock2.append(pairs[currPair])
                            del pairs[currPair]
                        else:
                            if bin(int(pairs[currPair][1]))[-7] == 1:    #Getting r1 based on sam flags
                                rBlock2.append(pairs[currPair])
                                rBlock.append(r)
                            else:
                                rBlock2.append(r)
                                rBlock.append(pairs[currPair])
                            del pairs[currPair]
                        
                    else:
                        pairs[key_mate + ":" + key] = r
                        continue
                    
                    #Keeping unmapped reads if users choose
                    if keep == 0:
                        if rBlock2[-1][2] == "*" or "N" in rBlock2[-1][9] or rBlock[-1][2] == "*" or "N" in rBlock[-1][9]:
                            cf += 2
                            rBlock.pop()
                            rBlock2.pop()
                            continue
                    elif rBlock[-1][2] == "*" and rBlock2[-1][2] == "*":
                        cu += 1
                        UM.write('\t'.join(rBlock[-1]))
                        UM.write('\t'.join(rBlock2[-1]))
                        rBlock.pop()
                        rBlock2.pop()
                        continue
                    #See if these new reads have an equal score to whats already in rBlock
                    if not AS == 0 and len(rBlock)>1:
                        if len(rBlock[-1])>AS and len(rBlock2[-1])>AS and len(rBlock[-2])>AS and \
                        len(rBlock2[-2])>AS and rBlock[-1][AS].startswith("AS") and \
                        rBlock2[-1][AS].startswith("AS") and rBlock[-2][AS].startswith("AS") and \
                        rBlock2[-2][AS].startswith("AS"):
                            AS1 = int(rBlock[-1][AS].split(':')[2])
                            AS11 = int(rBlock2[-1][AS].split(':')[2])
                            AS2 = int(rBlock[-2][AS].split(':')[2])
                            AS22 = int(rBlock2[-2][AS].split(':')[2])
                        else:
                            #Fix issues with alignment score column
                            #R1-1
                            if rBlock[-1][2] == "*":
                                AS1 = 0
                            else:
                                ASt = 0
                                for t in range(0,len(rBlock[-1])):
                                    if rBlock[-1][t].startswith("AS:"):
                                        ASt = t
                                AS1 = int(rBlock[-1][ASt].split(':')[2])
                            #R2-1
                            if rBlock2[-1][2] == "*":
                                AS11 = 0
                            else:
                                ASt = 0
                                for t in range(0,len(rBlock2[-1])):
                                    if rBlock2[-1][t].startswith("AS:"):
                                        ASt = t
                                AS11 = int(rBlock2[-1][ASt].split(':')[2])
                            #R1-2  
                            if rBlock[-2][2] == "*":
                                AS2 = 0
                            else:
                                ASt = 0
                                for t in range(0,len(rBlock[-2])):
                                    if rBlock[-2][t].startswith("AS:"):
                                        ASt = t
                                AS2 = int(rBlock[-2][ASt].split(':')[2])
                            #R2-2
                            if rBlock2[-2][2] == "*":
                                AS22 = 0
                            else:
                                ASt = 0
                                for t in range(0,len(rBlock2[-2])):
                                    if rBlock2[-2][t].startswith("AS:"):
                                        ASt = t
                                AS22 = int(rBlock2[-2][ASt].split(':')[2])
                        if AS1 + AS11 == AS2 + AS22:
                            continue
                        if AS1 + AS11 > AS2 + AS22:
                            rBlock = [rBlock[-1]]
                            rBlock2 = [rBlock2[-1]]
                            continue
                        if AS1 + AS11 < AS2 + AS22:
                            rBlock.pop()
                            rBlock2.pop()
                            continue
                            
            #If different read, assign it to a file
            else:
                if len(rBlock) == 1 and rBlock[-1][2] == "*" and rBlock2[-1][2] == "*":
                    cu += 1
                    UM.write('\t'.join(rBlock[-1]))
                    UM.write('\t'.join(rBlock2[-1]))
                if len(rBlock) == 1:
                    cu += 1
                    UM.write('\t'.join(rBlock[-1]))
                    UM.write('\t'.join(rBlock2[-1]))
                else:
                    for i in range(0,len(rBlock)):
                        MM.write('\t'.join(rBlock[i]))
                        MM.write('\t'.join(rBlock2[i]))
                rBlock = []
                rBlock2 = []
                pairs = {}
                #Getting read info
                r1chr = str(r[2])
                r1pos = str(r[3])
                key = r1chr + ";" + r1pos
                #Getting mate info
                r2chr = str(r[6])
                if r2chr == "=":
                    r2chr = str(r[2])
                if r2chr == "*":
                    r2chr = "*"
                    r2pos = "0"
                else:
                    r2pos = str(r[7])
                key_mate = r2chr + ";" + r2pos
                #Restarting dictionary
                pairs[key_mate + ":" + key] = r
                curRead = r[0]
                numR = numR + 1

    #Running function one more time for last read
    if b == 0:
        unpair = pairs.values()
        for i in range(0,len(rBlock)):
            B.write('\t'.join(rBlock[i]))
            B.write('\t'.join(rBlock2[i]))
        for k in pairs.values():
            B.write('\t'.join(k))
    if b == 1:
        if len(rBlock) == 1:
            UM.write('\t'.join(rBlock[-1]))
            UM.write('\t'.join(rBlock2[-1]))
            cu += 1
        #Put multimapped reads into a file
        if len(rBlock) > 1:
            for i in range(0,len(rBlock)):
                MM.write('\t'.join(rBlock[i]))
                MM.write('\t'.join(rBlock2[i]))
    
    UM.close()
    MM.close()
    if b == 0:
        B.close()
    os.remove(tempFile) #Removing old temp
    return [cu, cf]


#Assign reads (straight to dictionary for uniq and actual assign for multi-mapped)
def readAssignPE(rBlock, rBlock2, samOut, winSize, genLand, model, cnn_scores, rc, rmz, spliceD):
    random.seed(7)      #To make results reproducible
    ##Multi-mapped reads##
    ###CNN###
    scores_rc = []
    scores_nn = []
    allZ = True #seeing if all zero regions
    for i in rBlock:
        #Find closest 100 window, use that score instead if it's already been assigned, saves time
        if spliceD:
            pos = i[2]+";"+str(i[3])
        else:
            pos = i[2]+";"+str(round(int(i[3])/100)*100)
        if pos in cnn_scores:
            scores_nn.append(cnn_scores[pos])
            allZ = False
        else:
            countArray = getArray(i, winSize, genLand, spliceD)
            s = sum(countArray)
            if s > 0:
                allZ = False
            #Allocation options
            if rc == 1:
                if s == 0:
                    scores_rc.append(1)
                else:
                    scores_rc.append(s+1)
                continue
            if rc == 2:
                scores_rc.append(1)
                continue
            #Use no read score if zero region
            if s == 0:
                scores_nn.append(0.0012*(s+1))
            else:
                nn = predictPeak.predictNN(countArray, winSize, model)
                scores_nn.append(nn*(s+1))
                cnn_scores[pos] = (nn*(s+1))
    
    #Removing reads that mapped to all zero regions
    if allZ and rmz == 1:
        return
        
    ##Choosing which read to keep based on probabilities from nearby read counts
    #List clean up
    if len(scores_nn)>0:
        try:
            percs = scores_nn / np.sum(scores_nn)
        except:
            print("Issue with array")
            print(scores_nn)
            print(scores_rc, flush=True)
            sys.exit(0)
    else:
        percs = scores_rc / np.sum(scores_rc)

    #Picking read to allocate based on percentages
    np.random.seed(7)
    choice = np.random.choice(range(0,len(rBlock)), p=percs, replace=True, size=1)[0]
    
    #Adding to file
    rBlock[choice][-1] = rBlock[choice][-1].strip()
    rBlock2[choice][-1] = rBlock2[choice][-1].strip()
    if not allZ:
        rBlock[choice].append("ZA:Z:" + str(len(rBlock)) + "\n")
        rBlock2[choice].append("ZA:Z:" + str(len(rBlock)) + "\n")
    else:
        rBlock[choice].append("ZZ:Z:" + str(len(rBlock)) + "\n")
        rBlock2[choice].append("ZZ:Z:" + str(len(rBlock)) + "\n")
    samOut.write('\t'.join(rBlock[choice]))
    samOut.write('\t'.join(rBlock2[choice]))
    return


def parseMultiPE(tempFile, winSize, genLand, modelName, cnn_scores, rc, keep, rmz, maxa, spliceD):
    numLoc = [0,0] #Retain info on number of mapping sites
    #Getting trained CNN and making sure there is a compatible tensorflow installed
    if rc == 0:
        try:
            json_file = open(modelName+'.json', 'r')
            loaded_model_json = json_file.read()
            json_file.close()
            model = tf.keras.models.model_from_json(loaded_model_json)
            model.load_weights(modelName+'.h5')
            model = LiteModel.from_keras_model(model)
        except:
            print("Could not load Tensorflow model :( Allo was written with Tensorflow version 2.11")
            sys.exit(0)
    else:
        model = None

    #Exception that causes errors
    if os.stat(tempFile+"MM").st_size == 0:
        return numLoc
    
    #File with MM allocations
    AL = open(tempFile + "AL","w+")
    tempFile = tempFile + "MM"
    rBlock = []
    rBlock2 = []
    numR = 0 #Need to get the reads bordering on the thread cuts
    numL = 0    #Keep track of what line in a block
    curRead = '' #Keep track of current read name
    with open(tempFile) as f:
        for line in f:
            #Exception
            if not line.strip():
                break
                
            #Splitting columns
            r = line.split('\t')
            
            #For first read 
            if numR == 0:
                curRead = r[0]
                rBlock.append(r)
                numR = 1
                numL = numL + 1
                continue
            
           #Just adding to arrays if it's the same read
            if curRead == r[0]:
                if numR == 1 and numL == 1:
                    rBlock2.append(r)
                    numL = numL + 1
                    continue
                if numL%2 == 0:
                    rBlock.append(r)
                    numL = numL + 1  
                else:
                    #Exceptions that cause errors
                    if keep == 0:
                        if r[2] == "*" or "N" in r[9]:
                            rBlock.pop()
                            rBlock2.pop()
                    #Making sure pairs are correctly together
                    r2chr = rBlock[-1][6]
                    if r2chr == "=":
                        r2chr = rBlock[-1][2]
                    r2pos = rBlock[-1][7]
                    if r[2] == r2chr and r[3] == r2pos: 
                        rBlock2.append(r)
                        numL = numL + 1
                    else:
                        #Error if not sorted properly
                        print("Pairs not sorted properly. Issue with following entry.")
                        print(r[0])
                        sys.exit(0)
            #If different read, begin read assignment
            else:
                if maxa is not None and len(rBlock) > maxa:
                    continue
                readAssignPE(rBlock, rBlock2, AL, winSize, genLand, model, cnn_scores, rc, rmz, spliceD)
                numLoc[0] = (numLoc[0]*numLoc[1] + len(rBlock)) / (numLoc[1]+1)
                numLoc[1] = numLoc[1] + 1
                rBlock = []
                rBlock2 = []
                rBlock.append(r)
                curRead = r[0]
                numR = numR + 1
                numL = 1

    #For last read
    if maxa is None or len(rBlock) <= maxa:
        readAssignPE(rBlock, rBlock2, AL, winSize, genLand, model, cnn_scores, rc, rmz, spliceD)
        numLoc[0] = (numLoc[0]*numLoc[1] + len(rBlock)) / (numLoc[1]+1)
        numLoc[1] = numLoc[1] + 1

    os.remove(tempFile) #Removing old temp
    AL.close()
    
    return numLoc


#Class to speed up tensorflow prediction
#https://micwurm.medium.com/using-tensorflow-lite-to-speed-up-predictions-a3954886eb98
class LiteModel:
    @classmethod
    def from_file(cls, model_path):
        return LiteModel(tf.lite.Interpreter(model_path=model_path))
    
    @classmethod
    def from_keras_model(cls, kmodel):
        converter = tf.lite.TFLiteConverter.from_keras_model(kmodel)
        tflite_model = converter.convert()
        return LiteModel(tf.lite.Interpreter(model_content=tflite_model))
    
    def __init__(self, interpreter):
        self.interpreter = interpreter
        self.interpreter.allocate_tensors()
        input_det = self.interpreter.get_input_details()[0]
        output_det = self.interpreter.get_output_details()[0]
        self.input_index = input_det["index"]
        self.output_index = output_det["index"]
        self.input_shape = input_det["shape"]
        self.output_shape = output_det["shape"]
        self.input_dtype = input_det["dtype"]
        self.output_dtype = output_det["dtype"]
        
    def predict(self, inp):
        inp = inp.astype(self.input_dtype)
        count = inp.shape[0]
        out = np.zeros((count, self.output_shape[1]), dtype=self.output_dtype)
        for i in range(count):
            self.interpreter.set_tensor(self.input_index, inp[i:i+1])
            self.interpreter.invoke()
            out[i] = self.interpreter.get_tensor(self.output_index)[0]
        return out
    
    def predict_single(self, inp):
        inp = np.array([inp], dtype=self.input_dtype)
        self.interpreter.set_tensor(self.input_index, inp)
        self.interpreter.invoke()
        out = self.interpreter.get_tensor(self.output_index)
        return out[0]


