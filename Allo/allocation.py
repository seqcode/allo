#Lexi Morrissey, Mahony Lab @ Pennsylvania State University
#Last updated 09.05.2022
#Contains methods for read allocation procedure of Allo.

import predictPeak
import math
import random
import os
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.models import load_model, Model
import pickle
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import time 
import numpy as np
import sys


    
#Used to get counts in regions around multimapped reads
def getArray(read, winSize, genLand):
    array = []
    for k in range (int(read[3])-math.floor(winSize/2),int(read[3])+math.floor(winSize/2)):
        key = read[2] + ";" + str(k)
        #Seeing if current pos in genetic landscape
        if key in genLand:
            array.append(genLand[key])
        else:
            array.append(0)
    return array


#Assign reads (straight to dictionary for uniq and actual assign for multi-mapped)
def readAssign(rBlock, samOut, seed, winSize, genLand, model, cnn_scores, noread_score, rc, rmz):
    random.seed(seed)      #To make results reproducible
    
    ##Uniquely mapped reads##
    if len(rBlock) == 1:
        #Adding to file
        samOut.write('\t'.join(rBlock[0]))
        #Adding to genetic landscape
        key = rBlock[0][2] + ";" + str(int(rBlock[0][3]))
        if key in genLand:
            genLand[key] = genLand[key] + 1
        else:
            genLand[key] = 1
        return genLand

    ##Multi-mapped reads##
    ###CNN###
    scores = []
    basescores = []
    for i in rBlock:
        #Find closest 500 window, use that score instead if it's already been assigned, saves time
        pos = i[2]+str(round(int(i[3])/500)*500)
        if pos in cnn_scores:
            scores.append(cnn_scores[pos])
        else:
            countArray = getArray(i, winSize, genLand)
            s = sum(countArray)
            basescores.append(0)
            #Allocation options
            if rc == 1:
                if s == 0:
                    scores.append(1)
                else:
                    scores.append(s)
                continue
            if rc == 2:
                scores.append(1)
                continue
            #Use no read score if zero region
            if s == 0:
                scores.append(noread_score)
                cnn_scores[pos] = noread_score
            elif s <= 10:
                scores.append((noread_score+s)/2)
                cnn_scores[pos] = (noread_score+s)/2
            else:
                s = (predictPeak.predictNN(countArray, winSize, model)+s)/2
                scores.append(s)
                cnn_scores[pos] = s
    
    #Removing reads that mapped to all zero regions
    if sum(basescores) == 0 and rmz == 1:
        return genLand
        
    ##Choosing which read to keep based on probabilities from nearby read counts
    #Getting minimum value in list
    minval = np.amin([scores[i] for i in list(np.nonzero(scores)[0])])
    #Dividing all by minimum value so that nothing is less than 1
    counts = list(np.ceil(np.true_divide(scores,minval)).astype(int))
    
    #Procedure to allocate reads
    x = 0
    randN = random.randint(0,sum(counts))       #Used to randomly allocate reads
    for i in range(0, len(counts)):
        if randN in range(x,x+counts[i]):
            #Adding to file
            rBlock[i][-1] = rBlock[i][-1].strip()
            rBlock[i].append("ZA:Z:" + str(len(rBlock)) + "\n")
            samOut.write('\t'.join(rBlock[i]))
            '''
            key = rBlock[i][2] + ";" + str(int(rBlock[i][3]))
            if key in genLand:
              genLand[key] = genLand[key] + 1
            else:
              genLand[key] = 1
            ''' 
            return genLand
        x = x + counts[i]
    return genLand


#Getting uniquely mapped reads into temp files and putting them in the dictionary
def parseUniq(tempFile, seed, winSize, cnn_scores, AS, rc, keep):
    model = None
    noread_score = None
    rmz = None
    genLandCur = {}
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
                    continue
                if "N" in r[9]:
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
                        genLandCur = readAssign(rBlock, UM, seed, winSize, genLandCur, model, cnn_scores, noread_score, rc, rmz)
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
                        genLandCur = readAssign(rBlock, UM, seed, winSize, genLandCur, model, cnn_scores, noread_score, rc, rmz)
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
            genLandCur = readAssign(rBlock, UM, seed, winSize, genLandCur, model, cnn_scores, noread_score, rc, rmz)
    #Put multimapped reads into a file
        if len(rBlock) > 1:
            for i in range(0,len(rBlock)):
                MM.write('\t'.join(rBlock[i]))

    UM.close()
    MM.close()
    if b == 0:
        B.close()
    os.remove(tempFile) #Removing old temp
    return genLandCur


def parseMulti(tempFile, seed, winSize, genLand, modelName, cnn_scores, rc, keep, rmz, maxa):

    #Getting trained CNN
    try:
        json_file = open(modelName+'.json', 'r')
        loaded_model_json = json_file.read()
        json_file.close()
        model = tf.keras.models.model_from_json(loaded_model_json)
        model.load_weights(modelName+'.h5')
    except:
        print("Could not load Tensorflow model :( Allo was written with Tensorflow version 2.4.1")
        sys.exit(0)

    #Used for 0 count regions, dont keep predicting on them, waste of time
    countArray = []
    for i in range(0,winSize):
        countArray.append(0)
    noread_score = predictPeak.predictNN(countArray, winSize, model)

    #Exception that causes errors
    if os.stat(tempFile+"MM").st_size == 0:
        return
    
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
                    genLand = readAssign(rBlock, AL, seed, winSize, genLand, model, cnn_scores, noread_score, rc, rmz)
                    #Creating a new read block for the next read
                    rBlock = [] 
                    rBlock.append(r)

    #For last read
    if maxa is None or len(rBlock) <= maxa:
        genLand = readAssign(rBlock, AL, seed, winSize, genLand, model, cnn_scores, noread_score, rc, rmz)

    os.remove(tempFile) #Removing old temp
    AL.close()
    

#Getting uniquely mapped reads into temp files and putting them in the dictionary
def parseUniqPE(tempFile, seed, winSize, cnn_scores, AS, rc, keep, r2):
    model = None
    noread_score = None
    rmz = None
    if "border" in tempFile:
        b = 1
    else:
        b = 0
    genLandCur = {}
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
                            rBlock.pop()
                            rBlock2.pop()
                            continue
                    elif rBlock[-1][2] == "*" and rBlock2[-1][2] == "*":
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
                    UM.write('\t'.join(rBlock[-1]))
                    UM.write('\t'.join(rBlock2[-1]))
                if len(rBlock) == 1:
                    genLandCur = readAssignPE(rBlock, rBlock2, UM, seed, winSize, genLandCur, model, cnn_scores, noread_score, rc, rmz)
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
            genLandCur = readAssignPE(rBlock, rBlock2, UM, seed, winSize, genLandCur, model, cnn_scores, noread_score, rc, rmz)
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
    return genLandCur


#Assign reads (straight to dictionary for uniq and actual assign for multi-mapped)
def readAssignPE(rBlock, rBlock2, samOut, seed, winSize, genLand, model, cnn_scores, noread_score, rc, rmz):
    random.seed(seed)      #To make results reproducible
    
    ##Uniquely mapped reads##
    if len(rBlock) == 1:
        #Adding to file
        samOut.write('\t'.join(rBlock[0]))
        samOut.write('\t'.join(rBlock2[0]))
        #Adding to genetic landscape
        key = rBlock[0][2] + ";" + str(int(rBlock[0][3]))
        if key in genLand:
            genLand[key] = genLand[key] + 1
        else:
            genLand[key] = 1
        return genLand
        

    
    ##Multi-mapped reads##
    ###CNN###
    scores = []
    basescores = []
    for i in rBlock:
        #Find closest 500 window, use that score instead if it's already been assigned, saves time
        pos = i[2]+str(round(int(i[3])/500)*500)
        if pos in cnn_scores:
            scores.append(cnn_scores[pos])
        else:
            countArray = getArray(i, winSize, genLand)
            s = sum(countArray)
            basescores.append(s)
            #Allocation options
            if rc == 1:
                if s == 0:
                    scores.append(1)
                else:
                    scores.append(s)
                continue
            if rc == 2:
                scores.append(1)
                continue
            #Use no read score if zero region
            if s == 0:
                scores.append(noread_score)
                cnn_scores[pos] = noread_score
            elif s <= 10:
                scores.append((noread_score+s)/2)
                cnn_scores[pos] = (noread_score+s)/2
            else:
                s = (predictPeak.predictNN(countArray, winSize, model)+s)/2
                scores.append(s)
                cnn_scores[pos] = s
                
    #Removing reads that mapped to all zero regions
    if sum(basescores) == 0 and rmz == 1:
        return genLand
        
    ##Choosing which read to keep based on probabilities from nearby read counts
    #Getting minimum value in list
    minval = np.amin([scores[i] for i in list(np.nonzero(scores)[0])])
    #Dividing all by minimum value so that nothing is less than 1
    counts = list(np.ceil(np.true_divide(scores,minval)).astype(int))
    
    #Procedure to allocate reads
    x = 0
    randN = random.randint(0,sum(counts))       #Used to randomly allocate reads
    for i in range(0, len(counts)):
        if randN in range(x,x+counts[i]):
            #Adding to file
            rBlock[i][-1] = rBlock[i][-1].strip()
            rBlock2[i][-1] = rBlock2[i][-1].strip()
            rBlock[i].append("ZA:Z:" + str(len(rBlock)) + "\n")
            rBlock2[i].append("ZA:Z:" + str(len(rBlock)) + "\n")
            samOut.write('\t'.join(rBlock[i]))
            samOut.write('\t'.join(rBlock2[i]))
            return genLand
        x = x + counts[i]
    return genLand



def parseMultiPE(tempFile, seed, winSize, genLand, modelName, cnn_scores, rc, keep, rmz, maxa):

    #Getting trained CNN
    json_file = open(modelName+'.json', 'r')
    loaded_model_json = json_file.read()
    json_file.close()
    model = tf.keras.models.model_from_json(loaded_model_json)
    model.load_weights(modelName+'.h5')

    #Used for 0 count regions, dont keep predicting on them, waste of time
    countArray = []
    for i in range(0,winSize):
        countArray.append(0)
    noread_score = predictPeak.predictNN(countArray, winSize, model)

    #Exception that causes errors
    if os.stat(tempFile+"MM").st_size == 0:
        return
    
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
                genLand = readAssignPE(rBlock, rBlock2, AL, seed, winSize, genLand, model, cnn_scores, noread_score, rc, rmz)
                rBlock = []
                rBlock2 = []
                rBlock.append(r)
                curRead = r[0]
                numR = numR + 1
                numL = 1

    #For last read
    if maxa is None or len(rBlock) <= maxa:
        genLand = readAssignPE(rBlock, rBlock2, AL, seed, winSize, genLand, model, cnn_scores, noread_score, rc, rmz)

    os.remove(tempFile) #Removing old temp
    AL.close()
    
