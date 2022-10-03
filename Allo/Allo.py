#Lexi Morrissey, Mahony Lab @ Pennsylvania State University
#Last updated 9.06.2022
#Main method for Allo. Splits the sam files up and sends them to allocation procedure via multiprocessing package.

import sys

#Help message
if "-h" in sys.argv or len(sys.argv)<3:
    print("\nusage: Allo [-h]")
    print("Developed by Alexis Morrissey, Mahony Laboratory @ The Pennsylvania State University")
    print("Allo is a python program used allocate multi-mapped reads in ChIP-seq data.")
    print("For more pre-processing info and basic usage, please visit https://github.com/seqcode/allo" + "\n\n")
    print("##REQUIRED ARGUMENTS##")
    print("-seq      {se/pe}      Single-end or paired-end sequencing mode" + "\n\n")
    print("##OPTIONAL ARGUMENTS##")
    print("-o      {string}      Output file name")
    print("-m      {mixed/narrow}      Use CNN trained on either a narrow peak dataset or a dataset with mixed peaks, narrow by default")
    print("-t      {int}      Number of threads, 1 by default")
    print("-max      {int}      Maximum value for number of locations a read can map to")
    print("--keep-unmap      Keep unmapped reads and reads that include N in their sequence")
    print("--remove-zeros      Disregard multi-mapped reads that map to regions with 0 uniquely mapped reads (random assignment)")
    print("--r2      Use read 2 for allocation procedure instead of read 1 (only for paired-end sequencing)")
    print("--readcount      CNN will not be used in allocation, only read counts")
    print("--random      Reads will be randomly assigned (similar to Bowtie and BWA on default)")
    print("\n")
    sys.exit(0)

#Imports
from joblib import Parallel, delayed
import multiprocessing
import subprocess
import os
import allocation
import math
from random import randint
import shutil
import glob

#Function for concatenating files
def cat(files, outname):
    with open(outname,'wb') as wfd:
        for f in files:
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd)

##Main Method##          
if __name__ == '__main__':

    pos_args = ['-o','--mixed','-s','-t','-seq','--readcount','--random','--keep-unmap','--remove-zeros','-max','--r2','--ignore']
    for arg in sys.argv[2:]:
        if arg[0] == "-" and arg not in pos_args:
            print("Argument not recognized: " + arg)
            sys.exit(0)
    
    #Input name
    try:
        samfile = sys.argv[1]
    except:
        print("SAM file not given")
    #Output name
    if "-o" in sys.argv:
        outfile = sys.argv[sys.argv.index("-o")+1]
    else:
        outfile = samfile[:-4] + ".allo.sam"
    #Model and window size
    if "--mixed" in sys.argv:
        m = os.path.dirname(sys.argv[0])+"/mixed"
        winSize = 500
    else:
        m = os.path.dirname(sys.argv[0])+"/narrow"
        winSize = 500
    #Seed
    if "-s" in sys.argv:
        seed = int(sys.argv[sys.argv.index("-s")+1])
    else:
        seed = 7
    #Number of threads
    if "-t" in sys.argv:
        thr = int(sys.argv[sys.argv.index("-t")+1])
    else:
        thr = 1  
    #Single or paired end
    if "-seq" in sys.argv:
        if sys.argv[sys.argv.index("-seq")+1] == "se":
            seq = 1
        elif sys.argv[sys.argv.index("-seq")+1] == "pe":
            seq = 0
        else:
            print("ERROR: Sequencing type not recognized. Please specify single-end (se) or paired-end (pe).")
            sys.exit(0)
    else:
        print("ERROR: Sequencing type invalid or not specified.")
    #Options for probabilities
    if "--readcount" in sys.argv:
        rc = 1
    elif "--random" in sys.argv:
        rc = 2
    else:
        rc = 0
    #Keep unmapped reads
    if "--keep-unmap" in sys.argv:
        keep = 1
    else:
        keep = 0
    #Remove reads that only map to zero regions
    if "--remove-zeros" in sys.argv:
        rmz = 1
    else:
        rmz = 0
    #User defined max number of alignments
    if "-max" in sys.argv:
        maxa = int(sys.argv[sys.argv.index("-max")+1])
    else:
        maxa = None 
    #Use read 2 instead of read 1 for allocation
    if "--r2" in sys.argv:
        r2 = 1
    else:
        r2 = 0
    #Ignore allo warnings
    if "--ignore" in sys.argv:
        ig = 1
    else:
        ig = 0
    

    #Need list of file names
    ids = str(randint(0, 10000))
    tempList = []
    for i in range (0,thr):
        tempList.append("allo.temp." + ids + "." + str(i).zfill(2))
        
    #Getting header information and checking to see if file was sorted
    collate = 0
    f = open(samfile)
    header = open('header'+ids, 'a')
    for line in f:
        if line.startswith('@'):
            header.write(line)
            if "collate" in line:
                collate = 1
        else:
            break
    f.close()
    header.close()
    if collate == 0 and ig == 0:
        print("Warning: PG tag not present in header. File does not appear to be sorted.")
        print("Please use Samtools collate before running Allo. To ignore this warning use argument --ignore")
        sys.exit(0)
    elif collate == 0 and ig == 1:
        print("Suppressed warning: PG tag not present in header. File does not appear to be sorted.")
    
    #Get number of lines
    count = 0
    file = open(samfile, 'rb')
    while True:
        buffer = file.read(8192*1024)
        if not buffer: break
        count += buffer.count('\n'.encode())
    file.close()
    
    
    #Splitting up file
    size = int(count/len(tempList))
    f = open(samfile)
    linec = 0
    i = 0
    filen = open(tempList[0], 'a')
    for line in f:
        if linec <= size:
            filen.write(line)
            linec = linec + 1
        else:
            filen.write(line)
            filen.close()
            linec = 0
            i = i + 1
            filen = open(tempList[i], 'a')
    filen.close()

        
    #Getting column where alignment score is stored in sam (different for every aligner). If none listed assumes bowtie1 -m 1 -k x
    with open(tempList[-1]) as f:
        for line in f:
            if not line.startswith("@") and not line[2] == "*":
                l = line
                break
    AS = 0
    l = l.split('\t')
    for i in range(0,len(l)):
        if l[i].startswith("AS:"):
            AS = i

    ##Single-end mode##
    if seq == 1:
        #PHASE I: Parsing unique reads and putting into dictionary
        cnn_scores = {} #Used to keep some previous scores for places in the genome by the cnn to save time
        output = Parallel(n_jobs=thr)(delayed(allocation.parseUniq)(i, seed, winSize, cnn_scores, AS, rc, keep) for i in tempList)
        #Updating genetic landscape with uniquely mapped reads
        genLand = {}
        for d in output:
            genLand.update(d)
        #Need to run allocation on border reads
        if len(tempList) > 1:
            l = []
            for i in tempList:
                if os.path.exists(str(i)+"B"):
                    l.append(i + "B")
            cat(l,tempList[0]+".borders")
        else:
            os.rename(tempList[0] + 'B', tempList[0]+".borders")
            
        borderDict = allocation.parseUniq(tempList[0]+".borders", seed, winSize, cnn_scores, AS, rc, keep)
        d.update(borderDict)
        #Cat with first MM file 
        cat([tempList[0]+'MM', tempList[0]+'.bordersMM'], tempList[0]+'.catborders')
        #Fix name
        os.rename(tempList[0]+'.catborders', tempList[0]+'MM')
        #Cat with first UM file
        cat([tempList[0]+'UM', tempList[0]+'.bordersUM'], tempList[0]+'.catborders')
        #Fix name
        os.rename(tempList[0]+".catborders", tempList[0]+"UM")
        #Remove old files
        for filename in glob.glob('./*borders*'):
            os.remove(filename) 
        
        
        #PHASE II: Parsing multi-mapped reads
        output = []
        output = Parallel(n_jobs=thr)(delayed(allocation.parseMulti)(i, seed, winSize, genLand, m, cnn_scores, rc, keep, rmz, maxa) for i in tempList)
        
        
        #Final parsing and file clean up
        #Deleting temporary files and join allocated files
        l = ["header"+ids]
        for i in tempList:
            if os.path.exists(str(i)+"AL"):
                l.append(i + "AL")
            l.append(i + "UM")
        cat(l, outfile)
        
        #Remove all temporary files
        for filename in glob.glob('allo.temp.*'):
            os.remove(filename)
        
        os.remove("header"+ids)

        
    
    ##Paired-end mode##
    if seq == 0:
    
        #PHASE I: Parsing unique reads and putting into dictionary
        cnn_scores = {} #Used to keep some previous scores for places in the genome by the cnn to save time
        output = Parallel(n_jobs=thr)(delayed(allocation.parseUniqPE)(i, seed, winSize, cnn_scores, AS, rc, keep, r2) for i in tempList)
        #Updating genetic landscape with uniquely mapped reads
        genLand = {}
        for d in output:
            genLand.update(d)
        #Need to run allocation on border reads
        if len(tempList) > 1:
            l = []
            for i in tempList:
                if os.path.exists(str(i)+"B"):
                    l.append(i + "B")
            cat(l, tempList[0]+".borders")
        else:
            os.rename(tempList[0]+'B', tempList[0]+".borders")   
        borderDict = allocation.parseUniqPE(tempList[0]+".borders", seed, winSize, cnn_scores, AS, rc, keep, r2)
        d.update(borderDict)
        #Cat with first MM file 
        cat([tempList[0]+'MM', tempList[0]+'.bordersMM'], tempList[0]+'.catborders')
        #Fix name
        os.rename(tempList[0]+'.catborders', tempList[0]+'MM')
        #Cat with first UM file
        cat([tempList[0]+'UM', tempList[0]+'.bordersUM'], tempList[0]+'.catborders')
        #Fix name
        os.rename(tempList[0]+".catborders", tempList[0]+"UM")
       
        
        #PHASE II: Parsing multi-mapped reads
        output = []
        output = Parallel(n_jobs=thr)(delayed(allocation.parseMultiPE)(i, seed, winSize, genLand, m, cnn_scores, rc, keep, rmz, maxa) for i in tempList)
        
        
        #Final parsing and file clean up
        #Deleting temporary files and join allocated files
        l = ["header" + ids]
        for i in tempList:
            if os.path.exists(str(i)+"AL"):
                l.append(i + "AL")
            l.append(i + "UM")
        cat(l, outfile)
        
        for filename in glob.glob('allo.temp.*'):
            os.remove(filename)
        
        #Remove all temporary files
        os.remove('header'+ids)
