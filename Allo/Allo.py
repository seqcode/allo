#Lexi Morrissey, Mahony Lab @ Pennsylvania State University
#Last updated 03.06.2023
#Main method for Allo. Splits the sam files up and sends them to allocation procedure via multiprocessing package.

#Arguments
import argparse
parser = argparse.ArgumentParser(prog = 'Allo', \
                    description = 'Allo is a software that allocates multi-mapped reads in ChIP-seq data\n' \
                    'Developed by Alexis Morrissey, Mahony Laboratory @ The Pennsylvania State University', \
                    epilog= 'For more pre-processing info and basic usage, please visit https://github.com/seqcode/allo')
parser.add_argument('input')
parser.add_argument('-seq', type=str, nargs=1, help='Single-end or paired-end sequencing mode', \
                   choices=['pe','se'], required=True, dest='seq')
parser.add_argument('-o', type=str, nargs=1, help='Output file name', dest='outfile', default=None)
parser.add_argument('--mixed', help='Use CNN trained on a dataset with mixed peaks, narrow by default', action='store_true', default=None)
parser.add_argument('-t', type=int, nargs=1, help='Number of threads, 1 by default', dest='threads', default=None)
parser.add_argument('-max', type=int, nargs=1, help='Maximum value for number of locations a read can map', dest='maxlocations', default=None)
parser.add_argument('--keep-unmap', help='Keep unmapped reads and reads that include N in their sequence', action='store_true', default=None)
parser.add_argument('--remove-zeros', help='Disregard multi-mapped reads that map to regions with 0 uniquely mapped reads (random assignment)', action='store_true', default=None)
parser.add_argument('--r2', help='Use read 2 for allocation procedure instead of read 1 (only for paired-end sequencing)', action='store_true', default=None)
parser.add_argument('--readcount', help='CNN will not be used in allocation, only read counts', action='store_true', default=None)
parser.add_argument('--random', help='Reads will be randomly assigned (similar to Bowtie and BWA on default)', action='store_true', default=None)
parser.add_argument('--ignore', help='Ignore warnings about read sorting', action='store_true', default=None)
parser.add_argument('--parser', help='Ignore warnings about read sorting', action='store_true', default=None)
args = parser.parse_args()

#Imports
import sys
from joblib import Parallel, delayed
import multiprocessing
import subprocess
import os
import allocation
import math
from random import randint
import shutil
import glob
import pysam
import shutil

#Function for concatenating files
def cat(files, outname):
    with open(outname,'wb') as wfd:
        for f in files:
            with open(f,'rb') as fd:
                shutil.copyfileobj(fd, wfd)

##Main Method##          
if __name__ == '__main__':
    
    
    #Make a folder to store all temp files in allo
    ids = str(randint(0, 10000))
    curr_dir = os.getcwd()
    allo_dir = os.path.join(curr_dir, r'allo.'+ids)
    if not os.path.exists(allo_dir):
        os.makedirs(allo_dir)
    
    #Input name and bam conversion if necessary
    try:
        samfile = args.input
    except:
        print("SAM/BAM file not given or PySam does not recognize format of file.")
    #Convert from bam to sam if necessary
    if not ".sam" in samfile and not ".SAM" in samfile:
        pysam.set_verbosity(0) #Pysam outputs a warning about indexing
        infile = pysam.AlignmentFile(samfile, "rb")
        outfile = pysam.AlignmentFile(allo_dir+"/"+ids+".sam", "w", template=infile)
        for rec in infile:
            outfile.write(rec)
        samfile = allo_dir+"/"+ids+".sam"
    #Output name
    if args.outfile is not None:
        outfile = args.outfile[0]
    else:
        outfile = args.input[:-4] + ".allo.sam"
    #Model and window size
    if args.mixed is not None:
        m = os.path.dirname(sys.argv[0])+"/mixed"
        winSize = 500
    else:
        m = os.path.dirname(sys.argv[0])+"/narrow"
        winSize = 500
    #Number of threads
    if args.threads is not None:
        thr = args.threads[0]
    else:
        thr = 1
    #Single or paired end
    if args.seq[0] == "se":
        seq = 1
    elif args.seq[0] == "pe":
        seq = 0
    #Options for probabilities
    if args.readcount is not None:
        rc = 1
    elif args.random is not None:
        rc = 2
    else:
        rc = 0
    #Keep unmapped reads
    if args.keep_unmap is not None:
        keep = 1
    else:
        keep = 0
    #Remove reads that only map to zero regions
    if args.remove_zeros is not None:
        rmz = 1
    else:
        rmz = 0
    #User defined max number of alignments
    if args.maxlocations is not None:
        maxa = args.maxlocations[0]
    else:
        maxa = None
    #Use read 2 instead of read 1 for allocation
    if args.r2 is not None:
        r2 = 1
    else:
        r2 = 0
    #Ignore allo warnings
    if args.ignore is not None:
        ig = 1
    else:
        ig = 0
    #Use allo's parser function only
    if args.parser is not None:
        parse = 1
    else:
        parse = 0
    
    
    #Getting header information and checking to see if file was sorted
    collate = 0
    f = open(samfile)
    header = open(allo_dir+'/header', 'a')
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
        print("Warning: PG tag with collate not present in header. File does not appear to be sorted.")
        print("Please use Samtools collate before running Allo. To ignore this warning use argument --ignore")
        sys.exit(0)
    elif collate == 0 and ig == 1:
        print("Suppressed warning: PG tag not present in header. File does not appear to be sorted.")

    
    #Need list of file names
    tempList = []
    for i in range (0,thr):
        tempList.append(allo_dir+"/allo.temp." + ids + "." + str(i).zfill(2))
        
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
        output = Parallel(n_jobs=thr)(delayed(allocation.parseUniq)(i, winSize, cnn_scores, AS, rc, keep, parse) for i in tempList)
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
            
        borderDict = allocation.parseUniq(tempList[0]+".borders", winSize, cnn_scores, AS, rc, keep, parse)
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
        
        #Stopping if parsing only
        if parse == 1:
            l = [allo_dir+"/header"]
            for i in tempList:
                if os.path.exists(str(i)+"UM"):
                    l.append(i + "UM")
            cat(l, args.input[:-4]+".allo.UMR_only.sam")
            l = [allo_dir+"/header"]
            for i in tempList:
                if os.path.exists(str(i)+"MM"):
                    l.append(i + "MM")
            cat(l, args.input[:-4]+".allo.MMR_only.sam")
            shutil.rmtree(allo_dir)
            sys.exit(0)
        
        #PHASE II: Parsing multi-mapped reads
        output = []
        output = Parallel(n_jobs=thr)(delayed(allocation.parseMulti)(i, winSize, genLand, m, cnn_scores, rc, keep, rmz, maxa) for i in tempList)
        
        
        #Final parsing and file clean up
        #Deleting temporary files and join allocated files
        l = [allo_dir+"/header"]
        for i in tempList:
            if os.path.exists(str(i)+"AL"):
                l.append(i + "AL")
            l.append(i + "UM")
        cat(l, outfile)
        
        shutil.rmtree(allo_dir)

        
    
    ##Paired-end mode##
    if seq == 0:
    
        #PHASE I: Parsing unique reads and putting into dictionary
        cnn_scores = {} #Used to keep some previous scores for places in the genome by the cnn to save time
        output = Parallel(n_jobs=thr)(delayed(allocation.parseUniqPE)(i, winSize, cnn_scores, AS, rc, keep, r2, parse) for i in tempList)
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
        borderDict = allocation.parseUniqPE(tempList[0]+".borders", winSize, cnn_scores, AS, rc, keep, r2, parse)
        d.update(borderDict)
        #Cat with first MM file 
        cat([tempList[0]+'MM', tempList[0]+'.bordersMM'], tempList[0]+'.catborders')
        #Fix name
        os.rename(tempList[0]+'.catborders', tempList[0]+'MM')
        #Cat with first UM file
        cat([tempList[0]+'UM', tempList[0]+'.bordersUM'], tempList[0]+'.catborders')
        #Fix name
        os.rename(tempList[0]+".catborders", tempList[0]+"UM")
        
        #Stopping if parsing only
        if parse == 1:
            l = [allo_dir+"/header"]
            for i in tempList:
                if os.path.exists(str(i)+"UM"):
                    l.append(i + "UM")
            cat(l, args.input[:-4]+".allo.UMR_only.sam")
            l = [allo_dir+"/header"]
            for i in tempList:
                if os.path.exists(str(i)+"MM"):
                    l.append(i + "MM")
            cat(l, args.input[:-4]+".allo.MMR_only.sam")
            shutil.rmtree(allo_dir)
            sys.exit(0)
        
        #PHASE II: Parsing multi-mapped reads
        output = []
        output = Parallel(n_jobs=thr)(delayed(allocation.parseMultiPE)(i, winSize, genLand, m, cnn_scores, rc, keep, rmz, maxa) for i in tempList)
        
        
        #Final parsing and file clean up
        #Deleting temporary files and join allocated files
        l = [allo_dir+"/header"]
        for i in tempList:
            if os.path.exists(str(i)+"AL"):
                l.append(i + "AL")
            l.append(i + "UM")
        cat(l, outfile)
        
        shutil.rmtree(allo_dir)
