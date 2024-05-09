#!/usr/bin/env python
#Lexi Morrissey, Mahony Lab @ Pennsylvania State University
#Contains method for predicting whether area should receive multimapped reads via pre-trained CNN in Allo.

import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.models import load_model, Model
import os
tf.config.run_functions_eagerly(False)
import math
import sys
import numpy as np


def predictNN(counts, winSize, model):
    pic = np.zeros((100, 100), float)
    if sum(counts) == 0:
        pred = model(pic.reshape(-1,100,100), training=False)
        return pred.numpy()[0][0]
    binSize = math.floor(winSize/100)
    binTotal = math.floor(winSize/binSize)
    binned = []
    #Binning
    for i in range (0,binTotal):
        position = i*binSize
        if i == binTotal-1:
            binned.append(np.sum(counts[position:len(counts)]))
        else:
            binned.append(np.sum(counts[position:position+binSize]))
  
    binned = (99*(binned - np.min(binned))/np.ptp(binned)).astype(int) 
    for i in range(0,len(binned)):
        pic[binned[i],i] = 1
    pred = model.predict(pic.reshape(-1,100,100)) 
    try:
        #pred = model(pic.reshape(-1,100,100), training=False)
        pred = model.predict(pic.reshape(-1,100,100))
        #return pred.numpy()[0][0]
        return pred[0][0]

    except:
        print("Could not predict with Tensorflow model :( Allo was written with Tensorflow version 2.11", flush=True)
        sys.exit(0)
