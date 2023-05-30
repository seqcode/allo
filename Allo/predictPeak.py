#Lexi Morrissey, Mahony Lab @ Pennsylvania State University
#Last updated 03.01.2023
#Contains method for predicting whether area should receive multimapped reads via pre-trained CNN in Allo.

import os
import numpy as np
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.models import load_model, Model
import os
tf.config.run_functions_eagerly(False)
import math


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
            binned.append(np.sum(counts[position:position+5]))
  
    binned = (99*(binned - np.min(binned))/np.ptp(binned)).astype(int) 
    for i in range(0,len(binned)):
        pic[binned[i],i] = 1
        
    try:
        pred = model(pic.reshape(-1,100,100), training=False)
        return pred.numpy()[0][0]

    except:
        print("Could not predict with Tensorflow model :( Allo was written with Tensorflow version 2.11")
        sys.exit(0)



