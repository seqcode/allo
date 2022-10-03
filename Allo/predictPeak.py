#Lexi Morrissey, Mahony Lab @ Pennsylvania State University
#Last updated 11.08.2021
#Contains method for predicting whether area should receive multimapped reads via pre-trained CNN in Allo.

import time
import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras.models import load_model, Model
import math


def predictNN(countArray, winSize, model):
    #Empty window
    pic = np.zeros((100, 100), float)
    #Bins needed to make image 100x100
    binSize = math.floor(winSize / 100)
    countArray = np.array(countArray)
    #Binning
    binned = countArray.reshape(-1, binSize).mean(axis=1).astype(int).tolist()
    #Making image for non-zero areas
    if not np.sum(binned) == 0:
        binned = (99 * (binned - np.min(binned)) / np.ptp(binned)).astype(int)
        #Making image
        for i in range(0, len(binned)):
            pic[binned[i], i] = 1
    #Predicting with pre-trained CNN
    try:
        pred = model.predict(pic.reshape(-1, 100, 100))
        return pred[0][0]
    except:
        print("Could not predict with Tensorflow model :( Allo was written with Tensorflow version 2.4.1")
        sys.exit(0)



