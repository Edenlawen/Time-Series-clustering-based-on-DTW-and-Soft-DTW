import numpy as np
import os
import augmentation as aug
from data_loading import real_data_loading
from EC_compareCourbe import compute_indicateurComp


def compute_TimeWarp(data, nbpts):
    """Generates new points from the dataset given as an argument

    Args:
      - data: array
      - nbpts: sequence length
      
    Returns:
      - new_data: preprocessed data.
    """
    print("hello")
    
    
# Loading Data
x_train, max, min = real_data_loading('co2R', 467)

x_train = np.expand_dims(x_train,2)
y_train = np.expand_dims(x_train[0,365:],0)
x_train = np.expand_dims(x_train[0,:365],0)

tot = np.empty((1,), dtype=np.object).reshape((1, 1))

while(len(tot) < 20000) :
    time_warp = aug.time_warp(x_train)[0]
    res = compute_indicateurComp(x_train[0], time_warp)['Nombre de conditions vraies']
    if (res >= 7):
        denoR = time_warp * (max - min) + min
        denoR.flatten()
        tot = np.append(tot, denoR, axis=0)
tot.flatten()
tot = np.delete(tot,0)
