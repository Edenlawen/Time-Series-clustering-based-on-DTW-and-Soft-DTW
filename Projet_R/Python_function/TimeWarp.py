import time
import numpy as np
import os
import utils.augmentation as aug
from data_loading import real_data_loading
from EC_compareCourbe import compute_indicateurComp


def TimeWarp(data, nbpts, seq_len, condition, verbose=False):
    """Generates new points from the dataset given as an argument

    Args:
      - data: The original data converted to ndarray for numpy
      - nbpts: The number of points to be generated. The program will stop when it exceeds this number
      - seq_len : The size of the number of points desired when a sequence of points will be generated
      - condition: The number of conditions met to select a new dataset. The greater the number, the longer it may take. The maximum is 8
      - verbose

    Returns:
      - new_data: The original dataset and all those generated
    """
    if verbose:
        print("Start the Augmentation")

    # Preprocess the dataset
    if verbose:
        print("Preprocess the dataset")
    x_train, max, min = real_data_loading(data, (len(data)-1))
    denoR = x_train[0] * (max - min) + min
    denoR = np.expand_dims(denoR, 1)
    x_train = np.expand_dims(x_train[0, :int(seq_len)], 0)
    x_train = np.expand_dims(x_train, 2)

    if verbose:
        print("End of preprocess")

    array = np.empty((1,), dtype=np.object).reshape((1, 1))

    array = np.append(array, denoR, axis=0)

    i = 0
    compt = 0
    while (len(array) < nbpts+len(data)-1):
        if (verbose and i % 1000 == 0):
            print(i)
        time_warp = aug.time_warp(x_train)[0]
        res = compute_indicateurComp(x_train[0], time_warp)[
            'Nombre de conditions vraies']
        if (res >= condition):
            if verbose:
                print("New array found")
                compt += seq_len
                print(str((compt*100)/nbpts) + "%")
            denoR = time_warp * (max - min) + min
            denoR.flatten()
            array = np.append(array, denoR, axis=0)
        i += 1
    array.flatten()
    array = np.delete(array, 0)
    return array
