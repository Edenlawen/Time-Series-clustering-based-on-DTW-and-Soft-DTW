"""Modified code from Time-series Generative Adversarial Networks (TimeGAN) Codebase.

Reference: Jinsung Yoon, Daniel Jarrett, Mihaela van der Schaar, 
"Time-series Generative Adversarial Networks," 
Neural Information Processing Systems (NeurIPS), 2019.

Paper link: https://papers.nips.cc/paper/8789-time-series-generative-adversarial-networks

Last updated Date: April 24th 2020
Code author: Jinsung Yoon (jsyoon0823@gmail.com)

-----------------------------

data_loading.py

(0) MinMaxScaler: Min Max normalizer
(1) real_data_loading: Load and preprocess real data
"""

# Necessary Packages
import numpy as np


def MinMaxScaler(data):
    """Min Max normalizer.

    Args:
      - data: original data

    Returns:
      - norm_data: normalized data
    """
    numerator = data - np.min(data, 0)
    denominator = np.max(data, 0) - np.min(data, 0)
    norm_data = numerator / (denominator + 1e-7)
    return norm_data


def real_data_loading(array, seq_len):
    """Load and preprocess real-world datasets.

    Args:
      - array: array
      - seq_len: sequence length

    Returns:
      - data: preprocessed data.
      - max: The maximum value
      - min: The minimum value
    """

    max = np.max(array)
    min = np.min(array)
    # Normalize the data
    array = MinMaxScaler(array)

    # Preprocess the dataset
    temp_data = []
    # Cut data by sequence length
    for i in range(0, len(array) - seq_len):
        _x = array[i:i + seq_len]
        temp_data.append(_x)

    # Mix the datasets (to make it similar to i.i.d)
    idx = np.random.permutation(len(temp_data))
    data = []
    for i in range(len(temp_data)):
        data.append(temp_data[idx[i]])
    
    data = np.array(data)

    return data, max, min
