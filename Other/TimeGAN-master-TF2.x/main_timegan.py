"""Time-series Generative Adversarial Networks (TimeGAN) Codebase.

Reference: Jinsung Yoon, Daniel Jarrett, Mihaela van der Schaar, 
"Time-series Generative Adversarial Networks," 
Neural Information Processing Systems (NeurIPS), 2019.

Paper link: https://papers.nips.cc/paper/8789-time-series-generative-adversarial-networks

Last updated Date: April 24th 2020
Code author: Jinsung Yoon (jsyoon0823@gmail.com)

-----------------------------

main_timegan.py

(1) Import data
(2) Generate synthetic data
(3) Evaluate the performances in three ways
  - Visualization (t-SNE, PCA)
  - Discriminative score
  - Predictive score
"""

# Necessary packages
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from data_loading import real_data_loading, sine_data_generation
from timegan import timegan
import pandas as pd

import argparse
import numpy as np
import tensorflow as tf
from time import time
import warnings
warnings.filterwarnings("ignore")


# 1. TimeGAN model
# from timegan_upgraded import timegan

# 2. Data loading


def main(args):
    """Main function for timeGAN experiments.

    Args:
      - data_name: sine, stock, or energy
      - seq_len: sequence length
      - Network parameters (should be optimized for different datasets)
        - module: gru, lstm, or lstmLN
        - hidden_dim: hidden dimensions
        - num_layer: number of layers
        - iteration: number of training iterations
        - batch_size: the number of samples in each batch
      - metric_iteration: number of iterations for metric computation

    Returns:
      - ori_data: original data
      - generated_data: generated synthetic data
      - metric_results: discriminative and predictive scores
    """
    # Keep min&max for normalization
    # Lire le fichier CSV en sautant la première ligne
    df = pd.read_csv("data/co2R.csv", skiprows=1)

    # Obtenir la valeur maximale
    max_value = float(df.max().values)

    # Obtenir la valeur minimale
    min_value = float(df.min().values)

    # Afficher les résultats
    print("Valeur maximale :", max_value)
    print("Valeur minimale :", min_value)

    # Data loading
    if args.data_name in ['stock', 'energy', 'co2', 'co2V', 'co2R']:
        ori_data = real_data_loading(args.data_name, args.seq_len)
    elif args.data_name == 'sine':
        # Set number of samples and its dimensions
        no, dim = 10000, 5
        ori_data = sine_data_generation(no, args.seq_len, dim)

    ori_data = np.expand_dims(ori_data, 2)

    print(args.data_name + ' dataset is ready.')

    # Synthetic data generation by TimeGAN
    # Set newtork parameters
    parameters = dict()
    parameters['module'] = args.module
    parameters['hidden_dim'] = args.hidden_dim
    parameters['num_layer'] = args.num_layer
    parameters['iterations'] = args.iteration
    parameters['batch_size'] = args.batch_size

    generated_data = timegan(ori_data, parameters)
    print('Finish Synthetic Data Generation')

    ori_data = ori_data * (max_value - min_value) + min_value
    generated_data = generated_data * (max_value - min_value) + min_value

    liste = np.empty((0,))

    for i in range(generated_data.shape[0]):
        temp = ori_data[i]
        temp = np.squeeze(temp)
        temp2 = generated_data[i]
        temp2 = np.squeeze(temp2)
        res = np.corrcoef(temp, temp2)[0, 1]
        if (res >= 0.6):
            print(res)
            liste = np.concatenate((liste, temp2))

    # Spécifier le nom du fichier
    f2 = 'output/generated_data.csv'

    # Utiliser la fonction savetxt pour enregistrer le tableau dans le fichier CSV
    np.savetxt(f2, liste, delimiter=',')

    print("Le tableau a été enregistré dans le fichier CSV avec succès.")
    return ori_data, generated_data


if __name__ == '__main__':

    # Inputs for the main function
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '--data_name',
        choices=['stock', 'energy', 'co2', 'co2V', 'co2R'],
        default='co2R',
        type=str)
    parser.add_argument(
        '--seq_len',
        help='sequence length',
        default=365,
        type=int)
    parser.add_argument(
        '--module',
        choices=['gru', 'lstm', 'lstmLN'],
        default='gru',
        type=str)
    parser.add_argument(
        '--hidden_dim',
        help='hidden state dimensions (should be optimized)',
        default=24,
        type=int)
    parser.add_argument(
        '--num_layer',
        help='number of layers (should be optimized)',
        default=3,
        type=int)
    parser.add_argument(
        '--iteration',
        help='Training iterations (should be optimized)',
        default=10,
        type=int)
    parser.add_argument(
        '--batch_size',
        help='the number of samples in mini-batch (should be optimized)',
        default=64,
        type=int)
    parser.add_argument(
        '--metric_iteration',
        help='iterations of the metric computation',
        default=5,
        type=int)

    args = parser.parse_args()

    # Calls main function
    ori_data, generated_data = main(args)
