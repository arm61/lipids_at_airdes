import numpy as np


def data_cutoff(dataset, cutoff):
    dataset.y = dataset.y[:np.argmax(dataset.x > cutoff)]
    dataset.y_err = dataset.y_err[:np.argmax(dataset.x > cutoff)]
    dataset.x = dataset.x[:np.argmax(dataset.x > cutoff)]
    return dataset


def latex_asym(a, b, c):
    return '$' + str(a) + '^{+' + str(b) + '}_{-' + str(c) + '}$'

def latex_sym(a, b):
    return '$' + str(a) + '\pm{' + str(b) + '}$'
