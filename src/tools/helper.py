import numpy as np

def data_cutoff(dataset, cutoff):
    dataset.y = dataset.y[:np.argmax(dataset.x>cutoff)]
    dataset.y_err = dataset.y_err[:np.argmax(dataset.x>cutoff)]
    dataset.x = dataset.x[:np.argmax(dataset.x>cutoff)]
    return dataset
