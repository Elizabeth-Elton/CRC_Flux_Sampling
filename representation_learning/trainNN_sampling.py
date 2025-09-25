'''
train a neural network for representation learning

simulations should be formatted in a numpy array where the size of the
first dimension is the number of simulations

TRAINING METHODS
----------------
- triplet loss 
    - need to load labels
        labels is a 1d array corresponding to the simulations
        each simulation is labeled based on which group it corresponds to
        these would be the NN outputs if doing classification
        should NOT be one-hot-encoded. Simply number labels from 0 to nGroups-1
    - set training = 'triplet'
    - does NOT require distortion in DataGenerator(...)

- NT-Xent loss (from SimCLR)
    - does NOT require labels for DataGenerator(...)
    - requires distortion for DataGenerator(...)
        the default of 0.01 is probably sufficient
    set training = 'simclr'


NEURAL NETWORK TYPES
-------------------------
- for a basic multi-layer-perceptron, params = {..., 'structure': 'ffn', ...}
- for 1-dimensional convolutional network, params = {..., 'structure': 'c1d', ...}
- for 2-dimensional convolutional network, params = {..., 'structure': 'c2d', ...}

PARAMS
------
- cnnFilters
    - number of elements is the number of convolutional layers
    - the integer value of each element is the number of features extracted in each conv layer
- kernelSize is the size of the conv kernel
- fullyConnected
    - number of elements is the number of layers in the mlp
    - the integer value of each element is the number of neurons in each layer
    - the final element is the projection layer. If visualizing the representations, use 2 or 3
- hiddenActivation
    - activation function for the conv layers and mlp layers
- outputActivation
    - activation for the final projection layers. Leave as linear
- dropout
    - probability of randomly turning off a neuron connection during each training step
- learning_rate
    - the learning rate
- patience
    - number of epochs without improvement before stopping early
    - leaving this higher generally leads to a better training, but it takes longer
- batchNorm
    - True or False on using batchNormalization

'''

from genNN import createTrainedModel
from data_generator import DataGenerator
import os
import numpy as np

saveFld = 'data/results_samplingData'  # Folder to save results
n_nn = 25  # Neural network identifier

sims_list = []  # To temporarily store individual arrays
labels_list = []  # To temporarily store individual label arrays

# List of model names
model_names = [
    'base_kras_caf', 'caf_80_HPI', 'caf_80_HK', 'caf_60_GAPDH', 'caf_60_ALDO', 
    'caf_40_LDH', 'caf_20_LDH', 'caf_100_HPI', 'caf_100_G6PDH', 
    'caf_100_PYK', 'caf_100_ENO', 'caf_100_PGAM', 'caf_100_HK',
    'base_kras_crc', 'crc_80_HPI', 'crc_80_HK', 'crc_60_GAPDH', 'crc_60_ALDO', 
    'crc_40_LDH', 'crc_20_LDH', 'crc_100_HPI', 'crc_100_G6PDH', 
    'crc_100_PYK', 'crc_100_ENO', 'crc_100_PGAM', 'crc_100_HK'
]

# Assign unique numeric labels to each model
z = 0
for model_name in model_names:
    # File path for the current model
    file_path = f'data/sampling_results_rp/{model_name}.csv'
    # Load the data from the CSV file
    s = np.loadtxt(file_path, delimiter=',')  # Ensure the file exists and is formatted correctly
    
    # Debugging information
    print(f"Processing file: {model_name}, Shape: {s.shape}")
    
    # Create a label array for the current model
    label_array = np.full((s.shape[0],), z)  # Assign the same label to all simulations of this model
    labels_list.append(label_array)
    sims_list.append(s)
    
    z += 1  # Increment the label for the next model

# Combine all simulations and labels into single arrays
sims = np.vstack(sims_list)
labels = np.concatenate(labels_list)

# Debugging information
print(f"Total shape of sims: {sims.shape}")
print(f"Total shape of labels: {labels.shape}")
#print(f"Example labels: {labels[:5]}")  # Print a few example labels

# Check for mismatch between simulations and labels
if sims.shape[0] != labels.shape[0]:
    exit('Error: Mismatch between the number of simulations and labels.')

# Proceed with training
training = 'triplet'  # Use triplet loss for training
input_shape = sims[0].shape
batchSize = 15  # Adjust batch size as needed
data = DataGenerator(sims, batchSize, training, distortion=0.01, labels=labels)

# Define the parameters for the neural network
params = {
    'input_shape': input_shape,
    'structure': 'ffn',
    'fullyConnected': [64, 2],  # Adjust layers as needed
    'hiddenActivation': 'selu',
    'outputActivation': 'linear',
    'dropout': 0.0,
    'learning_rate': 0.5E-4,  # Default 
    'epochs': 1000,
    'patience': 1000,
    'batchNorm': False,
    'training': training
}

print("\nStarting training...\n")
projector, loss = createTrainedModel(params, data)

# Save the trained model and loss
os.makedirs(f'{saveFld}/projector', exist_ok=True)
os.makedirs(f'{saveFld}/loss', exist_ok=True)
projector.save(f'{saveFld}/projector/model_{n_nn}.keras')
np.savetxt(f'{saveFld}/loss/loss_{n_nn}.csv', loss, delimiter=',')