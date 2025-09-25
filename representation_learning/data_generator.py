import numpy as np
from tensorflow import keras
import random

class DataGenerator(keras.utils.Sequence):
    def __init__(self, sims, batchSize, training, distortion=0.01, labels=[]):
        self._training = training
        self._simSize = sims[0].shape
        self._sims = sims
        self._nSims = sims.shape[0]
        self._batchSize = batchSize
        self._drawOrder = []
        for i in range(self._sims.shape[0]):
                self._drawOrder.append(i)
                
        if self._training == 'simclr':
            self._distortion = distortion

        elif self._training == 'triplet':
            self._labels = labels
            self._groups = []
            for g in range(np.max(self._labels)+1):
                group = []
                for i in range(self._sims.shape[0]):
                    if self._labels[i] == g:
                        group.append(i)
                self._groups.append(group)

        else:
            print()
            exit('data_generator:init - Wrong training type')
        
        self.on_epoch_end()

    def __len__(self):
        return int(np.floor(self._nSims/self._batchSize))

    def __getitem__(self, index):
        # Generate data
        #X, y = self.__data_generation(index)

        #return X, y
        return self.__data_generation(index)

    def on_epoch_end(self):
        random.shuffle(self._drawOrder)

    def __data_generation(self, index):
        '''
        X = augmentations of simulations to be projected
        y = labels -> augmentations of the same simulation get the same label
        '''
        if self._training == 'simclr':
            dataSize = [self._batchSize]
            for i in self._simSize:
                dataSize.append(i)
            X = [np.zeros((dataSize)) for _ in range(2)]
            y = np.zeros((self._batchSize,))+1
            for i in range(self._batchSize):
                idx = self._drawOrder[index*self._batchSize + i]
                distSims = self.distortSimulation(self._sims[idx])
                X[0][i] = distSims[0]
                X[1][i] = distSims[1]
            return tuple(X), y #EE made X a tuple bc tensorflow was mad i gave it a list
        elif self._training == 'triplet':
            dataSize = [self._batchSize]
            for i in self._simSize:
                dataSize.append(i)
            X = [np.zeros(dataSize) for _ in range(3)]
            y = np.zeros((self._batchSize,)) + 1
            for i in range(self._batchSize):
                a_idx = self._drawOrder[index*self._batchSize + i]
                a_group = self._labels[a_idx]
                p_idx = a_idx
                while p_idx == a_idx:
                    p_idx = np.random.choice(self._groups[a_group])
                n_group = a_group
                while n_group == a_group:
                    n_group = np.random.randint(0,len(self._groups))
                n_idx = np.random.choice(self._groups[n_group])
                X[0][i] = self._sims[a_idx]
                X[1][i] = self._sims[p_idx]
                X[2][i] = self._sims[n_idx]
            return tuple(X), y #EE made X a tuple bc tensorflow was mad i gave it a list
        else:
            exit(':(')


    def distortSimulation(self, sim):
        '''
        Generates two distorted versions of the same simulation
        '''
        dSim0 = sim.copy()
        dSim1 = sim.copy()

        dSim0 = dSim0 + np.random.normal(0, self._distortion, size=dSim0.shape)
        dSim1 = dSim1 + np.random.normal(0, self._distortion, size=dSim1.shape)

        return [dSim0, dSim1]

        
