import tensorflow as tf 
from tensorflow import keras
from tensorflow.keras import layers
import numpy as np

def genFFN(params):
    '''
    creates a convolutional neural network

    cnnFilters - list of integers representing the number of neurons in each convolutional layer
              - length of list = number of layers
    kernelSize - list of tuples of kernel/filter sizes for each layer
    fullyConnected - list of integers representing the number of neurons in each fully connected layer
    '''

    input_shape = params['input_shape']
    fullyConnected = params['fullyConnected']
    hiddenActivation = params['hiddenActivation']
    outputActivation = params['outputActivation']
    dropout = params['dropout']
    batchNorm = params['batchNorm']

    nn_input = layers.Input(input_shape)
    x = layers.Dense(fullyConnected[0], activation=hiddenActivation)(nn_input)
    if batchNorm == True:
        x = layers.BatchNormalization()(x)
    x = layers.Dropout(dropout)(x)
    for i in range(1,len(fullyConnected)-1):
        x = layers.Dense(fullyConnected[i], activation=hiddenActivation)(x)
        if batchNorm == True:
            x = layers.BatchNormalization()(x)
        x = layers.Dropout(dropout)(x)
    x = layers.Dense(fullyConnected[-1], activation=outputActivation)(x)
    embedding_network = keras.Model(nn_input, x)

    return embedding_network 


def genC1D(params):
    '''
    creates a 1D convolutional neural network

    cnnFilters - list of integers representing the number of neurons in each convolutional layer
              - length of list = number of layers
    kernelSize - list of tuples of kernel/filter sizes for each layer
    fullyConnected - list of integers representing the number of neurons in each fully connected layer
    '''

    input_shape = params['input_shape']
    cnnFilters = params['cnnFilters']
    kernelSize = params['kernelSize']
    fullyConnected = params['fullyConnected']
    hiddenActivation = params['hiddenActivation']
    outputActivation = params['outputActivation']
    dropout = params['dropout']
    batchNorm = params['batchNorm']

    nn_input = layers.Input(input_shape)
    x = layers.Conv1D(cnnFilters[0], kernel_size=kernelSize, activation=hiddenActivation)(nn_input)
    x = layers.MaxPooling1D()(x)
    for i in range(1, len(cnnFilters)):
        x = layers.Conv1D(cnnFilters[i], kernel_size=kernelSize, activation=hiddenActivation)(x)
        x = layers.MaxPooling1D()(x)
    x = layers.Flatten()(x)
    for i in range(1,len(fullyConnected)-1):
        x = layers.Dense(fullyConnected[i], activation=hiddenActivation)(x)
        if batchNorm == True:
            x = layers.BatchNormalization()(x)
        x = layers.Dropout(dropout)(x)
    x = layers.Dense(fullyConnected[-1], activation=outputActivation)(x)
    embedding_network = keras.Model(nn_input, x)

    return embedding_network 


def genC2D(params):
    '''
    creates a 2D convolutional neural network

    cnnFilters - list of integers representing the number of neurons in each convolutional layer
              - length of list = number of layers
    kernelSize - list of tuples of kernel/filter sizes for each layer
    fullyConnected - list of integers representing the number of neurons in each fully connected layer
    '''

    input_shape = params['input_shape']
    cnnFilters = params['cnnFilters']
    kernelSize = params['kernelSize']
    fullyConnected = params['fullyConnected']
    hiddenActivation = params['hiddenActivation']
    outputActivation = params['outputActivation']
    dropout = params['dropout']
    batchNorm = params['batchNorm']

    nn_input = layers.Input(input_shape)
    x = layers.Conv2D(cnnFilters[0], kernel_size=kernelSize, activation=hiddenActivation)(nn_input)
    x = layers.MaxPooling2D()(x)
    for i in range(1, len(cnnFilters)):
        x = layers.Conv2D(cnnFilters[i], kernel_size=kernelSize, activation=hiddenActivation)(x)
        x = layers.MaxPooling2D()(x)
    x = layers.Flatten()(x)
    for i in range(1,len(fullyConnected)-1):
        x = layers.Dense(fullyConnected[i], activation=hiddenActivation)(x)
        if batchNorm == True:
            x = layers.BatchNormalization()(x)
        x = layers.Dropout(dropout)(x)
    x = layers.Dense(fullyConnected[-1], activation=outputActivation)(x)
    embedding_network = keras.Model(nn_input, x)

    return embedding_network 


def cosine_similarity(p1, p2):
    p1 = tf.math.l2_normalize(p1, axis=1)
    p2 = tf.math.l2_normalize(p2, axis=1)
    return tf.matmul(p1, p2, transpose_b=True)

def euclidean_similarity(p1, p2):
    # Get the dot product between all embeddings
    # shape (batch_size, batch_size)
    dot_product = tf.matmul(p1, tf.transpose(p2))

    # Get squared L2 norm for each embedding. We can just take the diagonal of `dot_product`.
    # This also provides more numerical stability (the diagonal of the result will be exactly 0).
    # shape (batch_size,)
    square_norm = tf.linalg.diag_part(dot_product)

    # Compute the pairwise distance matrix as we have:
    # ||a - b||^2 = ||a||^2  - 2 <a, b> + ||b||^2
    # shape (batch_size, batch_size)
    distances = tf.expand_dims(square_norm, 1) - 2.0 * dot_product + tf.expand_dims(square_norm, 0)

    # Because of computation errors, some distances might be negative so we put everything >= 0.0
    distances = tf.maximum(distances, 0.0)

    # Because the gradient of sqrt is infinite when distances == 0.0 (ex: on the diagonal)
    # we need to add a small epsilon where distances == 0.0
    mask = tf.cast(tf.equal(distances, 0.0), tf.float32)
    distances = distances + mask * 1e-16
    #distances = distances + 1e-16

    distances = tf.sqrt(distances)

    # Correct the epsilon added: set the distances on the mask to be exactly 0.0
    distances = distances * (1.0 - mask)

    # https://stats.stackexchange.com/questions/53068/euclidean-distance-score-and-similarity
    similarity = 1/(1 + distances)
    #similarity = 2*similarity - 1
    #distances = distances/tf.reduce_max(distances)
    #similarity = 1 - distances
    #similarity = 1/(1+distances)
    #similarity = 1/tf.math.exp(distances)

    return similarity


def contrastive_loss(projections, temperature=0.1):
    '''
    NT-Xent loss
    '''

    projections_1, projections_2 = projections
    similarities = (euclidean_similarity(projections_1, projections_2)/temperature)

    batchSize = tf.shape(projections_1)[0]
    contrastiveLabels = tf.range(batchSize)
    loss_1_2 = keras.losses.sparse_categorical_crossentropy(
        contrastiveLabels, similarities, from_logits=True
    )
    loss_2_1 = keras.losses.sparse_categorical_crossentropy(
        contrastiveLabels, tf.transpose(similarities), from_logits=True
    )

    return (loss_1_2 + loss_2_1)/2


def identity_loss(y_true, y_pred):
    return tf.math.reduce_mean(y_pred)


def genSCLR(params):
    if params['structure'] == 'ffn':
        embedding_network = genFFN(params)
    elif params['structure'] == 'c1d':
        embedding_network = genC1D(params)
    elif params['structure'] == 'c2d':
        embedding_network == genC2D(params)
    else:
        exit('wrong nn structure')
        
    input_shape = params['input_shape']

    input_a = layers.Input(input_shape)
    input_b = layers.Input(input_shape)

    network_a = embedding_network(input_a)
    network_b = embedding_network(input_b)

    loss = layers.Lambda(contrastive_loss)([network_a, network_b])
    simclr = keras.Model(inputs=[input_a, input_b], outputs=loss)

    simclr.compile(loss=identity_loss, optimizer=tf.keras.optimizers.Adam(learning_rate=params['learning_rate']))

    projector = keras.Model(inputs=input_a, outputs=network_a)
    projector.compile()

    return simclr, projector


def euclidean_distance(vects):
    x, y = vects
    sum_square = tf.math.reduce_sum(tf.math.square(x - y), axis=1, keepdims=True)
    return tf.math.sqrt(sum_square)


def triplet_loss(vects, margin=1):
    anchor, positive, negative = vects
    pos_dist = tf.math.reduce_sum(tf.math.square(anchor - positive), axis=1)
    neg_dist = tf.math.reduce_sum(tf.math.square(anchor - negative), axis=1)
    basic_loss = pos_dist - neg_dist + margin 
    loss = tf.math.maximum(basic_loss, 0.0)
    return loss


def genTL(params):
    if params['structure'] == 'ffn':
        embedding_network = genFFN(params)
    elif params['structure'] == 'c1d':
        embedding_network = genC1D(params)
    elif params['structure'] == 'c2d':
        embedding_network == genC2D(params)
    else:
        exit('wrong nn structure')
        
    learning_rate = params['learning_rate']
    input_shape = params['input_shape']

    input_a = layers.Input(input_shape)
    input_p = layers.Input(input_shape)
    input_n = layers.Input(input_shape)

    network_a = embedding_network(input_a)
    network_p = embedding_network(input_p)
    network_n = embedding_network(input_n)

    trip_loss = layers.Lambda(triplet_loss)([network_a, network_p, network_n])
    trainer = keras.Model(inputs=[input_a, input_p, input_n], outputs=trip_loss)
    trainer.compile(loss=identity_loss, optimizer=tf.keras.optimizers.Adam(learning_rate=learning_rate))

    projector = keras.Model(inputs=input_a, outputs=network_a)
    projector.compile()

    return trainer, projector


def createTrainedModel(params, data):
    if params['training'] == 'simclr':
        trainer, projector = genSCLR(params)
    elif params['training'] == 'triplet':
        trainer, projector = genTL(params)
    else:
        print()
        exit('genNN:createTrainedModel - wrong training type')

    stopEarly = tf.keras.callbacks.EarlyStopping(monitor='loss', patience=params['patience'])
    print('...fitting')
    history = trainer.fit(data,
                  epochs=params['epochs'],
                  shuffle=True,
                  callbacks=[stopEarly],
                  verbose=1)
    loss = history.history['loss']

    return projector, loss

def loadModel(fld):
    return keras.models.load_model(fld)


