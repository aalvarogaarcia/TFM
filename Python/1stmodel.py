import numpy as np
import tensorflow as tf
import matplotlib.pyplot as plt



class Matrix:
    def __init__(self, list):
        if type(list) is type("str"):
            self.list = np.loadtxt(list, delimiter=',')
        elif type(list) is type(np.zeros(1)):
            self.list = np.array(list)
            
        self.nrow = self.list.shape[0]
        self.ncolumn = self.list.shape[1]  
        
        
    def random_divise(self, prop):
        shuffle_matrix = self.list
        np.random.shuffle(shuffle_matrix)
        N = round(self.nrow*prop)
    
        matrix_1 = Matrix(shuffle_matrix[:N] )
        matrix_2 = Matrix(shuffle_matrix[N:])
        
        print(matrix_1.nrow)
        print(matrix_2.nrow)
        
        return (matrix_1, matrix_2)

class Network:
    def __init__(self, X, Y, loss):
        self.layer = tf.keras.layers.Dense(units = Y.ncolumn, input_shape = [X.ncolumn])
        self.modelo = tf.keras.Sequential([self.layer])
        
        self.modelo.compile(
            optimizer = tf.keras.optimizers.Adam(0.1),
            loss = loss
        )
        
        self.data = X.random_divise(0.3)
        self.label = Y.random_divise(0.3)
        self.train_data = self.data[0]
        self.train_label = self.label[0]
        self.test_data = self.data[1]
        self.test_label = self.label[1]
   
        
    def model_train(self, epochs):
        print('Comenzando entrenamiento...')
        background = self.modelo.fit(self.train_data, self.train_label, epochs = epochs, verbose = False)
        print('¡Entrenamiento completado!')
        return background
        
    def __str__(self):
        plt.xlabel('# epochs')
        plt.ylabel('Loss magnitude')
        background = self.model_train(1000)
        
        plt.plot(background.history['loss'])
    
    
    
    
X = Matrix ('C:/Users/agarm/Desktop/Repository/Data/datamatrix')
Y = Matrix ('C:/Users/agarm/Desktop/Repository/Data/labelmatrix')

network = Network (X, Y, 'mean_squared_error')

print(network)

