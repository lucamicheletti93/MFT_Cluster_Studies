import torch
import tensorflow as tf
from tensorflow.keras.models import Model
import matplotlib.pyplot as plt
import seaborn as sn
import pandas as pd
import numpy as np

class AutoEncoder(Model):
    def __init__(self, input_dim):
        super(AutoEncoder, self).__init__()
        self.encoder = tf.keras.Sequential([
                      tf.keras.layers.Dense(input_dim*5, activation="relu"),
                      tf.keras.layers.Dense(input_dim*4, activation="relu"),
                      tf.keras.layers.Dense(input_dim*3, activation="relu"),
                      tf.keras.layers.Dense(input_dim*2, activation="relu"),
                      tf.keras.layers.Dense(input_dim, activation="relu")
                  ])
        self.decoder = tf.keras.Sequential([
                       tf.keras.layers.Dense(input_dim, activation="relu"),
                       tf.keras.layers.Dense(input_dim*2, activation="relu"),
                       tf.keras.layers.Dense(input_dim*3, activation="relu"),
                       tf.keras.layers.Dense(input_dim*4, activation="relu"),
                       tf.keras.layers.Dense(input_dim*5, activation="relu"),
                       tf.keras.layers.Dense(input_dim, activation='relu')
                  ])
    def call(self, x):
        encoded = self.encoder(x)
        decoded = self.decoder(encoded)
        return decoded

    def get_label_pred(model, data, thr):
        reco = model(data)
        loss = tf.keras.losses.mae(data, reco)
        return tf.math.less(loss, thr)

    def get_stats(labels, predictions):
        print(f'Accuracy: {tf.keras.metrics.Accuracy()(labels, predictions)}')
        print(f'Precision: {tf.keras.metrics.Precision()(labels, predictions)}')
        print(f'Recall: {tf.keras.metrics.Recall()(labels, predictions)}')
        
    
class MLplots:
    def _plot_loss_vs_epoch(history, OutputDir):
        fig = plt.figure()
        plt.plot(history.history['loss'], label='train')
        plt.plot(history.history['val_loss'], label='test')
        plt.legend()
        plt.xlabel('Epoch')
        plt.ylabel('Loss')
        plt.title('Loss vs Epoch')
        plt.savefig(f'{OutputDir}/LossVsEpoch.pdf')
        
def plot_confusion_matrix(tp, fp, fn, tn, OutputDir, labels=['Anomaly', 'Normal'], title='Confusion Matrix'):
    df_cm = pd.DataFrame([[tn, fp], [fn, tp]], index=[i for i in labels], columns=[f'{i}Pred' for i in labels])
    fig = plt.figure()
    sn.heatmap(df_cm, cmap='Blues', annot=True, fmt='g')
    plt.title(title)
    plt.savefig(f'{OutputDir}/ConfusionMatrix.pdf')
