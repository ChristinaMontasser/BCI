import math
import numpy as np

class ERP:
    def __init__(self, signal_time_length, baselineTime = 0.5):
        #baselineTime is a precentage from the signal that describes baseline correction 
        self.baselineTime = baselineTime
        self.signal_time_length = signal_time_length

    def baseLine(self, signal):
        self.baseline = int(signal.shape[1]*((self.baselineTime)/(self.signal_time_length)))
        mean = np.mean(signal[:self.baseline], axis=1)
        signal = signal.T-mean
        return signal

    def analyse(self, signal):
        #Sampling rate or time, any of them is nessacry 
        self.signal = signal
        self.erp = []
        for i in range(len(signal)):
            self.signal[i] =self.baseLine(self.signal[i])
            self.erp.append(np.mean(self.signal[i][self.baseline:].T, axis=0))
        return self.erp

