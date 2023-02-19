"""
Python code that calls c++ shared object files
Postrual control entropy analsysis with Santos Dataset
100 Hz and COP data in m
"""
import ctypes
from pandas import read_csv
import wfdb
import matplotlib.pyplot as plt
from pandas import DataFrame, read_csv

class fastEntropy():
    def __init__(self):
        print('Initialize fastEntropy class and call libEntropy')
        # Load the shared library
        shannon_entropy_lib = ctypes.CDLL('./libEntropy.so')
        # Declare the C function
        self.shannon_entropy_func = shannon_entropy_lib.calculate_shannon_entropy 
        self.sample_entropy_func = shannon_entropy_lib.calculate_sample_entropy
        self.shannon_entropy_func.restype = ctypes.c_double
        self.sample_entropy_func.restype = ctypes.c_double
    def get_bds_info(self,path="/home/brianszekely/Desktop/ProjectsResearch/entropy_c++_python/data/BDSinfo.txt"):
        return read_csv(path, sep='\t', header=0)
    def extract_trials(self,df):
        """
        1. Extract trials that are open and closed eyes static two feet postural control
        2. Extract trials that have the same show type - minimal
        3. No disease 
        4. Healthy young and Healthy old
        """
        pass
    def get_data_participants(self,name="BDS00002"):
        record = wfdb.rdrecord(name, pn_dir='hbedb/1.0.0/')
        # record = wfdb.rdrecord('BDS00001', pn_dir='hbedb/1.0.0/')
        # signal = record.p_signal
        copX, copY = [], []
        for i in range(len(record.p_signal)):
            copX.append(record.p_signal[i][6]) #COPx
            copY.append(record.p_signal[i][7]) #COPy
        self.cop = DataFrame({"copX":copX,"copY":copY})
        # Mean Offset to orient the data around 0
        self.cop["copX"] = self.cop["copX"] - self.cop["copX"].mean()
        self.cop["copY"] = self.cop["copY"] - self.cop["copY"].mean()
        # # Put data in cm
        self.cop["copX"] = self.cop["copX"] * 100
        self.cop["copY"] = self.cop["copY"] * 100
    def shannon_ent(self):
        datax = self.cop['copX'].to_list()
        datay = self.cop['copY'].to_list()
        #Convert the input to C array
        data_array_x = (ctypes.c_double * len(datax))(*datax)
        data_array_y = (ctypes.c_double * len(datay))(*datay)
        self.sh_ent_x = self.shannon_entropy_func(data_array_x, len(datax))
        self.sh_ent_y = self.shannon_entropy_func(data_array_y, len(datay))
        print(self.sh_ent_x)
        print(self.sh_ent_y)
    def sample_entropy(self):
        datax = self.cop['copX'].to_list()
        datay = self.cop['copY'].to_list()
        #Convert the input to C array
        data_array_x = (ctypes.c_double * len(datax))(*datax)
        data_array_y = (ctypes.c_double * len(datay))(*datay)
        self.samp_en_x = self.sample_entropy_func(data_array_x, len(datax))
        self.samp_en_y = self.sample_entropy_func(data_array_y, len(datay))
        print(self.samp_en_x)
        print(self.samp_en_y)
    def plot_cop(self):
        fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(8, 10))
        ax[0].plot(self.cop['copX'],self.cop['copY'])
        ax[0].set_xlabel('copX')
        ax[0].set_ylabel('copY')
        ax[1].plot(self.cop.index,self.cop['copX'])
        ax[1].set_xlabel('Samples')
        ax[1].set_ylabel('copX')
        ax[2].plot(self.cop.index,self.cop['copY'])
        ax[2].set_xlabel('Samples')
        ax[2].set_ylabel('copY')
        plt.show()

def main():
    fast = fastEntropy()
    fast.get_data_participants()
    fast.shannon_ent()
    fast.sample_entropy()
    fast.plot_cop()
if __name__ == "__main__":
    main()