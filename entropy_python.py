"""
Python code that calls c++ shared object files
Postrual control entropy analsysis with Santos Dataset
100 Hz and COP data in m
"""
import ctypes
from pandas import read_csv
import wfdb
import matplotlib.pyplot as plt
from pandas import DataFrame, read_csv, concat
from tqdm import tqdm
from numpy import ones, zeros

class fastEntropy():
    def __init__(self):
        print('Initialize fastEntropy class and call libEntropy')
        # Load the shared library
        shannon_entropy_lib = ctypes.CDLL('./libEntropy.so')
        # Declare the C function
        self.shannon_entropy_func = shannon_entropy_lib.calculate_shannon_entropy 
        self.sample_entropy_func = shannon_entropy_lib.calculate_sample_entropy
        self.approximate_entropy_func = shannon_entropy_lib.calculate_approximate_entropy
        self.shannon_entropy_func.restype = ctypes.c_double
        self.sample_entropy_func.restype = ctypes.c_double
        self.approximate_entropy_func.restype = ctypes.c_double

    def get_bds_info(self,path="/home/brianszekely/Desktop/ProjectsResearch/entropy_c++_python/data/BDSinfo.txt"):
        self.bds_info =  read_csv(path, sep='\t', header=0)
        print(self.bds_info)

    def extract_trials(self):
        """
        1. Extract trials that are open and closed eyes static two feet postural control
        2. Extract trials that have the same shoe type - minimal
        3. No disease 
        4. Healthy young and Healthy old

        ['Flip-Flops' 'Walking shoes, Flip-Flops' 'Boots' 'Walking shoes'
        'Walking shoes with silicone insole' 'Casual shoes''Walking shoes, Flip-Flops, Casual shoe
        'Walking shoes, Casual shoes' s'
        'Walking shoes, Boots' 'Sandal' 'Dress shoes' 'Casual shoes, Dress shoes'
        'Casual shoes, Flip-Flops' 'Tennis with zero drop'
        'Walking shoes, Sandal' 'Sandal, Casual shoes' 'Casual shoes, Sandal'
        'Orthopedic shoes']

        """
        #Get all Open and Closed Trials
        filtered_df = self.bds_info[self.bds_info['Vision'].isin(['Open', 'Closed'])]
        #Get all Firm surfaces
        filtered_df = filtered_df[filtered_df['Surface'].isin(['Firm'])]
        #Remove any column that has Best and IPAQ in it
        filtered_columns = [col for col in filtered_df.columns if 'Best' not in col]
        filtered_df = filtered_df[filtered_columns]
        filtered_columns = [col for col in filtered_df.columns if 'IPAQ' not in col]
        filtered_df = filtered_df[filtered_columns]
        filtered_columns = [col for col in filtered_df.columns if 'FES' not in col]
        filtered_df = filtered_df[filtered_columns]
        filtered_columns = [col for col in filtered_df.columns if 'TMT' not in col]
        filtered_df = filtered_df[filtered_columns]
        filtered_df.drop(columns=['Nationality', 'SkinColor'],inplace=True)
        # Subset the DataFrame with the filtered columns
        #only include trials with minimal footwear
        filtered_df = filtered_df[filtered_df['Footwear'].isin(['Flip-Flops', 'Sandal',
                                                                'Walking shoes, Flip-Flops',
                                                                'Casual shoes, Flip-Flops',
                                                                'Sandal, Casual shoes',
                                                                'Casual shoes, Sandal',
                                                                'Walking shoes, Sandal'])]
        #Change str to int
        filtered_df['AgeGroup'] = filtered_df['AgeGroup'].replace({'Young': 0, 'Old': 1})
        self.saved_trials = [filtered_df['Trial'].values, filtered_df['AgeGroup'].values]

    def get_data_participants(self,name):
        """
        example: name="BDS00002"
        """
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
        return self.shannon_entropy_func(data_array_x, len(datax)),self.shannon_entropy_func(data_array_y, len(datay))

    def sample_entropy(self):
        datax = self.cop['copX'].to_list()
        datay = self.cop['copY'].to_list()
        #Convert the input to C array
        data_array_x = (ctypes.c_double * len(datax))(*datax)
        data_array_y = (ctypes.c_double * len(datay))(*datay)
        return self.sample_entropy_func(data_array_x, len(datax)), self.sample_entropy_func(data_array_y, len(datay))

    def approximate_entropy(self):
        datax = self.cop['copX'].to_list()
        datay = self.cop['copY'].to_list()
        #Convert the input to C array
        data_array_x = (ctypes.c_double * len(datax))(*datax)
        data_array_y = (ctypes.c_double * len(datay))(*datay)
        return self.approximate_entropy_func(data_array_x, len(datax)), self.approximate_entropy_func(data_array_y, len(datay))

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
    
    def run_analysis_all_trials(self):
        self.sampen_x_young, self.sampen_y_young, self.sh_x_young, self.sh_y_young, self.apen_x_young, self.apen_y_young  = [], [], [], [], [], []
        self.sampen_x_old, self.sampen_y_old, self.sh_x_old, self.sh_y_old, self.apen_x_old, self.apen_y_old  = [], [], [], [], [], []
        for trial, age in tqdm(zip(*self.saved_trials)):
            self.get_data_participants(trial)
            shannon_x, shannon_y = self.shannon_ent()
            samp_x, sampy = self.sample_entropy()
            ap_x, ap_y = self.approximate_entropy()
            #save in list
            if age == 0:
                self.sh_x_young.append(shannon_x)
                self.sh_y_young.append(shannon_y)
                self.sampen_x_young.append(samp_x)
                self.sampen_y_young.append(sampy)
                self.apen_x_young.append(ap_x)
                self.apen_y_young.append(ap_y)
            elif age == 1:
                self.sh_x_old.append(shannon_x)
                self.sh_y_old.append(shannon_y)
                self.sampen_x_old.append(samp_x)
                self.sampen_y_old.append(sampy)
                self.apen_x_old.append(ap_x)
                self.apen_y_old.append(ap_y)

        young_df = DataFrame({"shannon_x":self.sh_x_young,
                   "shannon_y":self.sh_y_young,
                   "sampen_x":self.sampen_x_young,
                   "sampen_y":self.sampen_y_young,
                   "apen_x":self.apen_x_young,
                   "apen_y":self.apen_y_young,
                   "label":zeros(len(self.sh_x_young))})
        old_df = DataFrame({"shannon_x":self.sh_x_old,
                   "shannon_y":self.sh_y_old,
                   "sampen_x":self.sampen_x_old,
                   "sampen_y":self.sampen_y_old,
                   "apen_x":self.apen_x_old,
                   "apen_y":self.apen_y_old,
                   "label":ones(len(self.sh_x_old))})
        ml_df = concat([young_df, old_df]) 
        ml_df.to_csv('ml_df.csv',index=False)

def main():
    fast = fastEntropy()
    fast.get_bds_info()
    fast.extract_trials()
    fast.run_analysis_all_trials()
    # fast.get_data_participants()
    # fast.shannon_ent()
    # fast.sample_entropy()
    # fast.plot_cop()

if __name__ == "__main__":
    main()