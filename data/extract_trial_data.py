import wfdb
import matplotlib.pyplot as plt
from pandas import DataFrame, read_csv

def get_bds_info(path="/home/brianszekely/Desktop/ProjectsResearch/entropy_c++_python/data/BDSinfo.txt"):
    return read_csv(path, sep='\t', header=0)
def extract_trials(df):
    """
    1. Extract trials that are open and closed eyes static two feet postural control
    2. Extract trials that have the same show type - minimal
    3. No disease 
    4. Healthy young and Healthy old
    """
    pass
def get_data_participants(name="BDS00001"):
    record = wfdb.rdrecord(name, pn_dir='hbedb/1.0.0/')
    # record = wfdb.rdrecord('BDS00001', pn_dir='hbedb/1.0.0/')
    signal = record.p_signal
    copX, copY = [], []
    for i in range(len(record.p_signal)):
        copX.append(record.p_signal[i][6]) #COPx
        copY.append(record.p_signal[i][7]) #COPy
    return DataFrame({"copX":copX,"copY":copY})
get_bds_info()
print(get_data_participants())
