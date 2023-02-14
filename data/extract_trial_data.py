import wfdb
import matplotlib.pyplot as plt
from pandas import DataFrame
def get_data_participants(name="BDS00001"):
    record = wfdb.rdrecord(name, pn_dir='hbedb/1.0.0/')
    # record = wfdb.rdrecord('BDS00001', pn_dir='hbedb/1.0.0/')
    signal = record.p_signal
    copX, copY = [], []
    for i in range(len(record.p_signal)):
        copX.append(record.p_signal[i][6]) #COPx
        copY.append(record.p_signal[i][7]) #COPy
    return DataFrame({"copX":copX,"copY":copY})
print(get_data_participants())
