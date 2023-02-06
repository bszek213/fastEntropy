#Python code that calls c++ shared object files
import ctypes
from pandas import read_csv
# import numpy as np
# Load the shared library
shannon_entropy_lib = ctypes.CDLL('./libEntropy.so')
# Declare the C function
shannon_entropy_func = shannon_entropy_lib.calculate_shannon_entropy 
# shannon_entropy_func.argtypes = [np.ctypeslib.ndpointer(dtype=np.double), ctypes.c_int]
shannon_entropy_func.restype = ctypes.c_double
#Load in data
df = read_csv('test_cop.csv')
# Convert the input to C array
data = df['COPx[cm]'].to_list()
data_array = (ctypes.c_double * len(data))(*data)
# Call the C function
result = shannon_entropy_func(data_array, len(data))
print(result)