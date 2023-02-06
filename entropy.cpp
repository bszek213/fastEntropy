#include <iostream>
#include <unordered_map>
#include <vector>
#include <cmath>
// compile with the following command input: "g++ -fPIC -shared -o libEntropy.so entropy.cpp -std=c++17"
double shannonEntropy(const std::vector<double>& data){
    //defining an unordered map that will contain a value along with how many times that value appeared 
    std::unordered_map<double,int> freqMap;
    for (double x: data){
        freqMap[x]++;
    }
    int n = data.size(); //get size of input data
    double entropy = 0; //initialize entropy variable
    for (auto& [x, count] : freqMap){
        double p = count / static_cast<double>(n); //get proba of a certain value
        entropy -= p * log2(p); //calculate shannon entropy
    }
    return entropy;
}

double MSE(const std::vector<double>& data){
    return 0;
}

double approximateEntropy(const std::vector<double>& data, int m, float r){
    return 0;
}

double sampleEntropy(const std::vector<double>& data,int m, float r){
    return 0;
}

double permEntropy(const std::vector<double>& data){
    return 0;
}

extern "C"
{
    double calculate_shannon_entropy(double* data, int data_length)
    {
        std::vector<double> data_vector(data, data + data_length);
        return shannonEntropy(data_vector);
    }
}