#include <iostream>
#include <unordered_map>
#include <vector>
#include <cmath>
// compile with the following command input: "g++ -fPIC -shared -o libEntropy.so entropy.cpp -std=c++17"
double shannonEntropy(const std::vector<double>& data){
    //Round all data to the tenths place - change this depending on input data
    std::vector<double> dataRound;
    for (double x : data) {
        double rounded = std::round(x * 10.0) / 10.0;
        dataRound.push_back(rounded);
    }
    //defining an unordered map that will contain a value along with how many times that value appeared 
    std::unordered_map<double,int> freqMap;
    for (double x: dataRound){
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

double calculatMean(const std::vector<double>& data) {
    double sum = 0;
    for (double x : data) {
        sum += x;
    }
    return sum / data.size();
}

double StandardDeviation(const std::vector<double>& data) {
    double mean = calculatMean(data);
    double sum_of_squares = 0;
    for (double x : data) {
        double deviation = x - mean;
        sum_of_squares += deviation * deviation;
    }
    return std::sqrt(sum_of_squares / (data.size() - 1));
}
double SampleEntropy_old(const std::vector<double>& data, int m, double r)
{
    int N = data.size();
    int Cm = 0, Cm1 = 0;
    double err = 0.0, sum = 0.0;
    double sd = StandardDeviation(data);
    err = sd * r;
    
    for (unsigned int i = 0; i < N - (m + 1) + 1; i++) {
        for (unsigned int j = i + 1; j < N - (m + 1) + 1; j++) {      
        bool eq = true;
        //m - length series
        for (unsigned int k = 0; k < m; k++) {
            if (std::abs(data[i+k] - data[j+k]) > err) {
            eq = false;
            break;
            }
        }
        if (eq) Cm++;
        
        //m+1 - length series
        int k = m;
        if (eq && std::abs(data[i+k] - data[j+k]) <= err)
            Cm1++;
        }
    }
    
    if (Cm > 0 && Cm1 > 0)
        return std::log((double)Cm / (double)Cm1);
    else
        return 0.0; 
  
}

double sampleEntropy_new(const std::vector<double>& data, int m, double r){
    //multiply r by the std
    double sd = StandardDeviation(data);
    double err = sd * r;
    int n = data.size();
    std::vector<double> last(m), curr(m);
    int count1 = 0, count2 = 0;

    for (int i = 0; i < n - m; i++) {
        for (int j = 0; j < m; j++) {
            last[j] = data[i + j];
            curr[j] = data[i + j + 1];
        }

        int matches1 = 0, matches2 = 0;
        for (int j = 0; j < m; j++) {
            for (int k = j + 1; k < m; k++) {
                if (std::abs(last[j] - last[k]) < err) {
                    matches1++;
                }
                if (std::abs(curr[j] - curr[k]) < err) {
                    matches2++;
                }
            }
        }
        count1 += matches1;
        count2 += matches2;
    }

    double e1 = std::log(static_cast<double>(count1) / (n - m));
    double e2 = std::log(static_cast<double>(count2) / (n - m - 1));
    return e1 - e2;
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
    double calculate_sample_entropy(double* data, int data_length)
    {
        std::vector<double> data_vector(data, data + data_length);
        return SampleEntropy_old(data_vector,4,0.25);
    }
}