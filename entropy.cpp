#include <iostream>
#include <unordered_map>
#include <vector>
#include <cmath>
// compile with the following command input: "g++ -fPIC -shared -o libEntropy.so entropy.cpp -std=c++17"
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

double MSE(std::vector<double>& x, int m, double r, int tau, double& e, int& A, int& B){
    // Coarse signal
    std::vector<double> y;
    for (int i = 0; i < x.size(); i += tau) {
        double sum = 0.0;
        for (int j = i; j < std::min(i + tau, (int)x.size()); j++) {
            sum += x[j];
        }
        y.push_back(sum / (double)tau);
    }

    // (m+1)-element sequences
    std::vector<std::vector<double>> X;
    for (int i = 0; i < y.size() - m - 1; i++) {
        std::vector<double> tmp;
        for (int j = i; j < i + m + 1; j++) {
            tmp.push_back(y[j]);
        }
        X.push_back(tmp);
    }

    // Matching (m+1)-element sequences
    A = 0;
    double stdev = StandardDeviation(x);
    for (int i = 0; i < X.size(); i++) {
        for (int j = i + 1; j < X.size(); j++) {
            bool match = true;
            for (int k = 0; k < m + 1; k++) {
                if (std::abs(X[i][k] - X[j][k]) > r * stdev) {
                    match = false;
                    break;
                }
            }
            if (match) {
                A++;
            }
        }
    }

    // Matching m-element sequences
    X.erase(X.begin() + m + 1, X.end());
    B = 0;
    for (int i = 0; i < X.size(); i++) {
        for (int j = i + 1; j < X.size(); j++) {
            bool match = true;
            for (int k = 0; k < m; k++) {
                if (std::abs(X[i][k] - X[j][k]) > r * stdev) {
                    match = false;
                    break;
                }
            }
            if (match) {
                B++;
            }
        }
    }

    // Take log
    if (A == 0 || B == 0) {
        e = NAN;
    } else {
        e = std::log(B / (double)A);
    }
    return e;
}
// Function to calculate the approximate entropy (ApEn)
double ApEn(const std::vector<double>& data, int m, double r) {
    int N = data.size();
    int count = 0;

    for (int i = 0; i <= N - m; ++i) {
        for (int j = i + 1; j <= N - m; ++j) {
            bool match = true;
            for (int k = 0; k < m; ++k) {
                if (std::fabs(data[i + k] - data[j + k]) > r) {
                    match = false;
                    break;
                }
            }
            if (match) {
                ++count;
            }
        }
    }

    double ApEn = -std::log(static_cast<double>(count) / (N - m + 1));
    return ApEn;
}

double SampleEntropy_old(const std::vector<double>& data, int m, double r)
// got the equation from this source: https://www.codeproject.com/Articles/27030/Approximate-and-Sample-Entropies-Complexity-Metric
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

// double sampleEntropy(const std::vector<double> data, int m, double r) 
// {
//     int n = data.size();
//     int count1 = 0, count2 = 0;
//     double maxdist1, maxdist2;

//     for (int i = 0; i < n - m + 1; i++) {
//         std::vector<double> temp1(m), temp2(m);

//         for (int j = 0; j < m; j++) {
//             temp1[j] = data[i + j];
//             temp2[j] = data[i + j + 1];
//         }

//         maxdist1 = maxdist2 = 0;

//         for (int j = 0; j < n - m + 1; j++) {
//             bool match1 = true, match2 = true;

//             for (int k = 0; k < m; k++) {
//                 if (abs(data[j + k] - temp1[k]) > r) {
//                     match1 = false;
//                 }
//                 if (abs(data[j + k] - temp2[k]) > r) {
//                     match2 = false;
//                 }
//             }

//             if (match1) {
//                 count1++;
//                 if (abs(data[j + m] - temp1[m - 1]) > maxdist1) {
//                     maxdist1 = abs(data[j + m] - temp1[m - 1]);
//                 }
//             }

//             if (match2) {
//                 count2++;
//                 if (abs(data[j + m] - temp2[m - 1]) > maxdist2) {
//                     maxdist2 = abs(data[j + m] - temp2[m - 1]);
//                 }
//             }
//         }
//     }

//     double sample_entropy = -log((count1 * 1.0) / (count2 * 1.0));
//     return sample_entropy;
// }

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
        return SampleEntropy_old(data_vector,2,0.2);
    }
    double calculate_approximate_entropy(double* data, int data_length)
    {
        std::vector<double> data_vector(data, data + data_length);
        return ApEn(data_vector,2,0.2);
    }
}