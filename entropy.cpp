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
}
double ApEn(std::vector<double> U, int m, double r) {
    auto maxdist = [](std::vector<double> x_i, std::vector<double> x_j) -> double {
        double max_val = 0.0;
        for (int k = 0; k < x_i.size(); ++k) {
            double diff = std::abs(x_i[k] - x_j[k]);
            if (diff > max_val) max_val = diff;
        }
        return max_val;
    };

    auto phi = [&](int m) -> double {
        int sz = U.size();
        std::vector<std::vector<double>> x(sz - m + 1, std::vector<double>(m, 0.0));
        for (int i = 0; i <= sz - m; ++i) {
            for (int j = 0; j < m; ++j) {
                x[i][j] = U[i+j];
            }
        }

        double C_sum = 0.0;
        for (int i = 0; i < sz - m + 1; ++i) {
            int count = 0;
            for (int j = 0; j < sz - m + 1; ++j) {
                if (maxdist(x[i], x[j]) <= r) {
                    count += 1;
                }
            }
            C_sum += std::log(static_cast<double>(count) / (sz - m + 1));
        }

        return (C_sum / (sz - m + 1));
    };

    int N = U.size();
    return std::abs(phi(m + 1) - phi(m));
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

double sampleEntropy_new(const std::vector<double>& L, int m, double r){
        // Sample entropy
    int N = L.size();
    double B = 0.0;
    double A = 0.0;

    // Split time series and save all templates of length m
    std::vector<std::vector<double>> xmi;
    for (int i = 0; i < N - m; i++) {
        std::vector<double> tmp;
        for (int j = i; j < i + m; j++) {
            tmp.push_back(L[j]);
        }
        xmi.push_back(tmp);
    }

    std::vector<std::vector<double>> xmj;
    for (int i = 0; i < N - m + 1; i++) {
        std::vector<double> tmp;
        for (int j = i; j < i + m; j++) {
            tmp.push_back(L[j]);
        }
        xmj.push_back(tmp);
    }

    // Save all matches minus the self-match, compute B
    for (auto xmii : xmi) {
        int count = 0;
        for (auto xmjji : xmj) {
            bool match = true;
            for (int i = 0; i < m; i++) {
                if (std::abs(xmii[i] - xmjji[i]) > r) {
                    match = false;
                    break;
                }
            }
            if (match) {
                count++;
            }
        }
        B += count - 1;
    }

    // Similar for computing A
    m += 1;
    std::vector<std::vector<double>> xm;
    for (int i = 0; i < N - m + 1; i++) {
        std::vector<double> tmp;
        for (int j = i; j < i + m; j++) {
            tmp.push_back(L[j]);
        }
        xm.push_back(tmp);
    }

    for (auto xmi : xm) {
        int count = 0;
        for (auto xmj : xm) {
            bool match = true;
            for (int i = 0; i < m; i++) {
                if (std::abs(xmi[i] - xmj[i]) > r) {
                    match = false;
                    break;
                }
            }
            if (match) {
                count++;
            }
        }
        A += count - 1;
    }

    // Return SampEn
    return -std::log(A / B);
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