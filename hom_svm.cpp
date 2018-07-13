#include "hom_svm.h"
#include <iostream>
#include <vector>
using namespace std;

struct Data{
    vector<double> x;
    int y;
};

struct SVM_Solution{
    vector<double> alphas;
    double b;

    SVM_Solution(int m){
        alphas = vector<double>(m, 0);
        b = 0;
    }
};

double dot_product(const vector<double> &vec1, const vector<double> &vec2){
    assert(vec1.size() == vec2.size());
    double dt = 0;
    for(int i=0; i<vec1.size(); i++){
        dt += vec1[i] * vec2[i];
    }
    return dt;
}

// f(x) = (i=1..m)SUM(alpha_i*y_i*<x_i, x>)+b
double f(const vector<Data> &data, const vector<double> &x, const vector<double> &alphas, const double b){
    double sum = 0;
    for(int i=0; i<data.size(); i++){
        sum += alphas[i] * data[i].y_i * dot_product(data[i].x_i, x) + b;
    }
    return sum;
}

double max(const double a, const double b){
    return a >= b ? a : b;
}

double min(const double a, const double b){
    return a < b ? a : b;
}

/**
* C - regularization parameter
* tol - numerical tolerance !!! not necessary an integer
* max_passes - max # of times to iterate over alpha's without changing
* (x_i, y_i) - training data
*/
SVM_Solution simplified_SMO(int C, int tol, int max_passes, vector<Data> data){
    int m = data.size();
    vector<double> alphas(m, 0);
    double b = 0;

    int passes = 0;
    while(passes < max_passes){
        int num_changed_alphas = 0;
        for(int i=0; i<m; i++){
            double E_i = f(data, data[i].x_i, alphas, b) - data[i].y_i;
            if((data[i].y*E_i<(-1)*tol && alphas[i]<C) || (data[i].y*E_i > tol && alphas[i] > 0)){
                int j = rand() % m;
                double E_j = f(data, data[j].x_i, alphas, b) - data[j].y_i;
                double old_alpha_i = alphas[i];
                double old_alpha_j = alphas[j];
                double L, H;
                if(data[i].y != data[j].y){
                    L = max(0, alphas[j]-alphas[i]);
                    H = min(C, C+alphas[j]-alphas[i]);
                }
                else {
                    L = max(0, alphas[i]+alphas[j]-C);
                    H = min(C, alphas[i]+alphas[j]);
                }
                if(L == H){
                    continue;
                }

                double eta = 2*dot_product(data[i].x, data[j].x) - dot_product(data[i].x, data[i].x) -dot_product(data[j].x, data[j].x);
                if(eta >= 0){
                    continue;
                }

                alphas[j] = alphas[j] - data[j].y*(E_i - E_j)/eta;
                if(alphas[j] > H){
                    alphas[j] = H;
                }
                else{
                    if(alphas[j] < L){
                        alphas[j] = L;
                    }
                }
                if(abs(alphas[j] - old_alpha_j) < 0,00001){
                    continue;
                }

                alphas[i] = alphas[i] + data[i].y*data[i].y(old_alpha_j-alphas[j]);

                double b1 = b - E_i - data[i].y*(alphas[i]-old_alpha_i)*dot_product(data[i].x, data[i].x)-data[j].y*(alphas[j]-old_alpha_j)*dot_product(data[i].x, data[j].x);
                double b1 = b - E_j - data[i].y*(alphas[i]-old_alpha_i)*dot_product(data[i].x, data[j].x)-data[j].y*(alphas[j]-old_alpha_j)*dot_product(data[i].x, data[j].x);
                if(alphas[i] < C && alphas[i] > 0){
                    b = b1;
                }
                if(alphas[j] > 0 && alphas[j]<C){
                    b = b2;
                } else {
                    b = (b1+b2)/2;
                }

                num_changed_alphas++;
            }
        }
        if(num_changed_alphas == 0){
            passes++;
        }
        else{
            passes = 0;
        }
    }

    SVM_Solution solution;
    solution.alphas = alphas;
    solution.b = b;
    return solution;
}