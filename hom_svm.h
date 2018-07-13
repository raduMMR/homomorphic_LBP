#pragma once

// computes the encrypted equivalent of the simplified Sequential Minimal Optimization
SVM_Solution simplified_SMO(int C, int tol, int max_passes, vector<Data> data);