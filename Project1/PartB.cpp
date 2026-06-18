#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include "LinearSystem.h"

int main() {
    std::ifstream file("machine.data");
    if (!file.is_open()) {
        std::cout << "Cannot open machine.data" << std::endl;
        return 1;
    }

    std::vector<double> MYCT, MMIN, MMAX, CACH, CHMIN, CHMAX, PRP;

    std::string line;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> cols;
        while (std::getline(ss, token, ','))
            cols.push_back(token);
        if (cols.size() < 10) continue;
        MYCT.push_back(std::stod(cols[2]));
        MMIN.push_back(std::stod(cols[3]));
        MMAX.push_back(std::stod(cols[4]));
        CACH.push_back(std::stod(cols[5]));
        CHMIN.push_back(std::stod(cols[6]));
        CHMAX.push_back(std::stod(cols[7]));
        PRP.push_back(std::stod(cols[8]));
    }
    file.close();

    int total = (int)PRP.size();
    int trainSize = (int)(total * 0.8);
    int testSize  = total - trainSize;

    std::cout << "Total: " << total << ", Train: " << trainSize << ", Test: " << testSize << std::endl;

    Matrix A(trainSize, 6);
    Vector b(trainSize);

    for (int i = 0; i < trainSize; i++) {
        A(i+1, 1) = MYCT[i];
        A(i+1, 2) = MMIN[i];
        A(i+1, 3) = MMAX[i];
        A(i+1, 4) = CACH[i];
        A(i+1, 5) = CHMIN[i];
        A(i+1, 6) = CHMAX[i];
        b(i+1)    = PRP[i];
    }

    OverDeterminedSystem solver(A, b);
    Vector x = solver.SolvePseudoInverse();

    std::cout << "\nModel parameters:" << std::endl;
    std::cout << "  x1 (MYCT)  = " << x(1) << std::endl;
    std::cout << "  x2 (MMIN)  = " << x(2) << std::endl;
    std::cout << "  x3 (MMAX)  = " << x(3) << std::endl;
    std::cout << "  x4 (CACH)  = " << x(4) << std::endl;
    std::cout << "  x5 (CHMIN) = " << x(5) << std::endl;
    std::cout << "  x6 (CHMAX) = " << x(6) << std::endl;

    double rmse = 0.0;
    for (int i = trainSize; i < total; i++) {
        double pred = x(1)*MYCT[i] + x(2)*MMIN[i] + x(3)*MMAX[i]
                    + x(4)*CACH[i] + x(5)*CHMIN[i] + x(6)*CHMAX[i];
        double err = pred - PRP[i];
        rmse += err * err;
    }
    rmse = std::sqrt(rmse / testSize);

    std::cout << "\nTest RMSE: " << rmse << std::endl;

    return 0;
}
