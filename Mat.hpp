#ifndef MAT_HPP
#define MAT_HPP

#include <iostream>
#include <vector>

#include "Dim.hpp"

namespace lab {

class Mat {
    
int _m, _n;
double** _data;
    
public:
    
    Mat(int m, int n);
    Mat(const std::vector<double>& vec);
    Mat(const double* mem, int m);
    Mat(const Mat& other);
    Mat();
    
    ~Mat();
    
    Dim size() const;
    void set(int i, int j, double v);
    const double& get(int i, int j) const;
    double* data(int i);
    double norm();
    void uniform_();
    
    Mat& operator=(const Mat& other);

};

std::ostream& operator<<(std::ostream& os, const Mat& obj);

bool operator==(const Mat& lhs, const Mat& rhs);
bool operator!=(const Mat& lhs, const Mat& rhs);

Mat operator+(const Mat& lhs, const Mat& rhs);
Mat operator-(const Mat& lhs, const Mat& rhs);
Mat operator*(const Mat& lhs, const Mat& rhs);
Mat operator*(const Mat& lhs, const double alpha);
    
Mat solve(const Mat& A, const Mat& b);

} // namespace lab

#endif // MAT_HPP