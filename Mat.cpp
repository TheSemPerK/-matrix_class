#include "Mat.hpp"

#include <random>
#include <stdexcept>

namespace lab {
    
Mat::Mat(int m, int n) :
    _m(m), _n(n), _data(0)
{
    
    if (_m <= 0 || _n <= 0) {
        throw std::invalid_argument("Mat(int,int), dim m/n <= 0");
    }
    
    _data = new double*[_m];
    if (_data == 0) {
        throw std::runtime_error("Mat(int,int), double** alloc err");
    } else {
        for (int i = 0; i < _m; ++i) {
            _data[i] = new double[_n];
            if (_data[i] == 0) {
                
                for (--i; i >= 0; --i) {
                    delete [] _data[i];
                }
                delete [] _data;
                
                throw std::runtime_error("Mat(int,int), double* alloc err");
            }
            
            for (int j = 0; j < _n; ++j) {
                _data[i][j] = 0.0;
            }
        }
    }
}

Mat::Mat(const std::vector<double>& vec) :
    _m(vec.size()), _n(1), _data(0)
{
    
    if (_m == 0) {
        return;
    }
    
    _data = new double*[_m];
    if (_data == 0) {
        throw std::runtime_error("Mat(const std::vector<double>&), double** alloc err");
    } else {
        for (int i = 0; i < _m; ++i) {
            _data[i] = new double[1];
            if (_data[i] == 0) {
                
                for (--i; i >= 0; --i) {
                    delete [] _data[i];
                }
                delete [] _data;
                
                throw std::runtime_error("Mat(const std::vector<double>&), double* alloc err");
            }
            
            _data[i][0] = vec[i];
        }
    }
}

Mat::Mat(const double* mem, int m) :
    _m(m), _n(1), _data(0)
{
    
    if (mem == 0) {
        throw std::invalid_argument("Mat(const double*,int), null array");
    }
    
    if (_m <= 0) {
        return;
    }
    
    _data = new double*[_m];
    if (_data == 0) {
        throw std::runtime_error("Mat(const double*,int), double** alloc err");
    } else {
        for (int i = 0; i < _m; ++i) {
            _data[i] = new double[1];
            if (_data[i] == 0) {
                
                for (--i; i >= 0; --i) {
                    delete [] _data[i];
                }
                delete [] _data;
                
                throw std::runtime_error("Mat(const double*,int), double* alloc err");
            }
            
            _data[i][0] = mem[i];
        }
    }
}

Mat::Mat(const Mat& other) :
    _m(other._m), _n(other._n), _data(0)
{
    
    if (other._data == 0 || _m == 0 || _n == 0) {
        return;
    }
    
    _data = new double*[_m];
    if (_data == 0) {
        throw std::runtime_error("Mat(const Mat&), double** alloc err");
    } else {
        for (int i = 0; i < _m; ++i) {
            _data[i] = new double[_n];
            if (_data[i] == 0) {
                
                for (--i; i >= 0; --i) {
                    delete [] _data[i];
                }
                delete [] _data;
                
                throw std::runtime_error("Mat(const Mat&), double* alloc err");
            }
            
            for (int j = 0; j < _n; ++j) {
                _data[i][j] = other._data[i][j];
            }
        }
    }
}

Mat::Mat() :
    _m(0), _n(0), _data(0)
{
    // ctor
}

Mat::~Mat() {
    
    if (_data != 0) {
        for (int i = 0; i < _m; ++i) {
            if (_data[i] != 0) {
                delete [] _data[i];
            }
        }
        delete [] _data;
    }
}

Dim Mat::size() const {
    return Dim{_m, _n};
}

void Mat::set(int i, int j, double v) {
    
    if (i < 0 || i >= _m || j < 0 || j >= _n) {
        throw std::out_of_range("Mat::set, i/j incorrect");
    }
    
    _data[i][j] = v;
}

const double& Mat::get(int i, int j) const {
    
    if (i < 0 || i >= _m || j < 0 || j >= _n) {
        throw std::out_of_range("Mat::get, i/j incorrect");
    }
    
    return _data[i][j];
}

double* Mat::data(int i) {
    
    if (i < 0 || i >= _m) {
        throw std::out_of_range("Mat::data, i incorrect");
    }
    
    return _data[i];
}

double Mat::norm() {
    double sum = 0.0;
    for (int i = 0; i < _m; ++i) {
        for (int j = 0; j < _n; ++j) {
            sum += _data[i][j] * _data[i][j];
        }
    }
    return sqrt(sum);
}

void Mat::uniform_() {
    std::mt19937 m_mersenneTwisterEngine;
    std::uniform_real_distribution<double> urd(0.0, 1.0);
    
    for (int i = 0; i < _m; ++i) {
        for (int j = 0; j < _n; ++j) {
            _data[i][j] = urd(m_mersenneTwisterEngine);
        }
    }
}

Mat& Mat::operator=(const Mat& other) {
    
    if (this != &other) {
        
        if (other._data == 0) {
                
            if (_data != 0) {
                for (int i = 0; i < _m; ++i) {
                    if (_data[i] != 0) {
                        delete [] _data[i];
                    }
                }
                delete [] _data;
            }
                
            _m = 0;
            _n = 0;
            _data = 0;
        
            return *this;
        }
        
        if (_m == other._m && _n == other._n) {
            for (int i = 0; i < _m; ++i) {
                for (int j = 0; j < _m; ++j) {
                    _data[i][j] = other._data[i][j];
                }
            }
        } else {
            
            if (_data != 0) {
                for (int i = 0; i < _m; ++i) {
                    if (_data[i] != 0) {
                        delete [] _data[i];
                    }
                }
                delete [] _data;
            }
        
            _m = other._m;
            _n = other._n;
            
            _data = new double*[_m];
            if (_data == 0) {
                throw std::runtime_error("TMat::operator=, double** alloc err");
            } else {
                for (int i = 0; i < _m; ++i) {
                    _data[i] = new double[_n];
                    if (_data[i] == 0) {
                        
                        for (--i; i >= 0; --i) {
                            delete [] _data[i];
                        }
                        delete [] _data;
                        
                        throw std::runtime_error("TMat::operator=, double* alloc err");
                    }
                
                    for (int j = 0; j < _n; ++j) {
                        _data[i][j] = other._data[i][j];
                    }
                }
            }
        }
    }
    
    return *this;
}

std::ostream& operator<<(std::ostream& os, const Mat& obj) {
    Dim dim(obj.size());
    for (int i = 0; i < dim.i; ++i) {
        for (int j = 0; j < dim.j; ++j) {
            os << obj.get(i, j) << " ";
        }
        os << "\n";
    }
    return os;
}

bool operator==(const Mat& lhs, const Mat& rhs) {
    
    Dim dim(lhs.size());
    if (!dim.eq(rhs.size())) {
        return false;
    }
    
    for (int i = 0; i < dim.i; ++i) {
        for (int j = 0; j < dim.j; ++j) {
            if (lhs.get(i, j) != rhs.get(i, j)) {
                return false;
            }
        }
    }
    
    return true;
}

bool operator!=(const Mat& lhs, const Mat& rhs) {
    
    Dim dim(lhs.size());
    if (!dim.eq(rhs.size())) {
        return true;
    }
    
    for (int i = 0; i < dim.i; ++i) {
        for (int j = 0; j < dim.j; ++j) {
            if (lhs.get(i, j) != rhs.get(i, j)) {
                return true;
            }
        }
    }
    
    return false;
}

Mat operator+(const Mat& lhs, const Mat& rhs) {
    
    Dim dim(lhs.size());
    if (!dim.eq(rhs.size())) {
        throw std::runtime_error("Mat::operator+, can't be done, diff. dimensions");
    }
    
    Mat mat(dim.i, dim.j);
    for (int i = 0; i < dim.i; ++i) {
        for (int j = 0; j < dim.j; ++j) {
            mat.set(i, j, lhs.get(i, j) + rhs.get(i, j));
        }
    }
    return mat;
}

Mat operator-(const Mat& lhs, const Mat& rhs) {
    Dim dim(lhs.size());
    if (!dim.eq(rhs.size())) {
        throw std::runtime_error("Mat::operator-, can't be done, diff. dimensions");
    }
    
    Mat mat(dim.i, dim.j);
    for (int i = 0; i < dim.i; ++i) {
        for (int j = 0; j < dim.j; ++j) {
            mat.set(i, j, lhs.get(i, j) - rhs.get(i, j));
        }
    }
    return mat;
}

Mat operator*(const Mat& lhs, const Mat& rhs) {
    
    Dim ldim(lhs.size());
    Dim rdim(rhs.size());
    
    if (ldim.j != rdim.i) {
        throw std::runtime_error("Mat::operator*, can't be done, diff. dimensions");
    }
    
    Dim dim{ldim.i, rdim.j};
    Mat mat(dim.i, dim.j);
    for (int i = 0; i < dim.i; ++i) {
        for (int j = 0; j < dim.j; ++j) {
            double v = 0.0;
            for (int k = 0; k < ldim.j; ++k) {
                v += lhs.get(i, k) * rhs.get(k, j);
            }
            mat.set(i, j, v);
        }
    }
    return mat;
}

Mat operator*(const Mat& lhs, const double alpha) {
    Dim dim(lhs.size());
    Mat mat(dim.i, dim.j);
    for (int i = 0; i < dim.i; ++i) {
        for (int j = 0; j < dim.j; ++j) {
            mat.set(i, j, lhs.get(i, j) * alpha);
        }
    }
    return mat;
}
    
Mat solve(const Mat& A, const Mat& b) {
    
    Dim dim(A.size());
    
    Mat E(dim.i, dim.j);
    Mat C(A);
    
    for (int i = 0; i < dim.i; ++i) {
        E.set(i, i, 1.0);
    }
    
    for (int i = 0; i < dim.i; ++i) {
            
        double mm = 1.0 / C.get(i,i);
        
        for (int k = 0; k < dim.i; ++k) {
            C.set(i,k, C.get(i,k) * mm);
        }
        
        for (int k = 0; k < dim.i; ++k) {
            E.set(i,k, E.get(i,k) * mm);
        }
        
        for (int j = 0; j < dim.i; ++j) {
            
            if (i == j) {
                continue;
            }
            
            mm = C.get(j,i);
            
            for (int k = 0; k < dim.i; ++k) {
                C.set(j,k, C.get(j,k) - mm * C.get(i,k));
            }
            
            for (int k = 0; k < dim.i; ++k) {
                E.set(j,k, E.get(j,k) - mm * E.get(i,k));
            }
        }
    }
    
    return E * b;
}

} // namespace lab