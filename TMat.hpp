#ifndef TMAT_HPP
#define TMAT_HPP

#include <iostream>
#include <vector>

#include <cstdlib>
#include <random>
#include <stdexcept>

#include "Dim.hpp"

namespace lab {

template <typename T>
class TMat {
    
int _m, _n;
T** _data;
    
public:
    
    TMat(int m, int n) :
        _m(m), _n(n), _data(0)
    {
        
        if (_m <= 0 || _n <= 0) {
            throw std::invalid_argument("TMat(int,int), dim i/j <= 0");
        }
        
        _data = (T**) malloc(sizeof(T**) * _m);
        if (_data == 0) {
            throw std::runtime_error("TMat(int,int), T** alloc err");
        } else {
            for (int i = 0; i < _m; ++i) {
                _data[i] = (T*) malloc(sizeof(T*) * _n);
                if (_data[i] == 0) {
                    
                    for (--i; i >= 0; --i) {
                        free(_data[i]);
                    }
                    free(_data);
                    
                    throw std::runtime_error("TMat(int,int), T* alloc err");
                }
                
                for (int j = 0; j < _n; ++j) {
                    _data[i][j] = T(0);
                }
            }
        }
    }

    TMat(const std::vector<T>& vec) :
        _m(vec.size()), _n(1), _data(0)
    {
        
        if (_m == 0) {
            return;
        }
        
        _data = (T**) malloc(sizeof(T**) * _m);
        if (_data == 0) {
            throw std::runtime_error("TMat(const std::vector<T>&), T** alloc err");
        } else {
            for (int i = 0; i < _m; ++i) {
                _data[i] = (T*) malloc(sizeof(T*));
                if (_data[i] == 0) {
                    
                    for (--i; i >= 0; --i) {
                        free(_data[i]);
                    }
                    free(_data);
                    
                    throw std::runtime_error("TMat(const std::vector<T>&), T* alloc err");
                }
                
                _data[i][0] = vec[i];
            }
        }
    }

    TMat(const T* mem, int m) :
        _m(m), _n(1), _data(0)
    {
        
        if (mem == 0) {
            throw std::invalid_argument("TMat(const T*,int), null array");
        }
        
        if (_m == 0) {
            return;
        }
        
        _data = (T**) malloc(sizeof(T**) * _m);
        if (_data == 0) {
            throw std::runtime_error("TMat(const T*,int), T** alloc err");
        } else {
            for (int i = 0; i < _m; ++i) {
                _data[i] = (T*) malloc(sizeof(T*));
                if (_data[i] == 0) {
                    
                    for (--i; i >= 0; --i) {
                        free(_data[i]);
                    }
                    free(_data);
                    
                    throw std::runtime_error("TMat(const T*,int), T* alloc err");
                }
                
                _data[i][0] = mem[i];
            }
        }
    }

    TMat(const TMat<T>& other) :
        _m(other._m), _n(other._n), _data(0)
    {
        
        if (other._data == 0 || _m == 0 || _n == 0) {
            return;
        }
        
        _data = (T**) malloc(sizeof(T**) * _m);
        if (_data == 0) {
            throw std::runtime_error("TMat(const TMat&), T** alloc err");
        } else {
            for (int i = 0; i < _m; ++i) {
                _data[i] = (T*) malloc(sizeof(T*) * _n);
                if (_data[i] == 0) {
                    
                    for (--i; i >= 0; --i) {
                        free(_data[i]);
                    }
                    free(_data);
                    
                    throw std::runtime_error("TMat(const TMat&), T* alloc err");
                }
                
                for (int j = 0; j < _n; ++j) {
                    _data[i][j] = other._data[i][j];
                }
            }
        }
    }

    TMat() :
        _m(0), _n(0), _data(0)
    {

    }

    ~TMat() {
        
        if (_data != 0) {
            for (int i = 0; i < _m; ++i) {
                if (_data[i] != 0) {
                    free(_data[i]);
                }
            }
            free(_data);
        }
    }
    
    Dim size() const {
        return Dim{_m, _n};
    }

    void set(int i, int j, T v) {
        
        if (i < 0 || i >= _m || j < 0 || j >= _n) {
            throw std::out_of_range("TMat::set, i/j incorrect");
        }
        
        _data[i][j] = v;
    }

    const T& get(int i, int j) const {
        
        if (i < 0 || i >= _m || j < 0 || j >= _n) {
            throw std::out_of_range("TMat::get, i/j incorrect");
        }
        
        return _data[i][j];
    }

    T* data(int i) {
        
        if (i < 0 || i >= _m) {
            throw std::out_of_range("TMat::data, i incorrect");
        }
        
        return _data[i];
    }

    double norm() {
        double sum = 0.0;
        for (int i = 0; i < _m; ++i) {
            for (int j = 0; j < _n; ++j) {
                sum += (double) _data[i][j] * _data[i][j];
            }
        }
        return sqrt(sum);
    }

    void uniform_() {
        std::mt19937 m_mersenneTwisterEngine;
        std::uniform_real_distribution<double> urd(0.0, 1.0);
        
        for (int i = 0; i < _m; ++i) {
            for (int j = 0; j < _n; ++j) {
                _data[i][j] = urd(m_mersenneTwisterEngine);
            }
        }
    }
    
    TMat<T>& operator=(const TMat<T>& other) {
        
        if (this != &other) {
            
            if (other._data == 0) {
                    
                if (_data != 0) {
                    for (int i = 0; i < _m; ++i) {
                        if (_data[i] != 0) {
                            free(_data[i]);
                        }
                    }
                    free(_data);
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
                            free(_data[i]);
                        }
                    }
                    free(_data);
                }
            
                _m = other._m;
                _n = other._n;
                
                _data = (T**) malloc(sizeof(T**) * _m);
                if (_data == 0) {
                    throw std::runtime_error("TMat::operator=, T** alloc err");
                } else {
                    for (int i = 0; i < _m; ++i) {
                        _data[i] = (T*) malloc(sizeof(T*) * _n);
                        if (_data[i] == 0) {
                            
                            for (--i; i >= 0; --i) {
                                free(_data[i]);
                            }
                            free(_data);
                            
                            throw std::runtime_error("TMat::operator=, T* alloc err");
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

};

template <typename T>
std::ostream& operator<<(std::ostream& os, const TMat<T>& obj) {
    Dim dim(obj.size());
    for (int i = 0; i < dim.i; ++i) {
        for (int j = 0; j < dim.j; ++j) {
            os << obj.get(i, j) << " ";
        }
        os << "\n";
    }
    return os;
}

template <typename T>
bool operator==(const TMat<T>& lhs, const TMat<T>& rhs) {
    
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

template <typename T>
bool operator!=(const TMat<T>& lhs, const TMat<T>& rhs) {
    
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

template <typename T>
TMat<T> operator+(const TMat<T>& lhs, const TMat<T>& rhs) {
    
    Dim dim(lhs.size());
    if (!dim.eq(rhs.size())) {
        throw std::runtime_error("TMat::operator+, can't be done, diff. dimensions");
    }
    
    TMat<T> mat(dim.i, dim.j);
    for (int i = 0; i < dim.i; ++i) {
        for (int j = 0; j < dim.j; ++j) {
            mat.set(i, j, lhs.get(i, j) + rhs.get(i, j));
        }
    }
    return mat;
}

template <typename T>
TMat<T> operator-(const TMat<T>& lhs, const TMat<T>& rhs) {
    Dim dim(lhs.size());
    if (!dim.eq(rhs.size())) {
        throw std::runtime_error("TMat::operator-, can't be done, diff. dimensions");
    }
    
    TMat<T> mat(dim.i, dim.j);
    for (int i = 0; i < dim.i; ++i) {
        for (int j = 0; j < dim.j; ++j) {
            mat.set(i, j, lhs.get(i, j) - rhs.get(i, j));
        }
    }
    return mat;
}

template <typename T>
TMat<T> operator*(const TMat<T>& lhs, const TMat<T>& rhs) {
    
    Dim ldim(lhs.size());
    Dim rdim(rhs.size());
    
    if (ldim.j != rdim.i) {
        throw std::runtime_error("TMat::operator*, can't be done, diff. dimensions");
    }
    
    Dim dim{ldim.i, rdim.j};
    TMat<T> mat(dim.i, dim.j);
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

template <typename T>
TMat<T> operator*(const TMat<T>& lhs, const double alpha) {
    Dim dim(lhs.size());
    TMat<T> mat(dim.i, dim.j);
    for (int i = 0; i < dim.i; ++i) {
        for (int j = 0; j < dim.j; ++j) {
            mat.set(i, j, lhs.get(i, j) * alpha);
        }
    }
    return mat;
}
        
template <typename T>
TMat<T> solve(const TMat<T>& A, const TMat<T>& b) {
    
    // A(m,n) * x(n,1) = b(m,1)
    // x = A^-1 * b
    
    if (std::is_floating_point<T>::value) {
        Dim dim(A.size());
        TMat<T> E(dim.i, dim.j);
        TMat<T> C(A);
        
        for (int i = 0; i < dim.i; ++i) {
            E.set(i, i, T(1));
        }
        
        for (int i = 0; i < dim.i; ++i) {
                
            T mm = 1.0 / C.get(i,i);
            
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
    
    return TMat<T>();
}

} // namespace lab

#endif // TMAT_HPP