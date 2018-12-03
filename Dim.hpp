#ifndef DIM_HPP_INCLUDED
#define DIM_HPP_INCLUDED

namespace lab {
    
struct Dim {
    
    int i, j;
    
    bool eq(const Dim& rhs) {
        return i == rhs.i && j == rhs.j;
    }
};

} // namespace lab

#endif // DIM_HPP_INCLUDED