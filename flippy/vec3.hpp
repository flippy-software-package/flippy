#ifndef EXPERIMENT_H_UTILITIES_VEC3_H_
#define EXPERIMENT_H_UTILITIES_VEC3_H_

#include <array>
#include <initializer_list>
#include <ostream>
#include <iostream>
#include <cmath>

namespace fp{
template<typename Type>
class vec3 : public std::array<Type, 3> {
public:

    vec3() = default;

    vec3(Type x, Type y, Type z)
    {
        (*this)[0] = x;
        (*this)[1] = y;
        (*this)[2] = z;
    };

    vec3(std::initializer_list<Type> const& lst)
    {
        (*this)[0] = *lst.begin();
        (*this)[1] = *(lst.begin() + 1);
        (*this)[2] = *(lst.begin() + 2);
    };

    void add(vec3<Type> const& v)
    {
        for (std::size_t i = 0; i<this->size(); ++i) { (*this)[i] += v[i]; }
    }

    void subtract(vec3<Type> const& v)
    {
        for (std::size_t i = 0; i<this->size(); ++i) { (*this)[i] -= v[i]; }
    }

    void scale(Type const& s)
    {
        for (std::size_t i = 0; i<this->size(); ++i) { (*this)[i] = s*((*this)[i]); }
    }

    Type dot(vec3<Type> const& v) const
    {
        Type res = 0.;
        for (std::size_t i = 0; i<this->size(); ++i) { res += (*this)[i]*v[i]; }
        return res;
    }

    static inline vec3<Type> cross(vec3<Type> const& a, vec3<Type> const& b)
    {
        vec3<Type> res;
        res[0] = a[1]*b[2] - a[2]*b[1];
        res[1] = a[2]*b[0] - a[0]*b[2];
        res[2] = a[0]*b[1] - a[1]*b[0];
        return res;
    }

    vec3<Type> cross(vec3<Type> const& other) const { return cross(*this, other); }

    Type norm() const { return sqrt(this->dot(*this)); }

    Type norm_square() const { return this->dot(*this); }

    friend std::ostream& operator<<(std::ostream& os, const vec3<Type>& obj)
    {
        bool cond(obj.size()<10);
        int display_size(cond ? obj.size()/2 : 4);
        os << "{ ";
        for (int i = 0; i<display_size; ++i) { os << obj[i] << " "; }
        if (not cond) { os << " ... "; }
        for (std::size_t i = display_size; i<obj.size(); ++i) { os << obj[i] << " "; }

        os << "}";

        return os;
    }

   bool operator==(vec3<Type> const& other) const =default;

    // mathematical operations
    friend vec3<Type> operator+(vec3<Type> lhs, vec3<Type> const& rhs)
    {
        lhs+=rhs;
        return lhs;
    }

    friend void operator+=(vec3<Type>& lhs, vec3<Type> const& rhs)
    {
        lhs.add(rhs);
    }

    friend vec3<Type> operator-(vec3<Type> lhs, vec3<Type> const& rhs)
    {
        lhs-=rhs;
        return lhs;
    }

    friend void operator-=(vec3<Type>& lhs, vec3<Type> const& rhs)
    {
        lhs.subtract(rhs);
    }

    friend vec3<Type> operator*(Type const& lhs, vec3<Type> rhs)
    {
        rhs.scale(lhs);
        return rhs;
    }

    friend vec3<Type> operator*(vec3<Type> lhs, Type const& rhs)
    {
        lhs.scale(rhs);
        return lhs;
    }

    friend void operator/=(vec3<Type>& lhs, Type const& rhs){
        lhs.scale((Type)1/rhs);
    }
    friend vec3<Type> operator/(vec3<Type> lhs, Type const& rhs)
    {
        lhs/=rhs;
        return lhs;
    }

    friend vec3<Type> operator-(vec3<Type> v)
    {
        for (std::size_t i = 0; i<v.size(); ++i) { v[i] = -v[i]; }
        return v;
    }

};
}

#endif //EXPERIMENT_H_UTILITIES_VEC3_H_
