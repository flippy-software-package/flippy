#ifndef FLIPPY_VEC3_HPP
#define FLIPPY_VEC3_HPP

#include <ostream>
#include <iostream>
#include <cmath>

namespace fp {
template<typename Type>
class vec3
{
public:
    Type x, y, z;

    void add(vec3<Type> const& v)
    {
        x += v.x;
        y += v.y;
        z += v.z;
    }

    void subtract(vec3<Type> const& v)
    {
        x -= v.x;
        y -= v.y;
        z -= v.z;
    }

    void scale(Type const& s)
    {
        x = s*x;
        y = s*y;
        z = s*z;
    }

    Type dot(vec3<Type> const& v) const
    {
        return x*v.x + y*v.y + z*v.z;
    }

    [[nodiscard]] constexpr std::size_t size() const { return 3; }

    static inline vec3<Type> cross(vec3<Type> const& a, vec3<Type> const& b)
    {
        vec3<Type> res;
        res.x = a.y*b.z - a.z*b.y;
        res.y = a.z*b.x - a.x*b.z;
        res.z = a.x*b.y - a.y*b.x;
        return res;
    }

    vec3<Type> cross(vec3<Type> const& other) const { return cross(*this, other); }

    Type norm() const { return std::sqrt(this->dot(*this)); }

    Type norm_square() const { return this->dot(*this); }

    friend std::ostream& operator<<(std::ostream& os, const vec3<Type>& obj)
    {
        os << "{" << obj.x << ',' << obj.y << ',' << obj.z << '}';
        return os;
    }

    friend bool operator==(vec3<Type> const& lhs, vec3<Type> const& rhs)
    {
        return (lhs.x==rhs.x) && (lhs.y==rhs.y) && (lhs.z==rhs.z);
    }

    // mathematical operations
    friend vec3<Type> operator+(vec3<Type> lhs, vec3<Type> const& rhs)
    {
        lhs.add(rhs);
        return lhs;
    }

    friend void operator+=(vec3<Type>& lhs, vec3<Type> const& rhs)
    {
        lhs.add(rhs);
    }

    friend vec3<Type> operator-(vec3<Type> lhs, vec3<Type> const& rhs)
    {
        lhs.subtract(rhs);
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

    friend vec3<Type> operator/(vec3<Type> lhs, Type rhs)
    {
        lhs.scale((Type) 1/rhs);
        return lhs;
    }

    Type& operator[](int idx)
    {
        switch (idx) {
            case 0:return x;
            case 1:return y;
            case 2:return z;
            default:std::cerr << idx << "is out of range for as vec3 index";
                exit(12);
        }
    }

    const Type& operator[](int idx) const
    {
        switch (idx) {
            case 0:return x;
            case 1:return y;
            case 2:return z;
            default:std::cerr << idx << "is out of range for as vec3 index";
                exit(12);
        }
    }

    friend vec3<Type> operator-(vec3<Type> v)
    {
        v.x = -v.x;
        v.y = -v.y;
        v.z = -v.z;
        return v;
    }

    friend vec3<Type>& operator-(vec3<Type>&& v)
    {
        v.x = -v.x;
        v.y = -v.y;
        v.z = -v.z;
        return v;
    }

};
}

#endif //FLIPPY_VEC3_HPP
