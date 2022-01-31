#ifndef FLIPPY_VEC3_HPP
#define FLIPPY_VEC3_HPP

#include <ostream>
#include <iostream>
#include <cmath>

namespace fp{

/**
 * Own implementation of a 3D vector.
 * To keep the external dependencies low, flippy implements it's own 3D vector class with basic functionality like dot product and cross product
 *
 * Example:
 * ```c++
 * fp::vec3<double> v1{1,0,0};
 * fp::vec3<double> v2{0,0,1};
 *
 *  assert(v1.dot(v2)==0);
 *  assert(v1.cross(v2).norm()==1);
 *  assert(((v1-v2)==fp::vec3<double>{1.,0.,-1.}));
 * ```
 */
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

    void scale(Type s)
    {
        x = s*x;
        y = s*y;
        z = s*z;
    }

    Type dot(vec3<Type> const& v) const
    {
        Type res = x*v.x + y*v.y + z*v.z;
        return res;
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

    void normalize(){
        /**
         * normalize vector in place.
         */
        *this= *this/this->norm();
    }

    friend std::ostream& operator<<(std::ostream& os, const vec3<Type>& obj)
    {
        os << "{" << obj.x << ',' << obj.y << ',' << obj.z << '}';
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

    template<typename Index>
    requires std::is_integral_v<Index>
    Type& operator[](Index idx)
    {
        switch (idx) {
            case 0:return x;
            case 1:return y;
            case 2:return z;
            default:std::cerr << idx << "is out of range for as vec3 index";
            exit(12);
        }
    }

    template<typename Index>
    requires std::is_integral_v<Index>
    const Type& operator[](Index idx) const
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
