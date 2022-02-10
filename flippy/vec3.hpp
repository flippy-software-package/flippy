#ifndef FLIPPY_VEC3_HPP
#define FLIPPY_VEC3_HPP

#include <ostream>
#include <iostream>
#include <cmath>
#include <concepts>

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
template<std::floating_point Real>
class vec3
{
public:
    Real x, y, z;

    void add(vec3<Real> const& v)
    {
        x += v.x;
        y += v.y;
        z += v.z;
    }

    void subtract(vec3<Real> const& v)
    {
        x -= v.x;
        y -= v.y;
        z -= v.z;
    }

    void scale(Real s)
    {
        x = s*x;
        y = s*y;
        z = s*z;
    }

    Real dot(vec3<Real> const& v) const
    {
        Real res = x*v.x + y*v.y + z*v.z;
        return res;
    }

    [[nodiscard]] constexpr std::size_t size() const { return 3; }

    static inline vec3<Real> cross(vec3<Real> const& a, vec3<Real> const& b)
    {
        vec3<Real> res;
        res.x = a.y*b.z - a.z*b.y;
        res.y = a.z*b.x - a.x*b.z;
        res.z = a.x*b.y - a.y*b.x;
        return res;
    }

    vec3<Real> cross(vec3<Real> const& other) const { return cross(*this, other); }

    Real norm() const { return std::sqrt(this->dot(*this)); }

    Real norm_square() const { return this->dot(*this); }

    void normalize(){
        /**
         * normalize vector in place.
         */
        *this= *this/this->norm();
    }

    friend std::ostream& operator<<(std::ostream& os, const vec3<Real>& obj)
    {
        os << "{" << obj.x << ',' << obj.y << ',' << obj.z << '}';
        return os;
    }

    bool operator==(vec3<Real> const& other) const =default;


    // mathematical operations
    friend vec3<Real> operator+(vec3<Real> lhs, vec3<Real> const& rhs)
    {
        lhs+=rhs;
        return lhs;
    }

    friend void operator+=(vec3<Real>& lhs, vec3<Real> const& rhs)
    {
        lhs.add(rhs);
    }

    friend vec3<Real> operator-(vec3<Real> lhs, vec3<Real> const& rhs)
    {
        lhs-=rhs;
        return lhs;
    }

    friend void operator-=(vec3<Real>& lhs, vec3<Real> const& rhs)
    {
        lhs.subtract(rhs);
    }

    friend vec3<Real> operator*(Real const& lhs, vec3<Real> rhs)
    {
        rhs.scale(lhs);
        return rhs;
    }

    friend vec3<Real> operator*(vec3<Real> lhs, Real const& rhs)
    {
        lhs.scale(rhs);
        return lhs;
    }

    friend void operator/=(vec3<Real>& lhs, Real const& rhs){
        lhs.scale((Real)1/rhs);
    }
    friend vec3<Real> operator/(vec3<Real> lhs, Real const& rhs)
    {
        lhs/=rhs;
        return lhs;
    }

    template<typename Index>
    requires std::is_integral_v<Index>
    Real& operator[](Index idx)
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
    const Real& operator[](Index idx) const
    {
        switch (idx) {
            case 0:return x;
            case 1:return y;
            case 2:return z;
            default:std::cerr << idx << "is out of range for as vec3 index";
                exit(12);
        }
    }

    friend vec3<Real> operator-(vec3<Real> v)
    {
        v.x = -v.x;
        v.y = -v.y;
        v.z = -v.z;
        return v;
    }

    friend vec3<Real>& operator-(vec3<Real>&& v)
    {
        v.x = -v.x;
        v.y = -v.y;
        v.z = -v.z;
        return v;
    }

};
}

#endif //FLIPPY_VEC3_HPP
