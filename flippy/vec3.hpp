#ifndef FLIPPY_VEC3_HPP
#define FLIPPY_VEC3_HPP
/**
 * @file
 * @brief Header file containing the definition and implementation a 3 dimensional vector class, with useful
 * mathematical operations like cross and dot products as member methods.
 */


#include <ostream>
#include <iostream>
#include <cmath>
#include "custom_concepts.hpp"


namespace fp{

/**
 * \brief Internal implementation of a 3D vector.
 *
 * !!! vec3 does not throw !!! This means that if you ask vec3 to divide a vector by 0 or more realistically if you
 * normalize a zero length vector vec3 will not check for the division by zero and will return a nan result!
 * Since vec3 is used everywhere in flippy, including in very expensive calculations, I decided to omit the security check
 * for the sake of speed.
 *
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
 *
 * @tparam Real @RealStub
 */

template<floating_point_number Real>
class vec3
{
public:

    Real x; //!< The x component of the vector.
    Real y; //!< The y component of the vector.
    Real z; //!< The z component of the vector.

    //! In place addition method.
    /**
     * Example:
     * ```c++
     * fp::vec3<double> v1{1,0,0};
     * fp::vec3<double> v2{0,0,1};
     * v1.add(v2);  // v1 will contain {1, 0, 1}
     * ```
     * @param v add this vector elementwise to the vector that is calling the *add* method.
     */
    void add(vec3<Real> const& v)
    {
        x += v.x;
        y += v.y;
        z += v.z;
    }

    //! In place subtraction method.
    /**
     * Example:
     * ```c++
     * fp::vec3<double> v1{2,0,0};
     * fp::vec3<double> v2{1,0,1};
     * v1.subtract(v2);  // v1 will contain {1, 0, -1}
     * ```
     * @param v subtract this vector elementwise from the vector that is calling the *subtract* method.
     */
    void subtract(vec3<Real> const& v)
    {
        x -= v.x;
        y -= v.y;
        z -= v.z;
    }
    //! Scale the vector by a real number s.
    /**
     * This function scales the vector in-place by the provided number `s`.
     * @param s multiplicative prefactor.
     */
    void scale(Real s)
    {
        x = s*x;
        y = s*y;
        z = s*z;
    }

    //! Calculate dot product with another vector.
    /**
     * Example:
     * @code{c++}
     * fp::vec3<double> v1{1,0,0};
     * fp::vec3<double> v2{2,0,1};
     * double res = v1.dot(v2);  // res will contain 2*1 + 0*0 + 0*1=2
     * @endcode
     * @param v the other vec3 vector
     * @return result of the dot product between the original vector and `v`.
     */
    Real dot(vec3<Real> const& v) const
    {
        Real res = x*v.x + y*v.y + z*v.z;
        return res;
    }

    //! Always returns 3.
    /**
     * This function always returns 3 since vec3 can only have three elements.
     * It was implemented for completeness, to make it more easy for vec3 to be used as a drop-in replacement for other vector types.
     * @return Size (number of elements) of vec3.
     */
    [[nodiscard]] constexpr std::size_t size() const { return 3; }

    //! Calculate cross product between two vectors.
    /**
     * A static method to calculate cross product between two vectors.
     * Example:
     * @code{c++}
     * fp::vec3<double> v1{1,0,0};
     * fp::vec3<double> v2{0,1,0};
     * fp::vec3<double> v3 = cross(v1, v2);  // v3 will contain {0,0,1}
     * @endcode
     * @param a first vector of the cross product
     * @param b second vector of the cross product
     * @return result of the cross product between the original vector and `v`.
     */
    static inline vec3<Real> cross(vec3<Real> const& a, vec3<Real> const& b)
    {
        vec3<Real> res;
        res.x = a.y*b.z - a.z*b.y;
        res.y = a.z*b.x - a.x*b.z;
        res.z = a.x*b.y - a.y*b.x;
        return res;
    }

    //! Calculate cross product with another vector.
    /**
     * Example:
     * @code{c++}
     * fp::vec3<double> v1{1,0,0};
     * fp::vec3<double> v2{0,1,0};
     * fp::vec3<double> v3 = v1.cross(v2);  // v3 will contain {0,0,1}
     * @endcode
     * @param other the other vec3 vector.
     * @return result of the cross product between the original vector and `other`.
     */
    vec3<Real> cross(vec3<Real> const& other) const { return cross(*this, other); }

    //! Returns the norm of the vector.
    /**
     * Example:
     * @code{c++}
     * fp::vec3<double> v{1,0,1};
     * double res = v.norm();  // res will contain 1,4142135624... i.e. sqrt(2)
     * @endcode
     * @return The euclidian norm of the vector.
     */
    Real norm() const { return std::sqrt(this->dot(*this)); }

    //! Returns the square of the norm of the vector.
    /**
     * Example:
     * @code{c++}
     * fp::vec3<double> v{1,0,1};
     * double res = v.norm_square();  // res will contain 2
     * @endcode
     * @return Square of the euclidian norm of the vector.
     */
    Real norm_square() const { return this->dot(*this); }

    //! Normalize the vector in place. And return a reference to the new normalized vector.
    /**
     * @warning If you normalize a zero length vector, you effectively
     * demand to divide by zero! this function will not do a security check
     * for you and will just return nan!
     * @return Reference to the normalized vector.
     */
    vec3<Real>const& normalize(){
        *this= *this/this->norm();
        return *this;
    }

    //! Streaming operator for easy printing of the vector.
    friend std::ostream& operator<<(std::ostream& os, const vec3<Real>& obj)
    {
        os << "{" << obj.x << ',' << obj.y << ',' << obj.z << '}';
        return os;
    }

    //! default equality operator.
    /**
     * @param other  vec3 on the right hand side of the comparison operator.
     * @return `true` if all elements of the compared vectors are equal and to `false` otherwise.
     */
    bool operator==(vec3<Real> const& other) const =default;


    //! Overloaded operator defined in terms of vec2::add.
    /**
     *
     * @param lhs left hand side of the `+` operator
     * @param rhs right hand side oif the `+` operator
     * @return equivalent to a new copy of `lhs.add(rhs)`.
     */
    friend vec3<Real> operator+(vec3<Real> lhs, vec3<Real> const& rhs)
    {
        lhs+=rhs;
        return lhs;
    }

    //! Overloaded operator defined in terms of vec3::add.
    /**
     * Equivalent to `lhs.add(rhs)`.
     * @param lhs left hand side of the `+=` operator
     * @param rhs right hand side oif the `+=` operator
     */
    friend void operator+=(vec3<Real>& lhs, vec3<Real> const& rhs)
    {
        lhs.add(rhs);
    }

    //! Overloaded operator defined in terms of vec3::subtract.
    /**
     *
     * @param lhs left hand side of the `-` operator
     * @param rhs right hand side oif the `-` operator
     * @return equivalent to a new copy of `lhs.subtract(rhs)`.
     */
    friend vec3<Real> operator-(vec3<Real> lhs, vec3<Real> const& rhs)
    {
        lhs-=rhs;
        return lhs;
    }

    //! Overloaded operator defined in terms of vec3::subtract.
    /**
     * Equivalent to `lhs.subtract(rhs)`.
     * @param lhs left hand side of the `-=` operator
     * @param rhs right hand side oif the `-=` operator
     */
    friend void operator-=(vec3<Real>& lhs, vec3<Real> const& rhs)
    {
        lhs.subtract(rhs);
    }

    //! Overloaded operator defined in terms of vec3::scale.
    /**
     * Left multiplication by a scalar `s*v`.
     * @param lhs left hand side of the `*` operator
     * @param rhs right hand side oif the `*` operator
     * @return equivalent to a new copy of `rhs.scale(lhs)`.
     */
    friend vec3<Real> operator*(Real const& lhs, vec3<Real> rhs)
    {
        rhs.scale(lhs);
        return rhs;
    }

    //! Overloaded operator defined in terms of vec3::scale.
    /**
     * Right multiplication by a scalar `v*s`.
     * @param lhs left hand side of the `*` operator
     * @param rhs right hand side oif the `*` operator
     * @return equivalent to a new copy of `lhs.scale(rhs)`.
     */
    friend vec3<Real> operator*(vec3<Real> lhs, Real const& rhs)
    {
        lhs.scale(rhs);
        return lhs;
    }

    //! Overloaded operator defined in terms of vec3::scale.
    /**
     * In place division by a scalar `v/s`,  equivalent to `lhs.scale(1/rhs)`.
     * @param lhs left hand side of the `/=` operator
     * @param rhs right hand side oif the `/=` operator
     * @warning for performance reasons, this function will not check for zero division!
     */
    friend void operator/=(vec3<Real>& lhs, Real const& rhs){
        lhs.scale(static_cast<Real>(1.)/rhs);
    }

    //! Overloaded operator defined in terms of vec3::scale.
    /**
     * Division by a scalar `v/s`.
     * @param lhs left hand side of the `/` operator
     * @param rhs right hand side oif the `/` operator
     * @return equivalent to a new copu of `lhs.scale(1/rhs)`.
     * @warning for performance reasons, this function will not check for zero division!
     */
    friend vec3<Real> operator/(vec3<Real> lhs, Real const& rhs)
    {
        lhs/=rhs;
        return lhs;
    }

    //! element access operator.
    /**
     * @tparam Index automatically deduced type of the index.
     * @param idx can only be 0 1 or 2. Any other number will cause the program to exit with an error.
     * @return for a vec3 v: v[1] returns v.x, v[2] returns v.y and v[3] returns v.z.
     *
     * @note: The use of the subscription operator might be slower than the direct access of the data member.
     */
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

    //! element access operator for constant environments.
    /**
     * @tparam Index automatically deduced type of the index.
     * @param idx can only be 0 1 or 2. Any other number will cause the program to exit with an error.
     * @return for a vec3 v: v[1] returns a constant reference to v.x, v[2] returns a constant reference to v.y and v[3] returns a constant reference to v.z.
     *
     * @note: The use of the subscription operator might be slower than the direct access of the data member.
     */
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

    //! Unary minus operator.
    /**
     *
     * @param v original vector.
     * @return A copy of -v the vector v itself stays unaffected.
     */
    friend vec3<Real> operator-(vec3<Real> v)
    {
        v.x = -v.x;
        v.y = -v.y;
        v.z = -v.z;
        return v;
    }

    //! Unary minus operator for rvalues.
    /**
     *
     * @param v an rvalue vec3 vector.
     * @return The rvalue vector `v` is moved into the function and `-v` is returned.
     */
    friend vec3<Real> operator-(vec3<Real>&& v)
    {
        v.x = -v.x;
        v.y = -v.y;
        v.z = -v.z;
        return v;
    }

};
}

#endif //FLIPPY_VEC3_HPP
