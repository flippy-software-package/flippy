#ifndef FLIPPY_VEC3_EXPERIMENTAL_HPP
#define FLIPPY_VEC3_EXPERIMENTAL_HPP
/**
 * @file
 * @brief Header file containing the definition and implementation a 3 dimensional vector class, with useful
 * mathematical operations like cross and dot products as member methods.
 */

#include "custom_concepts.hpp"
#include <cmath>
#include <iostream>
#include <ostream>
#include <valarray>

namespace fp::experimental {

/**
 * \brief Internal implementation of a 3D vector.
 *
 * !!! vec3 does not throw !!! This means that if you ask vec3 to divide a vector by 0 or more realistically if you
 * normalize a zero length vector vec3 will not check for the division by zero and will return a nan result!
 * Since vec3 is used everywhere in flippy, including in very expensive calculations, I decided to omit the security
 * check for the sake of speed.
 *
 * To keep the external dependencies low, flippy implements it's own 3D vector class with basic functionality like dot
 * product and cross product
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

template <floating_point_number Real> class vec3 {
    std::valarray<Real> data_;

    public:
    vec3() : data_(std::valarray<Real>(3)) {}
    vec3(Real x, Real y, Real z) : data_(std::valarray{x, y, z}) {}
    vec3(Real &&x, Real &&y, Real &&z) : data_(std::valarray{x, y, z}) {}
    //    vec3(Real const&x, Real const&y, Real const&z): data_(std::valarray{x, y, z}){}

    Real x() { return data_[0]; } //!< The x component of the vector.
    Real y() { return data_[1]; } //!< The y component of the vector.
    Real z() { return data_[2]; } //!< The z component of the vector.

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
    void add(const vec3<Real> &v) { data_ += v.data_; }

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
    void subtract(const vec3<Real> &v) { data_ -= v.data_; }
    //! Scale the vector by a real number s.
    /**
     * This function scales the vector in-place by the provided number `s`.
     * @param s multiplicative prefactor.
     */
    void scale(Real s) { data_ *= s; }

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
    Real dot(const vec3<Real> &v) const { return (data_ * v.data_).sum(); }

    //! Always returns 3.
    /**
     * This function always returns 3 since vec3 can only have three elements.
     * It was implemented for completeness, to make it more easy for vec3 to be used as a drop-in replacement for other
     * vector types.
     * @return Size (number of elements) of vec3.
     */
    [[nodiscard]] constexpr std::size_t size() const { return data_.size(); }

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
    static inline vec3<Real> cross(const vec3<Real> &a, const vec3<Real> &b) {
        vec3<Real>          res;
        std::valarray<Real> as = a.data_.cshift(1);
        std::valarray<Real> bs = b.data_.cshift(-1);

        res.data_ = as * bs;

        as = as.cshift(1);
        bs = bs.cshift(-1);

        res.data_ -= as * bs;

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
    vec3<Real> cross(const vec3<Real> &other) const { return cross(*this, other); }

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
    const vec3<Real> &normalize() {
        *this = *this / this->norm();
        return *this;
    }

    //! Streaming operator for easy printing of the vector.
    friend std::ostream &operator<<(std::ostream &os, const vec3<Real> &obj) {
        os << "{" << obj.x() << ',' << obj.y() << ',' << obj.z() << '}';
        return os;
    }

    //! default equality operator.
    /**
     * @param other  vec3 on the right hand side of the comparison operator.
     * @return `true` if all elements of the compared vectors are equal and to `false` otherwise.
     */
    bool operator==(const vec3<Real> &other) const = default;

    //! Overloaded operator defined in terms of vec2::add.
    /**
     *
     * @param lhs left hand side of the `+` operator
     * @param rhs right hand side oif the `+` operator
     * @return equivalent to a new copy of `lhs.add(rhs)`.
     */
    friend vec3<Real> operator+(vec3<Real> lhs, const vec3<Real> &rhs) {
        lhs += rhs;
        return lhs;
    }

    //! Overloaded operator defined in terms of vec3::add.
    /**
     * Equivalent to `lhs.add(rhs)`.
     * @param lhs left hand side of the `+=` operator
     * @param rhs right hand side oif the `+=` operator
     */
    friend void operator+=(vec3<Real> &lhs, const vec3<Real> &rhs) { lhs.add(rhs); }

    //! Overloaded operator defined in terms of vec3::subtract.
    /**
     *
     * @param lhs left hand side of the `-` operator
     * @param rhs right hand side oif the `-` operator
     * @return equivalent to a new copy of `lhs.subtract(rhs)`.
     */
    friend vec3<Real> operator-(vec3<Real> lhs, const vec3<Real> &rhs) {
        lhs -= rhs;
        return lhs;
    }

    //! Overloaded operator defined in terms of vec3::subtract.
    /**
     * Equivalent to `lhs.subtract(rhs)`.
     * @param lhs left hand side of the `-=` operator
     * @param rhs right hand side oif the `-=` operator
     */
    friend void operator-=(vec3<Real> &lhs, const vec3<Real> &rhs) { lhs.subtract(rhs); }

    //! Overloaded operator defined in terms of vec3::scale.
    /**
     * Left multiplication by a scalar `s*v`.
     * @param lhs left hand side of the `*` operator
     * @param rhs right hand side oif the `*` operator
     * @return equivalent to a new copy of `rhs.scale(lhs)`.
     */
    friend vec3<Real> operator*(const Real &lhs, vec3<Real> rhs) {
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
    friend vec3<Real> operator*(vec3<Real> lhs, const Real &rhs) {
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
    friend void operator/=(vec3<Real> &lhs, const Real &rhs) { lhs.scale((Real)1. / rhs); }

    //! Overloaded operator defined in terms of vec3::scale.
    /**
     * Division by a scalar `v/s`.
     * @param lhs left hand side of the `/` operator
     * @param rhs right hand side oif the `/` operator
     * @return equivalent to a new copu of `lhs.scale(1/rhs)`.
     * @warning for performance reasons, this function will not check for zero division!
     */
    friend vec3<Real> operator/(vec3<Real> lhs, const Real &rhs) {
        lhs /= rhs;
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
    template <typename Index>
        requires std::is_integral_v<Index>
    Real &operator[](Index idx) {
        assert(idx < 3);
        return data_[idx];
    }

    //! element access operator for constant environments.
    /**
     * @tparam Index automatically deduced type of the index.
     * @param idx can only be 0 1 or 2. Any other number will cause the program to exit with an error.
     * @return for a vec3 v: v[1] returns a constant reference to v.x, v[2] returns a constant reference to v.y and v[3]
     * returns a constant reference to v.z.
     *
     * @note: The use of the subscription operator might be slower than the direct access of the data member.
     */
    template <typename Index>
        requires std::is_integral_v<Index>
    const Real &operator[](Index idx) const {
        assert(idx < 3);
        return data_[idx];
    }

    //! Unary minus operator.
    /**
     *
     * @param v original vector.
     * @return A copy of -v the vector v itself stays unaffected.
     */
    friend vec3<Real> operator-(vec3<Real> v) {
        v.data_ = -v.data_;
        return v;
    }

    //! Unary minus operator for rvalues.
    /**
     *
     * @param v an rvalue vec3 vector.
     * @return The rvalue vector `v` is moved into the function and `-v` is returned.
     */
    friend vec3<Real> operator-(vec3<Real> &&v) {
        v.data_ = -v.data_;
        return v;
    }
};
} // namespace fp::experimental

#endif // FLIPPY_VEC3_HPP
