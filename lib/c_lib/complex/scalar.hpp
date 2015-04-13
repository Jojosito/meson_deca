#ifndef MESON_DECA__LIB__C_LIB__COMPLEX__SCALAR_HPP
#define MESON_DECA__LIB__C_LIB__COMPLEX__SCALAR_HPP

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <vector>


/*
 *  Introduce complex number operations in a STAN-friendly way.
 *
 *  DESCRIPTION
 *    See meson_deca/lib/c_lib/complex/complex.hpp
 *
 *  FUNCTIONS
 *    complex_scalar abs2(complex_scalar)
 *    complex_scalar complex(scalar, scalar)
 *    complex_scalar inverse(complex_scalar)
 *    complex_scalar one(scalar)
 *    complex_scalar one_i(scalar)
 *    complex_scalar mult(complex_scalar, complex_scalar)
 */


namespace complex {
  namespace scalar {

    /**
     * scalar abs2(complex_scalar)
     *
     * Square magnitude of a complex number.
     * Takes the vector (a, b), returns a**2 + b**2.
     *
     * @tparam T Scalar type
     */
    template <typename T>
    inline T
    abs2(const std::vector<T> &v) {
        return v[0] * v[0] + v[1] * v[1];
    }


    /**
     * complex_scalar complex(scalar, scalar)
     *
     * Complex number constructor.
     * Takes two scalars a,b returns the vector (a, b).
     *
     * @tparam T Scalar type
     */
    template <typename T>
    inline
    std::vector<T> complex(const T& re, const T& im) {
        std::vector<T> res(2);
        res[0] = re;
        res[1] = im;
        return res;
    }


    /**
     * complex_scalar inverse(complex_scalar)
     *
     * Inverse of a complex number.
     *
     * @tparam T Scalar type
     */
    template <typename T>
    inline
    std::vector<T> inverse(const std::vector<T>& y) {

        std::vector<T> res(2);
        T norm = y[0] * y[0] + y[1] * y[1];

        res[0] = y[0] / norm;
        res[1] = -y[1] / norm;
        return res;
    }


    /**
     * complex_scalar one(scalar)
     *
     * Returns the complex number 1.
     *
     * @param y Dummy scalar variable
     * @tparam T Scalar type
     */
    template <typename T>
    inline
    std::vector<T> one(const T& y) {
        std::vector<T> res(2);
        res[0] = 1.0;
        res[1] = 0.0;
        return res;
    }


    /**
     * complex_scalar one(scalar)
     *
     * Returns the imaginary unit 1*i.
     *
     * @param y Dummy scalar variable
     * @tparam T Scalar type
     */
    template <typename T>
    inline
    std::vector<T> one_i(const T& y) {
        std::vector<T> res(2);
        res[0] = 0.0;
        res[1] = 1.0;
        return res;
    }


    /**
     * complex_scalar mult(complex_scalar, complex_scalar)
     *
     * Multiplication of two complex numbers.
     *
     * Caveat: both arguments MUST have the same scalar type.
     * Probably should be rewritten to avoid that.
     *
     * @tparam T Scalar type
     */
    template <typename T>
    inline
    std::vector<T> mult(const std::vector<T> &v1, const std::vector<T> &v2) {
        std::vector<T> res(2);
        res[0] = v1[0] * v2[0] - v1[1] * v2[1];
        res[1] = v1[1] * v2[0] + v1[0] * v2[1];
        return res;
    }

  }
}
#endif