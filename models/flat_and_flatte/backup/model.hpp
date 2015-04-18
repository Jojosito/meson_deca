#ifndef MESON_DECA__LIB__C_LIB__MODEL_HPP
#define MESON_DECA__LIB__C_LIB__MODEL_HPP

#include <vector>

#include <stan/math/prim/mat/fun/Eigen.hpp>
#include <boost/math/tools/promotion.hpp>

#include <meson_deca/lib/c_lib/complex.hpp>
#include <meson_deca/lib/c_lib/structures/resonances.hpp>


// These variables should be adjusted manually
const int NUM_RES=3; // Number of PWA resonances
const int NUM_VAR=2; // Number of independent masses (e.g., 2 for 3-body-decay)

namespace stan {
  namespace math {

    /**
     * complex_scalar A_c(vector)
     *
     * Takes the data vector y as an argument, returns the corresponding PWA 
     * amplitude (complex number). The number res_id tells, which resonance
     * to use.
     *
     * @tparam T0__ Scalar type of the data vector
     */
    template <typename T0__>
    inline
    std::vector<typename boost::math::tools::promote_arg<T0__>::type>
    A_c(const int &res_id, const Eigen::Matrix<T0__, Eigen::Dynamic,1>& y) {

        switch (res_id) {
	// This resonance list must be adjusted manually
        case 1: return 1.0 * fct::valid(1.0*y(0,0), 0*y(1,0),
					particles::d,
					particles::pi,
					particles::pi,
					particles::pi);
	  /*return resonances::toy0_1000.value(y(0,0), y(1,0),
						       particles::d,
						       particles::pi,
						       particles::pi,
						       particles::pi);*/
	  
        case 2:  /*return complex::scalar::complex(y(1,0), 0*y(1,0));*/
	   return resonances::toy0_flatte.value(y(0,0), y(1,0),
						       particles::d,
						       particles::pi,
						       particles::pi,
						       particles::pi);

        case 3: return complex::scalar::one(y(0,0));

        default: {
            std::cout << "Fatal error: Unknown resonance occured.";
            return complex::scalar::one(y(0,0)); // Dummy return
        }
        }
    }


    /**
     * complex_vector A_cv(vector)
     *
     * Takes the data vector y as an argument, returns 
     * complex vector [A(1,y) ... A(NUM_RES, y)] of PWA amplitudes.
     */
    template <typename T0__>
    inline
    std::vector<Eigen::Matrix<typename boost::math::tools::promote_args<T0__>::type, Eigen::Dynamic, 1> >
    A_cv(const Eigen::Matrix<T0__, Eigen::Dynamic,1>& y) {

        typedef typename boost::math::tools::promote_args<T0__>::type T2;

        // Somewhat convoluted initialization of the return
        std::vector<Eigen::Matrix<T2, Eigen::Dynamic, 1> > res(2, (Eigen::Matrix<T2,Eigen::Dynamic,1> (NUM_RES)));
        for (int i = 0; i < NUM_RES; i++) {
            std::vector<T0__> tmp;
            tmp = A_c(i+1, y);
            res[0](i) = tmp[0];
            res[1](i) = tmp[1];
        }
        return res;
    }



    /**
     *
     * double f_model(vector A_y[2], vector theta[2]
     *
     * Takes two comlex vectors, returns |A_y * theta|^2
     *
     */
    template <typename T0, typename T1>
    typename boost::math::tools::promote_args<T0,T1>::type
    f_model(const std::vector<Eigen::Matrix<T0, Eigen::Dynamic, 1> >& A_r,
      const std::vector<Eigen::Matrix<T1, Eigen::Dynamic, 1> >& theta) {

      //typename boost::math::tools::promote_args<T0,T1>::type res = 0;

      return complex::scalar::abs2(
               complex::vector::sum(
                 complex::vector::mult(A_r, theta)));
    }


    /**
     *
     * double Norm(vector theta[2], matrix I[2])
     *
     * Takes complex vector theta and complex matrix I,
     * returns conj(theta)' * I * theta.
     */
    template <typename T0, typename T1>
    typename boost::math::tools::promote_args<T0,T1>::type
    Norm(const std::vector<Eigen::Matrix<T0, Eigen::Dynamic, 1> >& theta,
         const std::vector<Eigen::Matrix<T1, Eigen::Dynamic, Eigen::Dynamic> >& I) {

      typename boost::math::tools::promote_args<T0,T1>::type res = 0;
      // I * theta holder, real and imaginary part
      typename boost::math::tools::promote_args<T0,T1>::type tmp[2];

      for (int i = 0; i < NUM_RES; i++) {
      for (int j = 0; j < NUM_RES; j++) {
        // Complex multiplication
        tmp[0] = I[0](i,j) * theta[0](j) - I[1](i,j) * theta[1](j);
        tmp[1] = I[0](i,j) * theta[1](j) + I[1](i,j) * theta[0](j);
          // Keep only the real palt of the product; imaginary part
          // should be 0 (+- float calculation errors).
          // Note that Re(conj(a)*b) = a[0] * b[0] + a[1] * b[1]
        res = res + theta[0](i) * tmp[0] + theta[1](i) * tmp[1];
        }
      }

      return res;
    }


    /**
     * int num_resonances()
     *
     * Returns the number of resonances
     */
    inline int num_resonances() {
        return NUM_RES;
    }

    /**
     * int num_variables()
     *
     * Returns the number of variables
     */
    inline int num_variables() {
        return NUM_VAR;
    }
  }
}

#endif
