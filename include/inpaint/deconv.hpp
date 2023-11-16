#ifndef INPAINT_DECONV_HPP
#define INPAINT_DECONV_HPP

#include "inpaint/operators.hpp"
#include <vector>

namespace inpaint {

struct DeconvResult {
  Canvas image;
  std::vector<double> residual_norms;
};

Canvas apply_normal_operator(const Canvas& x,
                             const Kernel2D& kernel,
                             float lambda_reg = 0.0f);

DeconvResult deconvolve_steepest_descent(const Canvas& blurred,
                                         const Kernel2D& kernel,
                                         int max_iters = 25,
                                         double tol = 1e-6);

DeconvResult deconvolve_cg(const Canvas& blurred,
                           const Kernel2D& kernel,
                           float lambda_reg = 0.0f,
                           int max_iters = 25,
                           double tol = 1e-6);

DeconvResult deconvolve_anisotropic_cg(const Canvas& blurred,
                                       const Kernel2D& kernel,
                                       float lambda_reg = 0.045f,
                                       float kappa = 0.06f,
                                       int max_iters = 25,
                                       double tol = 1e-6);

} // namespace inpaint

#endif // INPAINT_DECONV_HPP
