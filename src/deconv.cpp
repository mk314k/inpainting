#include "inpaint/deconv.hpp"
#include <cmath>

namespace inpaint {

Canvas apply_normal_operator(const Canvas& x,
                             const Kernel2D& kernel,
                             float lambda_reg) {
  const Canvas mx = apply_forward_conv(x, kernel);
  const Canvas mt_mx = apply_adjoint_conv(mx, kernel);
  if (lambda_reg <= 0.0f) {
    return mt_mx;
  }
  return mt_mx + apply_laplacian(x) * lambda_reg;
}

DeconvResult deconvolve_steepest_descent(const Canvas& blurred,
                                         const Kernel2D& kernel,
                                         int max_iters,
                                         double tol) {
  Canvas x = blurred;
  const Canvas rhs = apply_adjoint_conv(blurred, kernel);
  std::vector<double> history;
  history.reserve(static_cast<std::size_t>(max_iters + 1));

  for (int k = 0; k < max_iters; ++k) {
    const Canvas r = rhs - apply_normal_operator(x, kernel, 0.0f);
    const double rho = dot_product(r, r);
    const double r_norm = std::sqrt(rho);
    history.push_back(r_norm);
    if (r_norm < tol) break;

    const Canvas mr = apply_forward_conv(r, kernel);
    const double denom = std::max(1e-12, dot_product(mr, mr));
    const float alpha = static_cast<float>(rho / denom);
    x = x + r * alpha;
  }

  return {x.clamped(), history};
}

DeconvResult deconvolve_cg(const Canvas& blurred,
                           const Kernel2D& kernel,
                           float lambda_reg,
                           int max_iters,
                           double tol) {
  Canvas x = blurred;
  const Canvas rhs = apply_adjoint_conv(blurred, kernel);
  Canvas r = rhs - apply_normal_operator(x, kernel, lambda_reg);
  Canvas d = r;
  double rho = dot_product(r, r);

  std::vector<double> history;
  history.reserve(static_cast<std::size_t>(max_iters + 1));
  history.push_back(std::sqrt(rho));

  for (int k = 0; k < max_iters; ++k) {
    if (std::sqrt(rho) < tol) break;
    const Canvas q = apply_normal_operator(d, kernel, lambda_reg);
    const double dq = std::max(1e-12, dot_product(d, q));
    const float alpha = static_cast<float>(rho / dq);

    x = x + d * alpha;
    r = r - q * alpha;

    const double rho_next = dot_product(r, r);
    history.push_back(std::sqrt(rho_next));
    if (std::sqrt(rho_next) < tol) break;

    const float beta = static_cast<float>(rho_next / std::max(1e-12, rho));
    d = r + d * beta;
    rho = rho_next;
  }

  return {x.clamped(), history};
}

DeconvResult deconvolve_anisotropic_cg(const Canvas& blurred,
                                       const Kernel2D& kernel,
                                       float lambda_reg,
                                       float kappa,
                                       int max_iters,
                                       double tol) {
  // Estimate edge structure tensor from Gaussian-smoothed observation
  const Canvas pre_smooth = imglib::gaussian_blur_separable(blurred, 1.0f);
  const AnisotropicField field = build_anisotropic_field(pre_smooth, kappa, 1.0f, 2.5f);

  auto apply_aniso_normal = [&](const Canvas& v) -> Canvas {
    const Canvas mv = apply_forward_conv(v, kernel);
    const Canvas mt_mv = apply_adjoint_conv(mv, kernel);
    if (lambda_reg <= 0.0f) return mt_mv;
    return mt_mv + apply_anisotropic_laplacian(v, field) * lambda_reg;
  };

  Canvas x = blurred;
  const Canvas rhs = apply_adjoint_conv(blurred, kernel);
  Canvas r = rhs - apply_aniso_normal(x);
  Canvas d = r;
  double rho = dot_product(r, r);

  std::vector<double> history;
  history.reserve(static_cast<std::size_t>(max_iters + 1));
  history.push_back(std::sqrt(rho));

  for (int k = 0; k < max_iters; ++k) {
    if (std::sqrt(rho) < tol) break;
    const Canvas q = apply_aniso_normal(d);
    const double dq = std::max(1e-12, dot_product(d, q));
    const float alpha = static_cast<float>(rho / dq);

    x = x + d * alpha;
    r = r - q * alpha;

    const double rho_next = dot_product(r, r);
    history.push_back(std::sqrt(rho_next));
    if (std::sqrt(rho_next) < tol) break;

    const float beta = static_cast<float>(rho_next / std::max(1e-12, rho));
    d = r + d * beta;
    rho = rho_next;
  }

  return {x.clamped(), history};
}

} // namespace inpaint
