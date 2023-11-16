#ifndef INPAINT_OPERATORS_HPP
#define INPAINT_OPERATORS_HPP

#include "imglib/imglib.hpp"
#include <vector>

namespace inpaint {

using imglib::BorderMode;
using imglib::Canvas;
using imglib::Kernel2D;

double dot_product(const Canvas& a, const Canvas& b);
double norm_l2(const Canvas& a);
double compute_psnr(const Canvas& ref, const Canvas& est);

Kernel2D laplacian_stencil();
Canvas apply_laplacian(const Canvas& src, BorderMode border = BorderMode::ClampEdge);

Canvas apply_forward_conv(const Canvas& src,
                          const Kernel2D& kernel,
                          BorderMode border = BorderMode::ClampEdge);

Canvas apply_adjoint_conv(const Canvas& src,
                          const Kernel2D& kernel,
                          BorderMode border = BorderMode::ClampEdge);

Canvas apply_masked(const Canvas& src, const Canvas& mask);

struct AnisotropicField {
  int width = 0;
  int height = 0;
  std::vector<float> d_xx;
  std::vector<float> d_xy;
  std::vector<float> d_yy;
};

AnisotropicField build_anisotropic_field(const Canvas& guide_img,
                                         float kappa = 0.08f,
                                         float sigma_grad = 1.0f,
                                         float sigma_window = 3.0f);

Canvas apply_anisotropic_laplacian(const Canvas& src,
                                   const AnisotropicField& field);

} // namespace inpaint

#endif // INPAINT_OPERATORS_HPP
