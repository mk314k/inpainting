#include "inpaint/operators.hpp"
#include <cmath>

namespace inpaint {

double dot_product(const Canvas& a, const Canvas& b) {
  imglib::ensure_same_shape(a, b);
  double sum = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) {
    sum += static_cast<double>(a.flat(i)) * static_cast<double>(b.flat(i));
  }
  return sum;
}

double norm_l2(const Canvas& a) {
  return std::sqrt(dot_product(a, a));
}

double compute_psnr(const Canvas& ref, const Canvas& est) {
  imglib::ensure_same_shape(ref, est);
  const Canvas r_clamped = ref.clamped();
  const Canvas e_clamped = est.clamped();
  double mse = 0.0;
  for (std::size_t i = 0; i < r_clamped.size(); ++i) {
    const double diff =
        static_cast<double>(r_clamped.flat(i)) - static_cast<double>(e_clamped.flat(i));
    mse += diff * diff;
  }
  mse /= static_cast<double>(std::max<std::size_t>(1, r_clamped.size()));
  if (mse < 1e-12) return 99.0;
  return 10.0 * std::log10(1.0 / mse);
}

Kernel2D laplacian_stencil() {
  return Kernel2D(3, 3, {0.0f, -1.0f, 0.0f, -1.0f, 4.0f, -1.0f, 0.0f, -1.0f, 0.0f});
}

Canvas apply_laplacian(const Canvas& src, BorderMode border) {
  return laplacian_stencil().convolve(src, border);
}

Canvas apply_forward_conv(const Canvas& src,
                          const Kernel2D& kernel,
                          BorderMode border) {
  return kernel.convolve(src, border);
}

Canvas apply_adjoint_conv(const Canvas& src,
                          const Kernel2D& kernel,
                          BorderMode border) {
  Kernel2D reflected(kernel.width(), kernel.height());
  for (int y = 0; y < kernel.height(); ++y) {
    for (int x = 0; x < kernel.width(); ++x) {
      reflected(kernel.width() - 1 - x, kernel.height() - 1 - y) = kernel(x, y);
    }
  }
  return reflected.convolve(src, border);
}

Canvas apply_masked(const Canvas& src, const Canvas& mask) {
  if (src.width() != mask.width() || src.height() != mask.height()) {
    throw std::invalid_argument("apply_masked requires matching width and height");
  }
  Canvas out(src.width(), src.height(), src.channels(), 0.0f);
  for (int y = 0; y < src.height(); ++y) {
    for (int x = 0; x < src.width(); ++x) {
      const float m = (mask.at(x, y, 0) >= 0.5f) ? 1.0f : 0.0f;
      for (int c = 0; c < src.channels(); ++c) {
        out.at(x, y, c) = m * src.at(x, y, c);
      }
    }
  }
  return out;
}

AnisotropicField build_anisotropic_field(const Canvas& guide_img,
                                         float kappa,
                                         float sigma_grad,
                                         float sigma_window) {
  const Canvas j_tensor =
      imglib::structure_tensor(guide_img, sigma_grad, sigma_window / std::max(0.1f, sigma_grad));
  const int w = guide_img.width();
  const int h = guide_img.height();

  AnisotropicField field;
  field.width = w;
  field.height = h;
  field.d_xx.resize(static_cast<std::size_t>(w) * h, 1.0f);
  field.d_xy.resize(static_cast<std::size_t>(w) * h, 0.0f);
  field.d_yy.resize(static_cast<std::size_t>(w) * h, 1.0f);

  const float k2 = std::max(1e-6f, kappa * kappa);
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      const float jxx = j_tensor.at(x, y, 0);
      const float jxy = j_tensor.at(x, y, 1);
      const float jyy = j_tensor.at(x, y, 2);

      const imglib::SymmetricEigen2x2 eig = imglib::solve_symmetric_2x2(jxx, jxy, jyy);
      const float g_normal = 1.0f / (1.0f + std::max(0.0f, eig.lambda_max) / k2);
      const float g_tangent = 1.0f;

      const float nx = std::cos(eig.major_angle);
      const float ny = std::sin(eig.major_angle);
      const float tx = -ny;
      const float ty = nx;

      const std::size_t idx = static_cast<std::size_t>(y) * w + x;
      field.d_xx[idx] = g_normal * nx * nx + g_tangent * tx * tx;
      field.d_xy[idx] = g_normal * nx * ny + g_tangent * tx * ty;
      field.d_yy[idx] = g_normal * ny * ny + g_tangent * ty * ty;
    }
  }
  return field;
}

Canvas apply_anisotropic_laplacian(const Canvas& src,
                                   const AnisotropicField& field) {
  const int w = src.width();
  const int h = src.height();
  const int ch = src.channels();
  Canvas out(w, h, ch, 0.0f);

  // Compute forward differences gx(x,y) = u(x+1,y) - u(x,y), gy(x,y) = u(x,y+1) - u(x,y)
  // and flux J = D * grad(u), then apply adjoint G^T J
  Canvas flux_x(w, h, ch, 0.0f);
  Canvas flux_y(w, h, ch, 0.0f);

  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      const std::size_t idx = static_cast<std::size_t>(y) * w + x;
      const float dxx = field.d_xx[idx];
      const float dxy = field.d_xy[idx];
      const float dyy = field.d_yy[idx];

      for (int c = 0; c < ch; ++c) {
        const float u0 = src.at(x, y, c);
        const float ux = src.sample(x + 1, y, c, BorderMode::ClampEdge) - u0;
        const float uy = src.sample(x, y + 1, c, BorderMode::ClampEdge) - u0;
        flux_x.at(x, y, c) = dxx * ux + dxy * uy;
        flux_y.at(x, y, c) = dxy * ux + dyy * uy;
      }
    }
  }

  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      for (int c = 0; c < ch; ++c) {
        const float div_x =
            flux_x.at(x, y, c) - ((x > 0) ? flux_x.at(x - 1, y, c) : 0.0f);
        const float div_y =
            flux_y.at(x, y, c) - ((y > 0) ? flux_y.at(x, y - 1, c) : 0.0f);
        out.at(x, y, c) = -(div_x + div_y);
      }
    }
  }
  return out;
}

} // namespace inpaint
