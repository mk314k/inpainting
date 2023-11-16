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

} // namespace inpaint
