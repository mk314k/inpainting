#include "inpaint/poisson.hpp"
#include <cmath>

namespace inpaint {

Canvas threshold_mask(const Canvas& raw_mask, float threshold) {
  Canvas mask(raw_mask.width(), raw_mask.height(), 1, 0.0f);
  for (int y = 0; y < raw_mask.height(); ++y) {
    for (int x = 0; x < raw_mask.width(); ++x) {
      mask.at(x, y, 0) = (raw_mask.at(x, y, 0) >= threshold) ? 1.0f : 0.0f;
    }
  }
  return mask;
}

Canvas naive_composite(const Canvas& background,
                       const Canvas& foreground,
                       const Canvas& mask) {
  imglib::ensure_same_shape(background, foreground);
  const Canvas bin_m = threshold_mask(mask);
  Canvas out(background.width(), background.height(), background.channels());
  for (int y = 0; y < out.height(); ++y) {
    for (int x = 0; x < out.width(); ++x) {
      const float m = bin_m.at(x, y, 0);
      for (int c = 0; c < out.channels(); ++c) {
        out.at(x, y, c) = m * foreground.at(x, y, c) + (1.0f - m) * background.at(x, y, c);
      }
    }
  }
  return out;
}

PoissonResult poisson_steepest_descent(const Canvas& background,
                                       const Canvas& foreground,
                                       const Canvas& mask,
                                       int max_iters,
                                       double tol) {
  imglib::ensure_same_shape(background, foreground);
  const Canvas bin_m = threshold_mask(mask);
  Canvas x = background;
  const Canvas b = apply_laplacian(foreground);
  Canvas r = apply_masked(b - apply_laplacian(x), bin_m);

  std::vector<double> history;
  history.reserve(static_cast<std::size_t>(max_iters + 1));

  for (int k = 0; k < max_iters; ++k) {
    const double rho = dot_product(r, r);
    const double r_norm = std::sqrt(rho);
    history.push_back(r_norm);
    if (r_norm < tol) break;

    const Canvas ar = apply_masked(apply_laplacian(r), bin_m);
    const double denom = std::max(1e-12, dot_product(r, ar));
    const float alpha = static_cast<float>(rho / denom);
    x = x + r * alpha;
    r = apply_masked(b - apply_laplacian(x), bin_m);
  }

  return {x.clamped(), history};
}

PoissonResult poisson_cg(const Canvas& background,
                         const Canvas& foreground,
                         const Canvas& mask,
                         int max_iters,
                         double tol,
                         const Canvas* initial_guess) {
  imglib::ensure_same_shape(background, foreground);
  const Canvas bin_m = threshold_mask(mask);
  const MaskIndexList active = extract_mask_indices(bin_m);
  if (active.empty()) {
    return {background.clamped(), {0.0}};
  }

  Canvas x = background;
  if (initial_guess != nullptr) {
    x = naive_composite(background, *initial_guess, bin_m);
  }

  const Canvas b = apply_sparse_laplacian(foreground, active);
  Canvas r = apply_masked(b - apply_sparse_laplacian(x, active), bin_m);
  Canvas d = r;
  double rho = masked_dot_product(r, r, active);

  std::vector<double> history;
  history.reserve(static_cast<std::size_t>(max_iters + 1));
  history.push_back(std::sqrt(rho));

  for (int k = 0; k < max_iters; ++k) {
    if (std::sqrt(rho) < tol) break;
    const Canvas q = apply_sparse_laplacian(d, active);
    const double dq = std::max(1e-12, masked_dot_product(d, q, active));
    const float alpha = static_cast<float>(rho / dq);

    masked_axpy(x, d, alpha, active);
    masked_axpy(r, q, -alpha, active);

    const double rho_next = masked_dot_product(r, r, active);
    history.push_back(std::sqrt(rho_next));
    if (std::sqrt(rho_next) < tol) break;

    const float beta = static_cast<float>(rho_next / std::max(1e-12, rho));
    for (const PixelCoord& pt : active) {
      for (int c = 0; c < d.channels(); ++c) {
        d.at(pt.x, pt.y, c) = r.at(pt.x, pt.y, c) + beta * d.at(pt.x, pt.y, c);
      }
    }
    rho = rho_next;
  }

  return {x.clamped(), history};
}

PoissonResult poisson_mixed_cg(const Canvas& background,
                               const Canvas& foreground,
                               const Canvas& mask,
                               int max_iters,
                               double tol) {
  imglib::ensure_same_shape(background, foreground);
  const int w = background.width();
  const int h = background.height();
  const int ch = background.channels();
  const Canvas bin_m = threshold_mask(mask);
  const MaskIndexList active = extract_mask_indices(bin_m);
  if (active.empty()) {
    return {background.clamped(), {0.0}};
  }

  // Construct structure-dominant mixed-gradient divergence G^T v*
  Canvas vx(w, h, ch, 0.0f);
  Canvas vy(w, h, ch, 0.0f);
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      for (int c = 0; c < ch; ++c) {
        const float bx =
            background.sample(x + 1, y, c, BorderMode::ClampEdge) - background.at(x, y, c);
        const float by =
            background.sample(x, y + 1, c, BorderMode::ClampEdge) - background.at(x, y, c);
        const float fx =
            foreground.sample(x + 1, y, c, BorderMode::ClampEdge) - foreground.at(x, y, c);
        const float fy =
            foreground.sample(x, y + 1, c, BorderMode::ClampEdge) - foreground.at(x, y, c);

        if (bx * bx + by * by > fx * fx + fy * fy) {
          vx.at(x, y, c) = bx;
          vy.at(x, y, c) = by;
        } else {
          vx.at(x, y, c) = fx;
          vy.at(x, y, c) = fy;
        }
      }
    }
  }

  Canvas b_mixed(w, h, ch, 0.0f);
  for (const PixelCoord& pt : active) {
    const int x = pt.x;
    const int y = pt.y;
    for (int c = 0; c < ch; ++c) {
      const float div_x = vx.at(x, y, c) - ((x > 0) ? vx.at(x - 1, y, c) : 0.0f);
      const float div_y = vy.at(x, y, c) - ((y > 0) ? vy.at(x, y - 1, c) : 0.0f);
      b_mixed.at(x, y, c) = -(div_x + div_y);
    }
  }

  Canvas x = background;
  Canvas r = apply_masked(b_mixed - apply_sparse_laplacian(x, active), bin_m);
  Canvas d = r;
  double rho = masked_dot_product(r, r, active);

  std::vector<double> history;
  history.reserve(static_cast<std::size_t>(max_iters + 1));
  history.push_back(std::sqrt(rho));

  for (int k = 0; k < max_iters; ++k) {
    if (std::sqrt(rho) < tol) break;
    const Canvas q = apply_sparse_laplacian(d, active);
    const double dq = std::max(1e-12, masked_dot_product(d, q, active));
    const float alpha = static_cast<float>(rho / dq);

    masked_axpy(x, d, alpha, active);
    masked_axpy(r, q, -alpha, active);

    const double rho_next = masked_dot_product(r, r, active);
    history.push_back(std::sqrt(rho_next));
    if (std::sqrt(rho_next) < tol) break;

    const float beta = static_cast<float>(rho_next / std::max(1e-12, rho));
    for (const PixelCoord& pt : active) {
      for (int c = 0; c < d.channels(); ++c) {
        d.at(pt.x, pt.y, c) = r.at(pt.x, pt.y, c) + beta * d.at(pt.x, pt.y, c);
      }
    }
    rho = rho_next;
  }

  return {x.clamped(), history};
}

} // namespace inpaint
