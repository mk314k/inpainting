#include "inpaint/inpaint.hpp"
#include <cmath>
#include <queue>
#include <vector>

namespace inpaint {

namespace {

enum class FmmState : unsigned char {
  Known = 0,
  Band = 1,
  Inside = 2
};

struct HeapNode {
  float time;
  int x;
  int y;
  bool operator>(const HeapNode& other) const { return time > other.time; }
};

float solve_eikonal_upwind(int x,
                           int y,
                           const Canvas& T,
                           const std::vector<FmmState>& states,
                           int w,
                           int h) {
  constexpr float kInf = 1e6f;
  float tx = kInf;
  float ty = kInf;

  if (x > 0 && states[static_cast<std::size_t>(y) * w + (x - 1)] == FmmState::Known) {
    tx = std::min(tx, T.at(x - 1, y, 0));
  }
  if (x + 1 < w && states[static_cast<std::size_t>(y) * w + (x + 1)] == FmmState::Known) {
    tx = std::min(tx, T.at(x + 1, y, 0));
  }
  if (y > 0 && states[static_cast<std::size_t>(y - 1) * w + x] == FmmState::Known) {
    ty = std::min(ty, T.at(x, y - 1, 0));
  }
  if (y + 1 < h && states[static_cast<std::size_t>(y + 1) * w + x] == FmmState::Known) {
    ty = std::min(ty, T.at(x, y + 1, 0));
  }

  const float min_t = std::min(tx, ty);
  if (min_t >= kInf * 0.5f) return kInf;

  const float diff = std::abs(tx - ty);
  if (diff < 1.0f) {
    return 0.5f * (tx + ty + std::sqrt(std::max(0.0f, 2.0f - diff * diff)));
  }
  return min_t + 1.0f;
}

void inpaint_fmm_pixel(int px,
                       int py,
                       Canvas& img,
                       const Canvas& T,
                       const std::vector<FmmState>& states,
                       int radius) {
  const int w = img.width();
  const int h = img.height();
  const int ch = img.channels();

  // Estimate normal grad(T) at (px, py)
  const float tx_l = (px > 0 && states[static_cast<std::size_t>(py) * w + px - 1] != FmmState::Inside)
                         ? T.at(px - 1, py, 0)
                         : T.at(px, py, 0);
  const float tx_r = (px + 1 < w && states[static_cast<std::size_t>(py) * w + px + 1] != FmmState::Inside)
                         ? T.at(px + 1, py, 0)
                         : T.at(px, py, 0);
  const float ty_u = (py > 0 && states[static_cast<std::size_t>(py - 1) * w + px] != FmmState::Inside)
                         ? T.at(px, py - 1, 0)
                         : T.at(px, py, 0);
  const float ty_d = (py + 1 < h && states[static_cast<std::size_t>(py + 1) * w + px] != FmmState::Inside)
                         ? T.at(px, py + 1, 0)
                         : T.at(px, py, 0);

  float gx_t = 0.5f * (tx_r - tx_l);
  float gy_t = 0.5f * (ty_d - ty_u);
  const float g_norm = std::sqrt(gx_t * gx_t + gy_t * gy_t);
  if (g_norm > 1e-6f) {
    gx_t /= g_norm;
    gy_t /= g_norm;
  } else {
    gx_t = 1.0f;
    gy_t = 0.0f;
  }

  std::vector<float> numer(static_cast<std::size_t>(ch), 0.0f);
  float denom = 0.0f;
  const float tp = T.at(px, py, 0);

  for (int dy = -radius; dy <= radius; ++dy) {
    const int qy = py + dy;
    if (qy < 0 || qy >= h) continue;
    for (int dx = -radius; dx <= radius; ++dx) {
      if (dx == 0 && dy == 0) continue;
      const int qx = px + dx;
      if (qx < 0 || qx >= w) continue;
      if (states[static_cast<std::size_t>(qy) * w + qx] != FmmState::Known) continue;

      const float rx = static_cast<float>(-dx);
      const float ry = static_cast<float>(-dy);
      const float d2 = rx * rx + ry * ry;
      const float d_len = std::sqrt(d2);

      const float dir = std::max(1e-4f, std::abs(rx * gx_t + ry * gy_t) / d_len);
      const float dst = 1.0f / d2;
      const float lev = 1.0f / (1.0f + std::abs(tp - T.at(qx, qy, 0)));
      const float w_pq = dir * dst * lev;

      for (int c = 0; c < ch; ++c) {
        float val = img.at(qx, qy, c);
        // Include first-order Taylor gradient when interior to KNOWN region
        if (qx > 0 && qx + 1 < w && qy > 0 && qy + 1 < h &&
            states[static_cast<std::size_t>(qy) * w + qx - 1] == FmmState::Known &&
            states[static_cast<std::size_t>(qy) * w + qx + 1] == FmmState::Known &&
            states[static_cast<std::size_t>(qy - 1) * w + qx] == FmmState::Known &&
            states[static_cast<std::size_t>(qy + 1) * w + qx] == FmmState::Known) {
          const float ix = 0.5f * (img.at(qx + 1, qy, c) - img.at(qx - 1, qy, c));
          const float iy = 0.5f * (img.at(qx, qy + 1, c) - img.at(qx, qy - 1, c));
          const float taylor = std::max(-0.15f, std::min(0.15f, ix * rx + iy * ry));
          val = std::max(0.0f, std::min(1.0f, val + taylor));
        }
        numer[static_cast<std::size_t>(c)] += w_pq * val;
      }
      denom += w_pq;
    }
  }

  if (denom > 1e-12f) {
    for (int c = 0; c < ch; ++c) {
      img.at(px, py, c) = numer[static_cast<std::size_t>(c)] / denom;
    }
  }
}

} // namespace

InpaintResult inpaint_harmonic_cg(const Canvas& corrupted,
                                  const Canvas& mask,
                                  int max_iters,
                                  double tol) {
  const Canvas bin_m = threshold_mask(mask);
  const Canvas zero_fg(corrupted.width(), corrupted.height(), corrupted.channels(), 0.0f);
  Canvas init_cold = corrupted;
  for (int y = 0; y < init_cold.height(); ++y) {
    for (int x = 0; x < init_cold.width(); ++x) {
      if (bin_m.at(x, y, 0) >= 0.5f) {
        for (int c = 0; c < init_cold.channels(); ++c) {
          init_cold.at(x, y, c) = 0.0f;
        }
      }
    }
  }
  const PoissonResult res =
      poisson_cg(init_cold, zero_fg, bin_m, max_iters, tol);
  return {res.image, Canvas(corrupted.width(), corrupted.height(), 1, 0.0f), res.residual_norms};
}

InpaintResult inpaint_telea_fmm(const Canvas& corrupted,
                                const Canvas& mask,
                                int radius) {
  const int w = corrupted.width();
  const int h = corrupted.height();
  const Canvas bin_m = threshold_mask(mask);

  Canvas out = corrupted;
  Canvas T(w, h, 1, 0.0f);
  std::vector<FmmState> states(static_cast<std::size_t>(w) * h, FmmState::Known);

  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      if (bin_m.at(x, y, 0) >= 0.5f) {
        states[static_cast<std::size_t>(y) * w + x] = FmmState::Inside;
        T.at(x, y, 0) = 1e6f;
        for (int c = 0; c < out.channels(); ++c) out.at(x, y, c) = 0.0f;
      }
    }
  }

  std::priority_queue<HeapNode, std::vector<HeapNode>, std::greater<HeapNode>> pq;
  const int dx4[4] = {-1, 1, 0, 0};
  const int dy4[4] = {0, 0, -1, 1};

  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      if (states[static_cast<std::size_t>(y) * w + x] != FmmState::Inside) continue;
      bool touches_known = false;
      for (int k = 0; k < 4; ++k) {
        const int nx = x + dx4[k];
        const int ny = y + dy4[k];
        if (nx >= 0 && nx < w && ny >= 0 && ny < h &&
            states[static_cast<std::size_t>(ny) * w + nx] == FmmState::Known) {
          touches_known = true;
          break;
        }
      }
      if (touches_known) {
        states[static_cast<std::size_t>(y) * w + x] = FmmState::Band;
        const float t_val = solve_eikonal_upwind(x, y, T, states, w, h);
        T.at(x, y, 0) = t_val;
        inpaint_fmm_pixel(x, y, out, T, states, std::max(1, radius));
        pq.push({t_val, x, y});
      }
    }
  }

  while (!pq.empty()) {
    const HeapNode curr = pq.top();
    pq.pop();
    const std::size_t curr_idx = static_cast<std::size_t>(curr.y) * w + curr.x;
    if (states[curr_idx] == FmmState::Known) continue;
    states[curr_idx] = FmmState::Known;

    for (int k = 0; k < 4; ++k) {
      const int nx = curr.x + dx4[k];
      const int ny = curr.y + dy4[k];
      if (nx < 0 || nx >= w || ny < 0 || ny >= h) continue;
      const std::size_t n_idx = static_cast<std::size_t>(ny) * w + nx;
      if (states[n_idx] == FmmState::Known) continue;

      const float t_new = solve_eikonal_upwind(nx, ny, T, states, w, h);
      if (states[n_idx] == FmmState::Inside) {
        T.at(nx, ny, 0) = t_new;
        inpaint_fmm_pixel(nx, ny, out, T, states, std::max(1, radius));
        states[n_idx] = FmmState::Band;
        pq.push({t_new, nx, ny});
      } else if (t_new < T.at(nx, ny, 0)) {
        T.at(nx, ny, 0) = t_new;
        pq.push({t_new, nx, ny});
      }
    }
  }

  // Normalize arrival time map to [0, 1] for visualization
  float max_t = 0.0f;
  for (int y = 0; y < h; ++y) {
    for (int x = 0; x < w; ++x) {
      if (bin_m.at(x, y, 0) >= 0.5f && T.at(x, y, 0) < 1e5f) {
        max_t = std::max(max_t, T.at(x, y, 0));
      }
    }
  }
  if (max_t > 1e-6f) {
    T = T / max_t;
  }

  return {out.clamped(), T, {}};
}

InpaintResult inpaint_hybrid_structure_cg(const Canvas& corrupted,
                                          const Canvas& mask,
                                          float kappa,
                                          int fmm_radius,
                                          int cg_iters,
                                          double tol) {
  const Canvas bin_m = threshold_mask(mask);

  // Stage 1: Telea Eikonal Fast Marching warm start
  const InpaintResult fmm = inpaint_telea_fmm(corrupted, bin_m, fmm_radius);

  // Stage 2: Build structure-tensor anisotropic conductivity field D(p)
  const AnisotropicField field = build_anisotropic_field(fmm.image, kappa, 1.0f, 3.0f);

  // Stage 3: Anisotropic Conjugate Gradient relaxation from FMM warm start
  Canvas x = fmm.image;
  Canvas r = apply_masked(Canvas(x.width(), x.height(), x.channels(), 0.0f) -
                              apply_anisotropic_laplacian(x, field),
                          bin_m);
  Canvas d = r;
  double rho = dot_product(r, r);

  std::vector<double> history;
  history.reserve(static_cast<std::size_t>(cg_iters + 1));
  history.push_back(std::sqrt(rho));

  for (int k = 0; k < cg_iters; ++k) {
    if (std::sqrt(rho) < tol) break;
    const Canvas q = apply_masked(apply_anisotropic_laplacian(d, field), bin_m);
    const double dq = std::max(1e-12, dot_product(d, q));
    const float alpha = static_cast<float>(rho / dq);

    x = x + d * alpha;
    r = apply_masked(r - q * alpha, bin_m);

    const double rho_next = dot_product(r, r);
    history.push_back(std::sqrt(rho_next));
    if (std::sqrt(rho_next) < tol) break;

    const float beta = static_cast<float>(rho_next / std::max(1e-12, rho));
    d = r + d * beta;
    rho = rho_next;
  }

  return {x.clamped(), fmm.arrival_time, history};
}

} // namespace inpaint
