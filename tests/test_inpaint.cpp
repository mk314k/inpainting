#include "inpaint/inpaint.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

namespace {

using imglib::Canvas;
using imglib::Kernel2D;

void test_adjoint_and_laplacian_spd() {
  Canvas x(12, 10, 1, 0.0f);
  Canvas y(12, 10, 1, 0.0f);
  for (int j = 0; j < 10; ++j) {
    for (int i = 0; i < 12; ++i) {
      x.at(i, j, 0) = std::sin(0.3f * i + 0.2f * j);
      y.at(i, j, 0) = std::cos(0.4f * i - 0.1f * j);
    }
  }

  Kernel2D k(3, 3, {0.05f, 0.1f, 0.15f, 0.2f, 0.1f, 0.05f, 0.1f, 0.15f, 0.1f});
  const Canvas kx = inpaint::apply_forward_conv(x, k, imglib::BorderMode::ZeroPad);
  const Canvas kty = inpaint::apply_adjoint_conv(y, k, imglib::BorderMode::ZeroPad);
  const double lhs = inpaint::dot_product(kx, y);
  const double rhs = inpaint::dot_product(x, kty);
  assert(std::abs(lhs - rhs) < 1e-4);

  // Verify <x, L x> >= 0
  const Canvas lx = inpaint::apply_laplacian(x, imglib::BorderMode::ZeroPad);
  assert(inpaint::dot_product(x, lx) > 0.0);
}

void test_deconvolution_cg_improves_psnr() {
  Canvas gt(24, 24, 3, 0.2f);
  gt.fill_rect(6, 6, 17, 17, 0.85f, 0.4f, 0.2f);

  const std::vector<float> gw = imglib::gaussian_kernel_2d(1.2f, 2.0f);
  const int kd = static_cast<int>(std::round(std::sqrt(gw.size())));
  const Kernel2D psf(kd, kd, gw);

  const Canvas blurred = inpaint::apply_forward_conv(gt, psf);
  const double blur_psnr = inpaint::compute_psnr(gt, blurred);

  const inpaint::DeconvResult sd = inpaint::deconvolve_steepest_descent(blurred, psf, 15);
  const inpaint::DeconvResult cg = inpaint::deconvolve_cg(blurred, psf, 1e-3f, 15);

  const double sd_psnr = inpaint::compute_psnr(gt, sd.image);
  const double cg_psnr = inpaint::compute_psnr(gt, cg.image);
  assert(sd_psnr > blur_psnr);
  assert(cg_psnr >= sd_psnr);
}

void test_poisson_cloning_convergence() {
  Canvas bg(20, 20, 3, 0.3f);
  Canvas fg(20, 20, 3, 0.7f);
  fg.fill_rect(8, 8, 12, 12, 0.9f, 0.9f, 0.9f);
  Canvas mask(20, 20, 1, 0.0f);
  mask.fill_rect(5, 5, 15, 15, 1.0f);

  const inpaint::PoissonResult res = inpaint::poisson_cg(bg, fg, mask, 80);
  assert(!res.residual_norms.empty());
  assert(res.residual_norms.back() < res.residual_norms.front() * 0.05);
  // Exterior must match bg exactly
  assert(std::abs(res.image.at(1, 1, 0) - 0.3f) < 1e-4f);
}

void test_inpainting_solvers_and_hybrid_warm_start() {
  Canvas gt(28, 28, 3, 0.0f);
  for (int y = 0; y < 28; ++y) {
    for (int x = 0; x < 28; ++x) {
      const float v = (x >= 14) ? 0.8f : 0.2f;
      gt.at(x, y, 0) = v;
      gt.at(x, y, 1) = 0.5f;
      gt.at(x, y, 2) = 1.0f - v;
    }
  }

  Canvas mask(28, 28, 1, 0.0f);
  mask.fill_rect(10, 10, 18, 18, 1.0f);
  Canvas dmg = gt;
  dmg.fill_rect(10, 10, 18, 18, 0.0f, 0.0f, 0.0f);

  const inpaint::InpaintResult harm = inpaint::inpaint_harmonic_cg(dmg, mask, 25);
  const inpaint::InpaintResult fmm = inpaint::inpaint_telea_fmm(dmg, mask, 3);
  const inpaint::InpaintResult hyb =
      inpaint::inpaint_hybrid_structure_cg(dmg, mask, 0.08f, 3, 25);

  const double dmg_psnr = inpaint::compute_psnr(gt, dmg);
  assert(inpaint::compute_psnr(gt, harm.image) > dmg_psnr + 4.0);
  assert(inpaint::compute_psnr(gt, fmm.image) > dmg_psnr + 6.0);
  assert(inpaint::compute_psnr(gt, hyb.image) >= inpaint::compute_psnr(gt, harm.image));

  // Verify sparse active-mask operators match full-grid operators on masked domain
  const inpaint::MaskIndexList active = inpaint::extract_mask_indices(mask);
  assert(active.size() == 81u);
  const Canvas full_lap = inpaint::apply_masked(inpaint::apply_laplacian(gt), mask);
  const Canvas sparse_lap = inpaint::apply_sparse_laplacian(gt, active);
  for (const inpaint::PixelCoord& pt : active) {
    assert(std::abs(full_lap.at(pt.x, pt.y, 0) - sparse_lap.at(pt.x, pt.y, 0)) < 1e-5f);
  }
}

} // namespace

int main() {
  test_adjoint_and_laplacian_spd();
  test_deconvolution_cg_improves_psnr();
  test_poisson_cloning_convergence();
  test_inpainting_solvers_and_hybrid_warm_start();
  std::cout << "All Inpainting unit tests passed (4/4 suites).\n";
  return 0;
}
