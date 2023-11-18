#include "inpaint/inpaint.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>

namespace {

using imglib::Canvas;
using imglib::Kernel2D;

Canvas build_reference_scene(int w = 128, int h = 128) {
  Canvas scene(w, h, 3, 0.0f);
  for (int y = 0; y < h; ++y) {
    const float yf = static_cast<float>(y) / static_cast<float>(h);
    for (int x = 0; x < w; ++x) {
      const float xf = static_cast<float>(x) / static_cast<float>(w);
      scene.at(x, y, 0) = 0.14f + 0.26f * yf;
      scene.at(x, y, 1) = 0.22f + 0.34f * yf;
      scene.at(x, y, 2) = 0.44f + 0.36f * (1.0f - yf);

      // Sun disc
      const float dx_s = static_cast<float>(x - 92);
      const float dy_s = static_cast<float>(y - 30);
      if (dx_s * dx_s + dy_s * dy_s <= 17.0f * 17.0f) {
        scene.at(x, y, 0) = 0.96f;
        scene.at(x, y, 1) = 0.82f;
        scene.at(x, y, 2) = 0.28f;
      }

      // Mountain ridges with sharp diagonal edges
      if (y >= static_cast<int>(66.0f + 0.48f * std::abs(static_cast<float>(x - 42)))) {
        scene.at(x, y, 0) = 0.24f;
        scene.at(x, y, 1) = 0.36f;
        scene.at(x, y, 2) = 0.29f;
      }
      if (y >= static_cast<int>(58.0f + 0.54f * std::abs(static_cast<float>(x - 84)))) {
        scene.at(x, y, 0) = 0.17f;
        scene.at(x, y, 1) = 0.26f;
        scene.at(x, y, 2) = 0.24f;
      }

      // High-contrast resolution checkerboard in foreground
      if (y >= 84 && y < 116 && x >= 18 && x < 66) {
        const bool white = (((x - 18) / 8) + ((y - 84) / 8)) % 2 == 0;
        scene.at(x, y, 0) = white ? 0.92f : 0.14f;
        scene.at(x, y, 1) = white ? 0.90f : 0.15f;
        scene.at(x, y, 2) = white ? 0.85f : 0.20f;
      }

      // Bright ring target
      const float r_ring =
          std::sqrt(static_cast<float>((x - 94) * (x - 94) + (y - 96) * (y - 96)));
      if (r_ring >= 8.5f && r_ring <= 14.5f) {
        scene.at(x, y, 0) = 0.90f;
        scene.at(x, y, 1) = 0.32f;
        scene.at(x, y, 2) = 0.24f;
      }
      (void)xf;
    }
  }
  return scene;
}

int run_benchmark_suite(const std::string& out_dir) {
  const Canvas gt = build_reference_scene(128, 128);
  imglib::save_png(gt, out_dir + "/exp1_gt.png");

  // 1. Deconvolution benchmark (Steepest Descent vs Conjugate Gradient)
  const std::vector<float> k_weights = imglib::gaussian_kernel_2d(1.6f, 2.5f);
  const int k_dim = static_cast<int>(std::round(std::sqrt(k_weights.size())));
  const Kernel2D psf(k_dim, k_dim, k_weights);

  const Canvas blurred = inpaint::apply_forward_conv(gt, psf);
  const inpaint::DeconvResult sd_res = inpaint::deconvolve_steepest_descent(blurred, psf, 20);
  const inpaint::DeconvResult cg_res = inpaint::deconvolve_cg(blurred, psf, 2e-4f, 20);

  imglib::save_png(blurred, out_dir + "/exp1_blur.png");
  imglib::save_png(sd_res.image, out_dir + "/exp1_sd.png");
  imglib::save_png(cg_res.image, out_dir + "/exp1_cg.png");

  std::cout << std::fixed << std::setprecision(2);
  std::cout << "[Exp 1] Blur PSNR: " << inpaint::compute_psnr(gt, blurred)
            << " dB | SD (20): " << inpaint::compute_psnr(gt, sd_res.image)
            << " dB | CG (20): " << inpaint::compute_psnr(gt, cg_res.image) << " dB\n";

  // 2. Noisy Deconvolution (Unregularized vs Tikhonov-Laplacian Regularized CG)
  Canvas noisy_blur = blurred;
  std::mt19937 rng(42u);
  std::normal_distribution<float> noise_dist(0.0f, 0.032f);
  for (std::size_t i = 0; i < noisy_blur.size(); ++i) {
    noisy_blur.flat(i) = std::max(0.0f, std::min(1.0f, noisy_blur.flat(i) + noise_dist(rng)));
  }

  const inpaint::DeconvResult unreg_res = inpaint::deconvolve_cg(noisy_blur, psf, 0.0f, 18);
  const inpaint::DeconvResult reg_res = inpaint::deconvolve_cg(noisy_blur, psf, 0.045f, 18);
  const inpaint::DeconvResult aniso_reg_res =
      inpaint::deconvolve_anisotropic_cg(noisy_blur, psf, 0.045f, 0.06f, 18);

  imglib::save_png(noisy_blur, out_dir + "/exp2_noisy.png");
  imglib::save_png(unreg_res.image, out_dir + "/exp2_unreg.png");
  imglib::save_png(reg_res.image, out_dir + "/exp2_reg.png");
  imglib::save_png(aniso_reg_res.image, out_dir + "/exp2_aniso_reg.png");

  std::cout << "[Exp 2] Noisy PSNR: " << inpaint::compute_psnr(gt, noisy_blur)
            << " dB | Unreg CG: " << inpaint::compute_psnr(gt, unreg_res.image)
            << " dB | Isotropic Reg CG: " << inpaint::compute_psnr(gt, reg_res.image)
            << " dB | Anisotropic Structure CG: " << inpaint::compute_psnr(gt, aniso_reg_res.image)
            << " dB\n";

  // 3. Poisson Seamless Cloning
  Canvas bg(128, 128, 3, 0.0f);
  Canvas fg(128, 128, 3, 0.0f);
  Canvas clone_mask(128, 128, 1, 0.0f);
  for (int y = 0; y < 128; ++y) {
    const float yf = static_cast<float>(y) / 128.0f;
    for (int x = 0; x < 128; ++x) {
      bg.at(x, y, 0) = 0.12f + 0.45f * yf;
      bg.at(x, y, 1) = 0.38f + 0.22f * (1.0f - yf) + 0.03f * std::sin(y * 0.35f);
      bg.at(x, y, 2) = 0.68f - 0.22f * yf;

      fg.at(x, y, 0) = 0.76f;
      fg.at(x, y, 1) = 0.62f;
      fg.at(x, y, 2) = 0.40f;

      const float rx = (static_cast<float>(x) - 64.0f) / 30.0f;
      const float ry = (static_cast<float>(y) - 64.0f) / 24.0f;
      const float d = std::sqrt(rx * rx + ry * ry);
      if (d <= 0.92f) clone_mask.at(x, y, 0) = 1.0f;
      if (d <= 0.72f) {
        fg.at(x, y, 0) = 0.86f;
        fg.at(x, y, 1) = 0.25f;
        fg.at(x, y, 2) = 0.18f;
      }
      if (d <= 0.44f) {
        fg.at(x, y, 0) = 0.98f;
        fg.at(x, y, 1) = 0.84f;
        fg.at(x, y, 2) = 0.25f;
      }
      if (d <= 0.20f) {
        fg.at(x, y, 0) = 0.16f;
        fg.at(x, y, 1) = 0.22f;
        fg.at(x, y, 2) = 0.38f;
      }
    }
  }

  const Canvas naive_comp = inpaint::naive_composite(bg, fg, clone_mask);
  const inpaint::PoissonResult pois_res = inpaint::poisson_cg(bg, fg, clone_mask, 160);
  const inpaint::PoissonResult mixed_res = inpaint::poisson_mixed_cg(bg, fg, clone_mask, 160);

  imglib::save_png(fg, out_dir + "/exp3_fg.png");
  imglib::save_png(bg, out_dir + "/exp3_bg.png");
  imglib::save_png(naive_comp, out_dir + "/exp3_naive.png");
  imglib::save_png(pois_res.image, out_dir + "/exp3_poisson.png");
  imglib::save_png(mixed_res.image, out_dir + "/exp3_mixed.png");

  // 4. Inpainting Comparison (Cold Harmonic CG vs Telea FMM vs Hybrid Structure CG)
  Canvas inp_mask(128, 128, 1, 0.0f);
  inp_mask.stroke_line(18, 22, 110, 106, 1.0f);
  inp_mask.stroke_line(19, 22, 111, 106, 1.0f);
  inp_mask.stroke_line(18, 23, 110, 107, 1.0f);
  inp_mask.fill_rect(26, 36, 102, 42, 1.0f);
  inp_mask.fill_rect(30, 74, 96, 80, 1.0f);

  Canvas damaged = gt;
  for (int y = 0; y < 128; ++y) {
    for (int x = 0; x < 128; ++x) {
      if (inp_mask.at(x, y, 0) >= 0.5f) {
        damaged.at(x, y, 0) = 0.98f;
        damaged.at(x, y, 1) = 0.16f;
        damaged.at(x, y, 2) = 0.52f;
      }
    }
  }

  const inpaint::InpaintResult harm_10 = inpaint::inpaint_harmonic_cg(damaged, inp_mask, 10);
  const inpaint::InpaintResult harm_40 = inpaint::inpaint_harmonic_cg(damaged, inp_mask, 40);
  const inpaint::InpaintResult fmm_res = inpaint::inpaint_telea_fmm(damaged, inp_mask, 3);
  const inpaint::InpaintResult hyb_10 =
      inpaint::inpaint_hybrid_structure_cg(damaged, inp_mask, 0.08f, 3, 10);
  const inpaint::InpaintResult hyb_40 =
      inpaint::inpaint_hybrid_structure_cg(damaged, inp_mask, 0.08f, 3, 40);

  imglib::save_png(damaged, out_dir + "/exp4_damaged.png");
  imglib::save_png(harm_40.image, out_dir + "/exp4_harmonic.png");
  imglib::save_png(fmm_res.image, out_dir + "/exp4_fmm.png");
  imglib::save_png(hyb_40.image, out_dir + "/exp4_hybrid.png");

  std::cout << "[Exp 4] Inpaint PSNR — Cold Harmonic CG (10): "
            << inpaint::compute_psnr(gt, harm_10.image)
            << " dB | Cold Harmonic CG (40): " << inpaint::compute_psnr(gt, harm_40.image)
            << " dB | Telea FMM (1 pass): " << inpaint::compute_psnr(gt, fmm_res.image)
            << " dB | Hybrid FMM-Structure CG (10): " << inpaint::compute_psnr(gt, hyb_10.image)
            << " dB | Hybrid FMM-Structure CG (40): " << inpaint::compute_psnr(gt, hyb_40.image)
            << " dB\n";

  return 0;
}

void print_usage() {
  std::cout
      << "Usage:\n"
      << "  inpaint benchmark <output_dir>\n"
      << "  inpaint deconv    <blurred.png> <out.png> [sigma=1.6] [lambda=0.02] [iters=25]\n"
      << "  inpaint poisson   <bg.png> <fg.png> <mask.png> <out.png> [iters=160]\n"
      << "  inpaint harmonic  <damaged.png> <mask.png> <out.png> [iters=160]\n"
      << "  inpaint fmm       <damaged.png> <mask.png> <out.png> [radius=3]\n"
      << "  inpaint hybrid    <damaged.png> <mask.png> <out.png> [iters=60]\n";
}

} // namespace

int main(int argc, char** argv) {
  if (argc < 3) {
    print_usage();
    return 0;
  }

  const std::string cmd = argv[1];
  try {
    if (cmd == "benchmark") {
      return run_benchmark_suite(argv[2]);
    }
    if (cmd == "deconv" && argc >= 4) {
      const Canvas y = imglib::load_png(argv[2]);
      const float sigma = (argc >= 5) ? std::stof(argv[4]) : 1.6f;
      const float lamb = (argc >= 6) ? std::stof(argv[5]) : 0.02f;
      const int iters = (argc >= 7) ? std::stoi(argv[6]) : 25;
      const std::vector<float> kw = imglib::gaussian_kernel_2d(sigma, 2.5f);
      const int kdim = static_cast<int>(std::round(std::sqrt(kw.size())));
      const Kernel2D psf(kdim, kdim, kw);
      const inpaint::DeconvResult res = inpaint::deconvolve_cg(y, psf, lamb, iters);
      imglib::save_png(res.image, argv[3]);
      return 0;
    }
    if (cmd == "poisson" && argc >= 6) {
      const Canvas bg = imglib::load_png(argv[2]);
      const Canvas fg = imglib::load_png(argv[3]);
      const Canvas mask = imglib::load_png(argv[4]);
      const int iters = (argc >= 7) ? std::stoi(argv[6]) : 160;
      const inpaint::PoissonResult res = inpaint::poisson_cg(bg, fg, mask, iters);
      imglib::save_png(res.image, argv[5]);
      return 0;
    }
    if ((cmd == "harmonic" || cmd == "fmm" || cmd == "hybrid") && argc >= 5) {
      const Canvas dmg = imglib::load_png(argv[2]);
      const Canvas mask = imglib::load_png(argv[3]);
      if (cmd == "harmonic") {
        const int iters = (argc >= 6) ? std::stoi(argv[5]) : 160;
        imglib::save_png(inpaint::inpaint_harmonic_cg(dmg, mask, iters).image, argv[4]);
      } else if (cmd == "fmm") {
        const int rad = (argc >= 6) ? std::stoi(argv[5]) : 3;
        imglib::save_png(inpaint::inpaint_telea_fmm(dmg, mask, rad).image, argv[4]);
      } else {
        const int iters = (argc >= 6) ? std::stoi(argv[5]) : 60;
        imglib::save_png(
            inpaint::inpaint_hybrid_structure_cg(dmg, mask, 0.08f, 3, iters).image, argv[4]);
      }
      return 0;
    }
  } catch (const std::exception& ex) {
    std::cerr << "inpaint error: " << ex.what() << "\n";
    return 1;
  }

  print_usage();
  return 1;
}
