#ifndef INPAINT_SOLVER_HPP
#define INPAINT_SOLVER_HPP

#include "inpaint/deconv.hpp"
#include "inpaint/operators.hpp"
#include "inpaint/poisson.hpp"

namespace inpaint {

struct InpaintResult {
  Canvas image;
  Canvas arrival_time;
  std::vector<double> residual_norms;
};

InpaintResult inpaint_harmonic_cg(const Canvas& corrupted,
                                  const Canvas& mask,
                                  int max_iters = 200,
                                  double tol = 1e-6);

InpaintResult inpaint_telea_fmm(const Canvas& corrupted,
                                const Canvas& mask,
                                int radius = 3);

InpaintResult inpaint_hybrid_structure_cg(const Canvas& corrupted,
                                          const Canvas& mask,
                                          float kappa = 0.08f,
                                          int fmm_radius = 3,
                                          int cg_iters = 80,
                                          double tol = 1e-6);

} // namespace inpaint

#endif // INPAINT_SOLVER_HPP
