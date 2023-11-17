#ifndef INPAINT_POISSON_HPP
#define INPAINT_POISSON_HPP

#include "inpaint/operators.hpp"
#include <vector>

namespace inpaint {

struct PoissonResult {
  Canvas image;
  std::vector<double> residual_norms;
};

Canvas threshold_mask(const Canvas& raw_mask, float threshold = 0.5f);

Canvas naive_composite(const Canvas& background,
                       const Canvas& foreground,
                       const Canvas& mask);

PoissonResult poisson_steepest_descent(const Canvas& background,
                                       const Canvas& foreground,
                                       const Canvas& mask,
                                       int max_iters = 200,
                                       double tol = 1e-6);

PoissonResult poisson_cg(const Canvas& background,
                         const Canvas& foreground,
                         const Canvas& mask,
                         int max_iters = 200,
                         double tol = 1e-6,
                         const Canvas* initial_guess = nullptr);

PoissonResult poisson_mixed_cg(const Canvas& background,
                               const Canvas& foreground,
                               const Canvas& mask,
                               int max_iters = 200,
                               double tol = 1e-6);

} // namespace inpaint

#endif // INPAINT_POISSON_HPP
