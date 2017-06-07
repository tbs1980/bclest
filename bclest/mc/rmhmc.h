#ifndef BCLEST_MC_RMHMC_H_
#define BCLEST_MC_RMHMC_H_

#include "chain.h"
#include <Eigen/CXX11/Tensor>
#include <Eigen/Dense>
#include <cstddef>
#include <functional>
#include <random>
#include <stdexcept>

namespace bclest {
namespace mc {

template <class Scalar> class ReimannHMC {
  static_assert(std::is_floating_point<RealScalar>::value,
                "ReimannHMC should have a floating point type");

public:
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, 1> Vector;
  typedef Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Matrix;
  typedef Eigen::Tensor<Scalar, 3> ThreeTensor;
  typedef typename std::function<Scalar(Vector const &)> LogPost;
  typedef typename std::function<Vector(Vector const &)> GradLogPost;
  typedef typename std::function<Matrix(Vector const &)> MetricTensorLogPost;
  typedef typename std::function<ThreeTensor(Vector const &)>
      DerivMetricTensorLogPost;
  typedef typename bclest::mc::Chain<Scalar> Chain;
  ReimannHMC(LogPost &log_post, GradLogPost &grad_log_post,
             MetricTensorLogPost &mtr_tnsr_log_post,
             DerivMetricTensorLogPost &deriv_mtr_tensr_log_post,
             std::size_t const num_dims, real_scalar_t const max_epsilon,
             std::size_t const max_leap_frog_steps,
             std::size_t const max_fixed_point_steps)
      : log_post_(log_post), grad_log_post_(grad_log_post),
        mtr_tnsr_log_post_(mtr_tnsr_log_post),
        deriv_mtr_tensr_log_post_(deriv_mtr_tensr_log_post),
        num_dims_(num_dims), max_epsilon_(max_epsilon),
        max_leap_frog_steps_(max_leap_frog_steps),
        max_fixed_point_steps_(max_fixed_point_steps), acc_rate_(0) {
    if (num_dims_ == std::size_t(0)) {
      throw std::invalid_argument("number of dimensions should be > 0.");
    }
    if (max_leap_frog_steps_ == std::size_t(0)) {
      throw std::invalid_argument("number of leapfrog steps should be > 0.");
    }
    if (max_fixed_point_steps_ == std::size_t(0)) {
      throw std::invalid_argument("number of fixed-point steps should be > 0.");
    }
    if (max_epsilon_ <= Scalar(0.) and max_epsilon_ > Scalar(1.)) {
      throw std::invalid_argument("0 <= max-epsiton < 1 not satisfied.");
    }
  }

  template <class RNG>
  Chain Run(std::size_t const num_samples, Vector const &start_point,
            RNG &rng) {
    typedef std::uniform_real_distribution<Scalar> UniformRealDistribution;
    typedef std::uniform_int_distribution<std::size_t> UniformIntDistribution;
    typedef std::normal_distribution<Scalar> NormalDistribution;
    if (num_samples == std::size_t(0)) {
      throw std::invalid_argument("num_samples should be > 0.");
    }
    if (start_point.rows() != num_dims_) {
      throw std::invalid_argument("rows in the start point != num_dims_.");
    }
    Chain rmhmc_chain(num_samples, num_dims_);
    UniformRealDistribution uni_real_dist;
    UniformIntDistribution uni_int_dist(std::size_t(1),
                                        max_leap_frog_steps_ + std::size_t(1));
    NormalDistribution norm_dist;
    Vector q_1(start_point);
    std::size_t num_accepted(0);
    std::size_t num_rejected(0);
    while (num_accepted < num_samples) {
      Vector q_0(q_1);
      Matrix G = mtr_tnsr_log_post_(q_0);
      Eigen::LLT<Matrix> LLTOfG(G);
      if (LLTOfG.info() != Eigen::Success) {
        throw std::runtime_error("Metric tensor not positive definite.");
      }
      Matrix chol_G = LLTOfG.matrixL();
      Scalar const step_size = max_epsilon_ * uni_real_dist(rng);
      std::size_t const num_leap_frog_steps = uni_int_dist(rng);
      Vector p_0(num_dims_);
      for (std::size_t i = 0; i < m_num_dims; ++i) {
        p_0(dim_i) = norm_dist(rng);
      }
      p_0 = p_0.transpose() * chol_G;
      Scalar const delta_h = stomer_verlet(
          num_leap_frog_steps, m_num_fixed_point_steps, step_size, ) {}
    }
  }

  ~ReimannHMC() {}

private:
  Scalar StomerVerlet(std::size_t const num_leap_frog_steps,
                      std::size_t const num_fixed_point_steps,
                      Scalar const step_size, LogPost &log_post,
                      GradLogPost &grad_log_post,
                      MetricTensorLogPost &mtr_tnsr_log_post,
                      DerivMetricTensorLogPost &drv_mtr_tnsr_log_post,
                      Vector &p_new, Vector &x_new) {
    Scalar const c_e = 1e-4;
    std::size_t const num_dims = p_new.rows();

    Matrix G = mtr_tnsr_log_post_(x_new); // FIXME this can be stored elsewhere
    Scalar det_G = compute_determinant<real_scalar_t>(G);
    Matrix chol_G <
  }

  LogPost &log_post_;
  GradLogPost &grad_log_post_;
  MetricTensorLogPost &mtr_tnsr_log_post_;
  DerivMetricTensorLogPost &deriv_mtr_tensr_log_post_;
  std::size_t const num_dims_;
  Scalar const max_epsilon_;
  std::size_t const max_leap_frog_steps_;
  std::size_t const max_fixed_point_steps_;
  RealScalar acc_rate_;
};

} // end namespace bclest
} // end namespace mc

#endif // BCLEST_MC_RMHMC_H_