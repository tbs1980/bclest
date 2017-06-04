#ifndef BCLEST_MC_CHAIN_H_
#define BCLEST_MC_CHAIN_H_

#include <cstddef>
#include <stdexcept>
#include <type_traits>

#include <Eigen/Dense>

namespace bclest {
namespace mc {

template <class RealScalar> class Chain {
  static_assert(std::is_floating_point<RealScalar>::value,
                "Chain should be a floating point type");
  typedef Eigen::Matrix<RealScalar, Eigen::Dynamic, Eigen::Dynamic> RealMatrix;
  typedef Eigen::Matrix<RealScalar, Eigen::Dynamic, 1> RealVector;

public:
  explicit Chain(std::size_t const num_samples, std::size_t const num_dims)
      : num_samples_(num_samples), num_dims_(num_dims),
        samples_(num_samples, num_dims), weights_(num_samples) {
    if (num_samples == std::size_t(0)) {
      throw std::length_error("num_samples should be > 0");
    }
    if (num_dims == std::size_t(0)) {
      throw std::length_error("num_dims should be > 0");
    }
  }

  ~Chain() {}

  inline std::size_t NumSamples() const { return num_samples_; }

  inline std::size_t NumDims() const { return num_dims_; }

  inline RealMatrix const &GetSamples() const { return samples_; }

  inline RealVector const &GetWeights() const { return weights_; }

  void SetSample(std::size_t const sample_id, RealVector const &sample,
                 RealScalar const weight) {
    if (sample_id >= num_samples_) {
      throw std::length_error("sample_id should be < num_samples_");
    }
    if (std::size_t(sample.rows()) != num_dims_) {
      throw std::length_error("rows in samples should be = num_dims_");
    }
    samples_.row(sample_id) = sample;
    weights_(sample_id) = weight;
  }

private:
  std::size_t num_samples_;
  std::size_t num_dims_;
  RealMatrix samples_;
  RealVector weights_;
};

} // end namespace bclest
} // end namespace mc

#endif // BCLEST_MC_CHAIN_H_