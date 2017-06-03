#ifndef BCLEST_MC_CHAIN_H_
#define BCLEST_MC_CHAIN_H_

#include <cstddef>
#include <stdexcept>
#include <type_traits>

namespace bclest {
namespace mc {

template <class RealScalar> class Chain {
  static_assert(std::is_floating_point<RealScalar>::value,
                "Chain should be a floating point type");

public:
  explicit Chain(std::size_t const num_samples, std::size_t const num_dims)
      : num_samples_(num_samples), num_dims_(num_dims) {
    if (num_samples == std::size_t(0)) {
      throw std::logic_error("num_samples should be > 0");
    }
    if (num_dims == std::size_t(0)) {
      throw std::logic_error("num_dims should be > 0");
    }
  }

  ~Chain() {}

  inline std::size_t NumSamples() const { return num_samples_; }

  inline std::size_t NumDims() const { return num_dims_; }

private:
  std::size_t num_samples_;
  std::size_t num_dims_;
};

} // end namespace bclest
} // end namespace mc

#endif // BCLEST_MC_CHAIN_H_