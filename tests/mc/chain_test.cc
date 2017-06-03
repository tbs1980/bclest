#include <bclest/mc/chain.h>
#include <gtest/gtest.h>

TEST(ChainTest, ChainConstructor) {
  using namespace bclest::mc;
  std::size_t const num_samples(1000);
  std::size_t const num_dims(100);
  Chain<double> chain(num_samples, num_dims);
  ASSERT_EQ(chain.NumSamples(), num_samples);
  ASSERT_EQ(chain.NumDims(), num_dims);
}