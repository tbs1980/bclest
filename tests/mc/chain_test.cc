#include <bclest/mc/chain.h>
#include <gtest/gtest.h>
#include <stdexcept>

#include <Eigen/Dense>

TEST(ChainTest, ChainConstructor) {
  using namespace bclest::mc;
  std::size_t const num_samples(1000);
  std::size_t const num_dims(100);
  Chain<double> chain(num_samples, num_dims);
  ASSERT_EQ(chain.NumSamples(), num_samples);
  ASSERT_EQ(chain.NumDims(), num_dims);
}

TEST(ChainTest, ChainConstructorNumSamplesZero) {
  using namespace bclest::mc;
  std::size_t const num_samples(0);
  std::size_t const num_dims(100);
  try {
    Chain<double> chain(num_samples, num_dims);
    // if the line below executes, the test should fail.
    // as we expect an exception to be thrown
    ASSERT_NE(chain.NumSamples(), num_samples);
  } catch (std::logic_error &lgerr) {
  }
}

TEST(ChainTest, ChainConstructorNumDimesZero) {
  using namespace bclest::mc;
  std::size_t const num_samples(1000);
  std::size_t const num_dims(0);
  try {
    Chain<double> chain(num_samples, num_dims);
    // if the line below executes, the test should fail.
    // as we expect an exception to be thrown
    ASSERT_NE(chain.NumDims(), num_dims);
  } catch (std::length_error &lgerr) {
  }
}

TEST(ChainTest, ChainSetSample) {
  using namespace bclest::mc;
  using namespace Eigen;
  std::size_t const num_samples(1000);
  std::size_t const num_dims(10);
  Chain<double> chain(num_samples, num_dims);
  VectorXd const sample = VectorXd::Random(num_dims);
  std::size_t const sample_id(10);
  double const weight(0.5);
  chain.SetSample(sample_id, sample, weight);
  for (std::size_t i = 0; i < num_dims; ++i) {
    ASSERT_EQ(chain.GetSamples().row(sample_id)(i), sample(i));
  }
  ASSERT_EQ(chain.GetWeights()(sample_id), weight);
}

TEST(ChainTest, ChainToCSV) {
  using namespace bclest::mc;
  using namespace Eigen;
  std::size_t const num_samples(1000);
  std::size_t const num_dims(10);
  Chain<double> chain(num_samples, num_dims);
  for (std::size_t i = 0; i < num_samples; ++i) {
    double const weight(0.5 * i);
    VectorXd const sample = VectorXd::Random(num_dims);
    chain.SetSample(i, sample, weight);
  }
  std::string chain_file_name("/tmp/chain_test.csv");
  chain.ToCSV(chain_file_name);
}