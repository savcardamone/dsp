/**
 * @file test_signal.cpp
 * @author Salvatore Cardamone (sav.cardamone@gmail.com)
 * @brief Unit test for the Signal class.
 */
#define BOOST_TEST_DYN_LINK

// Standard Library Inclusions
// External Library Inclusions
#include <boost/test/unit_test.hpp>
// Project Inclusions
#include "dsp/fourier.hpp"
#include "dsp/signal.hpp"

/**
 * @brief Check the brute force Fourier Transform works as expected.
 */
BOOST_AUTO_TEST_CASE( Vandermonde )
{

  // Make a the signal a single cosine
  Dsp::Signal<double> sig(
    [](const size_t idx){return std::cos(2*M_PI*static_cast<double>(idx)/8);}, 16, 32
  );
  Dsp::Fourier::Vandermonde<double,16> dft;

  std::cout << "Fourier Coefficients: " << dft.apply(sig) << std::endl;
  std::cout << "Inverse Fourier: " << dft.inverse(sig) << std::endl;

}
