/**
 * @file test_signal.cpp
 * @author Salvatore Cardamone (sav.cardamone@gmail.com)
 * @brief Unit test for the Signal class.
 */
#define BOOST_TEST_DYN_LINK

// Standard Library Inclusions
#include <vector>
// External Library Inclusions
#include <boost/test/unit_test.hpp>
// Project Inclusions
#include "dsp/signal.hpp"
#include "dsp/overlap.hpp"

static std::vector<float> taps = { 0.1, 0.2, 0.3 };

/**
 * @brief Check the Valid mode of the convolution. This is the only mode currently
 *        supported with the Dsp::Signal; the user should pad the signal at construction
 *        if they wish to invoke other convolutional modes.
 */
BOOST_AUTO_TEST_CASE( Valid )
{

  Dsp::Signal<float> sig(
    [](const size_t idx) {return static_cast<float>(idx+1);}, 5, 5
  );
  std::vector<float> golden = { 1.0, 1.6, 2.2 };

  auto convolution = Dsp::convolve(sig, taps, Dsp::OverlapMode::Valid);
  BOOST_CHECK_EQUAL_COLLECTIONS(
    convolution.first, convolution.second, golden.begin(), golden.end()
  );

}

/**
 * @brief Dsp::Signal does not support insertion (yet?), so we require a
 *        runtime_error exception be thrown from the TempMeta::insert with
 *        the Full convolution mode.
 */
BOOST_AUTO_TEST_CASE( Full )
{

  Dsp::Signal<float> sig(
    [](const size_t idx) {return static_cast<float>(idx+1);}, 5, 5
  );

  BOOST_CHECK_THROW(
    Dsp::convolve(sig, taps, Dsp::OverlapMode::Full), std::runtime_error
  );

}

/**
 * @brief Dsp::Signal does not support insertion (yet?), so we require a
 *        runtime_error exception be thrown from the TempMeta::insert with
 *        the Same convolution mode.
 */
BOOST_AUTO_TEST_CASE( Same )
{

  Dsp::Signal<float> sig(
    [](const size_t idx) {return static_cast<float>(idx+1);}, 5, 5
  );

  BOOST_CHECK_THROW(
    Dsp::convolve(sig, taps, Dsp::OverlapMode::Same), std::runtime_error
  );

}