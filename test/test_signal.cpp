/**
 * @file test_signal.cpp
 * @author Salvatore Cardamone (sav.cardamone@gmail.com)
 * @brief Unit test for the Signal class.
 */
#define BOOST_TEST_DYN_LINK

// Standard Library Inclusions
#include <array>
#include <vector>
// External Library Inclusions
#include <boost/test/unit_test.hpp>
// Project Inclusions
#include "dsp/signal.hpp"

/**
 * @brief Check the constructors work as expected.
 */
BOOST_AUTO_TEST_CASE( Construction )
{

  // Check the Signal can be initialised with a statically allocated container
  std::array<double,16> temp_static;
  for (size_t ival=0; ival<16; ++ival)
      temp_static[ival] = static_cast<float>(ival);
  Dsp::Signal<double,16> sig_static(temp_static, 16, 16);
  BOOST_TEST(sig_static.size() == temp_static.size());
  BOOST_CHECK_EQUAL_COLLECTIONS(
    sig_static.begin(), sig_static.end(), temp_static.begin(), temp_static.end()
  );

  // Check the Signal can be initialised with a dynamically allocated container
  std::vector<double> temp_dynamic(16);
  for (size_t ival=0; ival<16; ++ival)
      temp_dynamic[ival] = static_cast<float>(ival);
  Dsp::Signal<double> sig_dynamic(temp_dynamic, 16, 16);
  BOOST_TEST(sig_dynamic.size() == temp_dynamic.size());
  BOOST_CHECK_EQUAL_COLLECTIONS(
    sig_dynamic.begin(), sig_dynamic.end(), temp_dynamic.begin(), temp_dynamic.end()
  );

  // Check the Signal can be initialised with a function
  Dsp::Signal<double> sig_function([](const size_t idx){return static_cast<double>(idx);}, 16, 16);
  BOOST_TEST(sig_function.size() == temp_dynamic.size());
  BOOST_CHECK_EQUAL_COLLECTIONS(
    sig_function.begin(), sig_function.end(), temp_dynamic.begin(), temp_dynamic.end()
  );

}
