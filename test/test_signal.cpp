#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>
#include <vector>
#include <array>

#include "signal.hpp"

BOOST_AUTO_TEST_CASE(Construction)
{

    constexpr size_t static_size = 16;
    std::array<std::complex<double>,static_size> temp_static;
    std::vector<float> temp_dynamic(static_size);

    for (size_t ival=0; ival<static_size; ++ival) {
        temp_static[ival] = static_cast<float>(ival);
        temp_dynamic[ival] = static_cast<float>(ival);
    }

  Dsp::Signal<std::complex<double>,16> sig_static_fail(temp_static, 16, 16);
  Dsp::Signal<float,Eigen::Dynamic> sig_dynamic(temp_dynamic, 16, 16);
  std::cout << sig_static_fail;

}
