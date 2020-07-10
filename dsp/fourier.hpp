#ifndef __DSP_FOURIER_HPP
#define __DSP_FOURIER_HPP

#include <complex>
#include <Eigen/Dense>
#include "dsp/signal.hpp"

namespace Dsp::Fourier {

/**
 * @brief Brute force Fourier transforming using the Vandermonde matrix.
 *        Quadratic scaling in signal length. Wouldn't advise using this for
 *        anything other than verification.
 * @todo If Size is a compile-time constant, make the Vandermonde matrix at
 *       compile-time.
 * @tparam Precision Numerical datatype to use as the template argument to std::complex.
 * @tparam Size      Size of the Vandermonde matrix.
 */
template <typename Precision = float, int Size = Eigen::Dynamic>
struct Vandermonde {

  /**
   * @brief Construct a new Vandermonde object.
   *
   * @param dft_length Runtime size of the Vandermonde matrix. Defaults to the statically
   */
  Vandermonde(const int& dft_length = Size)
    : dft_matrix_(dft_length, dft_length)  {

    // If the Vandermonde matrix is dynamically-allocated, the user must provide a valid
    // runtime size to instantiate it with
    if (Size == Eigen::Dynamic && dft_length == Size)
      throw std::runtime_error("Vandermonde matrix is dynamically-allocated; runtime argument is -1.");
    // If Size indicates the DFT matrix is to be statically-allocated and the runtime
    // dft_length argument exceeds this value, we can't use the DFT matrix to do a
    // Fourier transform
    if (Size > 0 && dft_length > Size)
      throw std::runtime_error("DFT length can't exceed statically-allocated size of Vandermonde matrix.");

    // The base argument of omega
    std::complex<Precision> arg(0, -2*M_PI/dft_length);

    // We only need the upper triangle of the Vandermonde matrix
    for (int i_row=0; i_row<dft_length; ++i_row) {
      for (int i_col=0; i_col<dft_length; ++i_col) {
        dft_matrix_(i_row,i_col) = std::exp(arg*static_cast<Precision>(i_row*i_col));
      }
    }

  }

  /**
   * @brief Apply the Fourier transform to a Signal. Sets the Signal's fourier member.
   * @tparam SignalDatatype Underlying datatype of the samples.
   * @tparam SignalSize     Length of the Signal.
   * @param sig The Signal we're Fourier transforming.
   * @return auto The Fourier coefficients of the Signal.
   */
  template <typename SignalDatatype, int SignalSize>
  auto apply(Signal<SignalDatatype, SignalSize>& sig) {

    sig.fourier() = dft_matrix_ * sig.data();
    return sig.fourier();

  }

  /**
   * @brief Apply the inverse Fourier transform to a Signal's Fourier coefficients.
   *
   *        Note that we don't assign the result to the Signal's data member. The Signal
   *        could be real-valued but the inverse Fourier transform will return a complex-
   *        valued signal, albeit the imaginary parts should all be zero.
   * @tparam SignalDatatype Underlying datatype of the samples.
   * @tparam SignalSize     Length of the Signal.
   * @param sig The Signal we're inverse-Fourier transforming.
   * @return auto The inverse Fourier transform of the Signals' Fourier coefficients.
   */
  template <typename SignalDatatype, int SignalSize>
  auto inverse(Signal<SignalDatatype, SignalSize>& sig) {

    return dft_matrix_.adjoint() * sig.fourier() / sig.size();

  }

private:
  Eigen::Matrix<std::complex<Precision>, Size, Size> dft_matrix_;

};

} // namespace Dsp::Fourier

#endif /* #ifndef __DSP_FOURIER_HPP */