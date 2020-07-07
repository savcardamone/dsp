/**
 * @file signal.hpp
 * @author Salvatore Cardamone (sav.cardamone@gmail.com)
 * @brief Wrapper class to hold a discretely sampled signal.
 */
#ifndef __DSP_SIGNAL_HPP
#define __DSP_SIGNAL_HPP

// Standard Library Inclusions
#include <complex>
#include <ostream>
#include <algorithm>
#include <string_view>
#include <type_traits>
// External Library Inclusions
#include <Eigen/Dense>
#include "tempmeta/utility.hpp"
#include "tempmeta/numerical.hpp"
// Project Inclusions


namespace Dsp {

/**
 * @brief
 * @tparam Datatype Underlying datatype of the signal. Defaults to complex real
 *                  (single precision floating point).
 * @tparam Size     Number of samples in the signal. If this is known at compile
 *                  time we can use a stattically-allocated buffer to store the
 *                  signal, else default to dynamic allocation.
   */
template <typename Datatype = std::complex<float>, int Size = Eigen::Dynamic>
struct Signal {

  /**
   * @brief The type of the contents of the container.
   */
  using value_type = Datatype;
  /**
   * @brief Numerical precision of the underlying representation of the
   *        Signal samples.
   */
  using working_precision = typename TempMeta::numerical_precision<Datatype>::type;
  /**
   * @brief Underlying container which stores the signal samples.
   *        1D Eigen matrix so we can do linear algebra with it.
   */
  using container = typename Eigen::Matrix<Datatype, Size, 1>;
  /**
   * @brief Iterator type so we can use the Signal as an iterable container.
   */
  using iterator = typename container::iterator;
  /**
   * @brief Reverse iterator type so we can use the Signal as an iterable
   *        container.
   */
  using reverse_iterator = typename std::reverse_iterator<iterator>;

  /**
   * @brief Signal class constructor. Used for dynamically-allocated signals.
   *        Initialise the samples from a pre-existing source (array or
   *        container).
   *
   *        Note that Eigen guarantees a uniform API for Matrix creation, so
   *        even if Signal::container is statically allocated, we can still
   *        call its constructor in the initialiser list without a
   *        compilation error.
   * @todo  Should std::move the resources of the initial container rather
   *        than duplicate contents.
   * @tparam Init        Container of the samples we're initialising with.
   * @param  c           Samples to initialise the signal with.
   * @param  num_samples Number of samples the signal should be initialised
   *                     with.
   * @param  sample_rate Acquisition sample rate, in units of Samples / Unit
   *                     Time.
   */
  template <typename Init>
  Signal(Init& c, const size_t num_samples, const uint32_t sample_rate)
    : data_(num_samples), sample_rate_(sample_rate) {

    // Assign the signal samples
    for (size_t isample=0; isample<num_samples; ++isample) {
      data_[isample] = c[isample];
    }

  }

  /**
   * @brief Signal class constructor. Initialise the samples from a function.
   * @param f           Function which takes a sample index and returns the
   *                    sample value.
   * @param  num_samples Number of samples the signal should be initialised
   *                     with.
   * @param  sample_rate Acquisition sample rate, in units of Samples / Unit
   *                     Time.
   */
  Signal(std::function<value_type(size_t)> f, const size_t num_samples,
         const uint32_t sample_rate)
    : data_(num_samples), sample_rate_(sample_rate) {

    for (size_t isample=0; isample<num_samples; ++isample) {
      data_[isample] = f(isample);
    }

  }

  /**
   * @brief Return the size of the underlying signal container.
   * @retval size_t The container size.
   */
  size_t size() const { return data_.size(); }

  /**
   * @brief Resize the Signal to hold a different number of samples.
   * @param num_samples The number of samples the Signal should now hold.
   */
  void resize(const size_t num_samples) {

    if constexpr (Size == Eigen::Dynamic) {
      data_.conservativeResize(num_samples, 1);
    } else {
      throw std::runtime_error("The Signal is of fixed size: can't resize.");
    }

  }

  /**
   * @brief Getter for the sample rate of the Signal in units of Samples / Unit
   *        Time.
   * @return uint32_t Acquisition sample rate, in units of Samples / Unit Time.
   */
  uint32_t sample_rate() const { return sample_rate_; }

  /**
   * @brief Return an iterator to the beginning of the signal.
   * @retval begin iterator.
   */
  iterator begin() { return data_.begin(); }

  /**
   * @brief Return an iterator to the end of the signal.
   * @retval end iterator.
   */
  iterator end() { return data_.end(); }

  /**
   * @brief Return a reverse iterator to the end of the signal.
   * @retval rbegin iterator.
   */
  reverse_iterator rbegin() { return std::make_reverse_iterator(end()); }

  /**
   * @brief Return a reverse iterator to the beginning of the signal.
   * @retval rbegin iterator.
   */
  reverse_iterator rend() { return std::make_reverse_iterator(begin()); }

  /**
   * @brief Const indexing overload for getting samples from the Signal.
   * @param idx Sample index to return.
   * @retval Value of the sample.
   */
  value_type operator[](const size_t &idx) const { return *(begin() + idx); }

  /**
   * @brief Indexing overload for setting samples in the Signal.
   * @param idx Sample index to return.
   * @retval Reference to the sample.
   */
  value_type& operator[](const size_t &idx) { return *(begin() + idx); }

  /**
   * @brief Compute the frequency resolution of the Signal in units of inverse
   *        Unit Time.
   * @return working_precision The frequency resolution of the signal.
   */
  working_precision resolution() const {
    return static_cast<working_precision>(sample_rate()) / size();
  }

  /**
   * @brief Stream overload for printing Signal to an output stream.
   * @param os  The output stream to output to.
   * @param sig The Signal object to output.
   * @return std::ostream& The modified output stream.
   */
  friend std::ostream& operator<<(std::ostream& os, Signal& sig) {

    // If the output stream is std::cout, then we print some
    // condensed information about the signal
    if (&os == &std::cout) {
      sig.print_to_stdout(os);
    // Otherwise it's to file, so we'll print everything we can about the
    // object in xml format
    } else {
      sig.print_to_xml(os);
    }

    return os;

  }

private:
  container data_;
  uint32_t sample_rate_;

  /**
   * @brief Print the object in xml format to an output stream.
   * @param os The output stream we print to.
   */
  void print_to_xml(std::ostream& os) {

    // Set up some precision and numerical widths for the output stream
    os.precision(8);
    os << std::fixed << "<?xml version=\"1.0\"?>" << "\n\n";

    // Generic information about the Signal
    std::string signal_type = TempMeta::is_complex<value_type>() ? "complex" : "real";
    os << "<Signal num_samples=\"" << size() << "\" "
       << "sample_rate=\"" << sample_rate() << "\" "
       << "type=\"" << signal_type << "\""
       << ">\n";

    // Print the Signal samples and the associated sample times
    os << "\t<Samples>\n";
    for (uint32_t isample=0; isample<size(); ++isample) {

      auto val = (*this)[isample];
      os << "\t\t<Sample t=\""
         << static_cast<working_precision>(isample) / sample_rate()
         << "\"> ";

      // Got to handle whether the Signal comprises complex or real samples
      // and print accordingly
      if constexpr (TempMeta::is_complex<value_type>()) {
        os << std::real(val) << "," << std::imag(val);
      } else {
        os << val;
      }
      os << " </Sample>\n";

    }
    os << "\t</Samples>\n";

    // Finalise the xml
    os << "</Signal>\n";

  }

  /**
   * @brief Print the object in condensed format to std::cout.
   * @param os The output stream we print to.
   */
  void print_to_stdout(std::ostream& os) {

    if constexpr(Size == Eigen::Dynamic) {
      os << "Signal Data is dynamically allocated: ";
    } else {
      os << "Signal Data is statically allocated: ";
    }
    os << "Supports " << size() << " samples.\n";

    os << "Signal Datatype: " << TempMeta::type_name<value_type>() << '\n';
    os << "Fourier Datatype: " << TempMeta::type_name<std::complex<working_precision>>() << '\n';

  }

};

} // namespace Dsp

#endif /* #ifndef __DSP_SIGNAL_HPP */
