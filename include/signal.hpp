#ifndef __DSP_SIGNAL_HPP
#define __DSP_SIGNAL_HPP

// Standard Library Inclusions
#include <complex>
#include <algorithm>
#include <string_view>
#include <type_traits>
// External Library Inclusions
#include <Eigen/Dense>
#include <TempMeta/utility>
#include <TempMeta/numerical>
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
   * @brief Numerical precision of the underlying representation of the
   *        Signal samples.
   */
  using precision = typename TempMeta::numerical_precision<Datatype>::type;
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
  Signal(std::function<Datatype(size_t)> f, const size_t num_samples,
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
  Datatype operator[](const size_t &idx) const { return *(begin() + idx); }

  /**
   * @brief Indexing overload for setting samples in the Signal.
   * @param idx Sample index to return.
   * @retval Reference to the sample.
   */
  Datatype& operator[](const size_t &idx) { return *(begin() + idx); }

  /**
   * @brief Stream overload for printing Signal to an output stream.
   * @param os  The output stream to output to.
   * @param sig The Signal object to output.
   * @return std::ostream& The modified output stream.
   */
  friend std::ostream& operator<<(std::ostream& os, const Signal& sig) {

    if constexpr(Size == Eigen::Dynamic) {
      os << "Signal Data is dynamically allocated: ";
    } else {
      os << "Signal Data is statically allocated: ";
    }
    os << "Supports " << sig.size() << " samples.\n";

    os << "Signal Datatype: " << TempMeta::type_name<Datatype>() << '\n';
    os << "Fourier Datatype: " << TempMeta::type_name<std::complex<precision>>() << '\n';

    return os;

  }

private:
  container data_;
  uint32_t sample_rate_;

};

} // namespace Dsp

#endif /* #ifndef __DSP_SIGNAL_HPP */
