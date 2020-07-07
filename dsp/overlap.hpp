/**
 * @file overlap.hpp
 * @author Salvatore Cardamone (sav.cardamone@gmail.com)
 * @brief Overlapping functions -- convolution and correlation of signals.
 *
 *        Convolution and correlation, while not related, have far too similar a
 *        kernel to warrant factoring into separate implementations; convolution can
 *        be done with a reverse iterator over the signal while correlation can be done
 *        with a forward iterator. We'll let the compiler work out how best to implement.
 *
 *        Note that all of these functions are destructive -- they'll overwrite the
 *        LHS operand of the convolution/correlation operator.
 */
#ifndef __DSP_OVERLAP_HPP
#define __DSP_OVERLAP_HPP

// Standard Library Inclusions
#include <cstdlib>
#include <utility>
#include <cassert>
#include <iterator>
// External Library Inclusions
#include "tempmeta/container.hpp"
#include "tempmeta/numerical.hpp"
// Project Inclusions


namespace Dsp {

/**
 * @brief Enumeration of the modes supported by the overlapping.
 *
 *        For Full and Same modes, padding will be added to the signal
 *        to accommodate filter delay. Valid is the only mode which
 *        won't resize the signal container. It's recommended to
 *        initialise the signal container with the appropriate padding
 *        and use the Valid mode if performance is required. If Full
 *        or Same are chosen, the signal container must support an
 *        insert method.
 *
 * @param Full Include the filter delay in the overlap.
 * @param Valid No boundary effects from the filter delay.
 * @param Same Bounday effect at the end.
 */
enum class OverlapMode
{
    /**
     * Length of the result is signal_length + taps_length - 1, where
     * signal_length is the signal length before padding.
     *
     *        Signal: [p p 0 1 2 3 4 p p]
     *        Iter 1: [a b c]             ---> [* p 0 1 2 3 4 p p]
     *        Iter 2:   [a b c]           ---> [* * 0 1 2 3 4 p p]
     *           ...
     *        Iter 7:             [a b c] ---> [* * * * * * * p p]
     *                                          |           |
     *        Return:                         begin        end
     */
    Full,
    /**
     * Length of the result is signal_length - taps_length + 1.
     *
     *        Signal: [0 1 2 3 4]
     *        Iter 1: [a b c]     ---> [* 1 2 3 4]
     *        Iter 2:   [a b c]   ---> [* * 2 3 4]
     *        Iter 3:     [a b c] ---> [* * * 3 4]
     *                                  |   |
     *        Return:                 begin end
     */
    Valid,
    /**
     * Length of the result is signal_length.
     *
     *        Signal: [0 1 2 3 4 p p]
     *        Iter 1: [a b c]         ---> [* 1 2 3 4 p p]
     *        Iter 2:   [a b c]       ---> [* * 2 3 4 p p]
     *           ...
     *        Iter 5:         [a b c] ---> [* * * * * p p]
     *                                      |       |
     *        Return:                     begin    end
     */
    Same
};

/**
 * @brief  Main compute kernel for the overlapping of a signal with a filter.
 * @param  sig_start Iterator to the "start" of the signal.
 * @param  sig_end Iterator to the "end" of the filter.
 * @param  filter_start Iterator to the "start" of the filter taps.
 * @param  filter_end Iterator to the "end" of the filter taps.
 * @param  overlap The forward iterator for the destination of the result.
 * @retval None.
 */
template <class Iter, class Jter>
void overlap_kernel(Iter sig_start, Iter sig_end,
                    const Jter filter_start, const Jter filter_end,
                    Iter overlap)
{

    // Iterate over our signal between the argument bounds
    for (auto sig = sig_start; sig != sig_end; ++sig)
    {

        // Zero our accumulator
        typename std::iterator_traits<Iter>::value_type temp = 0;
        // Duplicate our signal iterator so we can do the overlap
        auto sig_val = sig;

        // Accumulate our overlap over signal and filter
        for (auto filter_val = filter_start; filter_val != filter_end; ++filter_val)
        {
            temp += *filter_val * *sig_val++;
        }

        // Dump the accumulator to the relevant bit of destination and increment
        *overlap++ = temp;
    }
}

/**
 * @brief Perform the convolution, a * b, and store the result in a.
 * @param  a The LHS of the convolution operator. Note len(a) >= len(b).
 * @param  b The RHS of the convolution operator.
 * @param  mode Overlap mode.
 * @retval Pair of iterators to the beginning and end of the convolution
 *         as per the overlap mode. The result overwrites a.
 */
template <class ContainerA, class ContainerB>
auto convolve(ContainerA& a, const ContainerB& b, const OverlapMode& mode)
    -> std::pair<typename ContainerA::iterator, typename ContainerA::iterator>
{

    // Since we're doing the convolution in place for ContainerA, we can't
    // store the result if ContainerA has a smaller datatype than ContainerB
    // e.g. if ContainerB is a std::complex<T> and ContainerA is a T.
    static_assert(sizeof(typename ContainerA::value_type) >= sizeof(typename ContainerB::value_type));
    // Seems pretty pointless to support a filter that's longer than the signal
    // May as well remove this kind of edge case
    assert(a.size() > b.size());

    auto filter_delay_length = b.size() - 1;

    switch (mode)
    {

    case OverlapMode::Full:
        // Pad a with zeros at the beginning to hold the prologue of the convolution
        TempMeta::insert(a, std::begin(a), filter_delay_length, 0);
        // Pad a with zeros at the end to hold the epilogue of the convolution
        TempMeta::insert(a, std::end(a), filter_delay_length, 0);
        break;

    case OverlapMode::Valid:
        // We have nothing to do here
        break;

    case OverlapMode::Same:
        // Pad a with zeros at the end to hold the epilogue of the convolution
        TempMeta::insert(a, std::begin(a), filter_delay_length, 0);
        break;

    default:
        // Return the signal unmodified if no supported mode was chosen
        return std::make_pair(a.begin(), a.end());

    }

    // Execute the convolutional kernel
    //overlap_kernel(a.rbegin(), a.rend() - filter_delay_length, b.begin(), b.end(), a.rbegin());
    overlap_kernel(a.rbegin(), a.rend() - filter_delay_length, b.begin(), b.end(), a.rbegin());
    // Return iterators to the valid convolutional result
    return std::make_pair(a.begin() + filter_delay_length, a.end());

}

/**
 * @brief Perform the correlate, a x b, and store the result in a.
 * @param  a The LHS of the correlation operator. Note len(a) >= len(b).
 * @param  b The RHS of the correlation operator.
 * @param  mode Overlap mode.
 * @retval Pair of iterators to the beginning and end of the correlation
 *         as per the overlap mode. The result overwrites a.
 */
template <class ContainerA, class ContainerB>
auto correlate(ContainerA& a, const ContainerB& b, const OverlapMode& mode)
    -> std::pair<typename ContainerA::iterator, typename ContainerA::iterator>
{

    // Since we're doing the correlation in place for ContainerA, we can't
    // store the result if ContainerA has a smaller datatype than ContainerB
    // e.g. if ContainerB is a std::complex<T> and ContainerA is a T.
    static_assert(sizeof(ContainerA::value_type) >= sizeof(ContainerB::value_type));
    // Seems pretty pointless to support a filter that's longer than the signal
    // May as well remove this kind of edge case
    assert(a.size() > b.size());

    // Resolve at compile-time: if the RHS of the correlation operator is complex,
    // take the complex-conjugate in-place
    if constexpr(TempMeta::is_complex<typename ContainerB::value_type>())
    {
        for (auto &val : b) val = std::conj(val);
    }

    auto filter_delay_length = b.size() - 1;

    switch (mode)
    {

    case OverlapMode::Full:
        // Pad a with zeros at the beginning to hold the prologue of the correlation
        TempMeta::insert(a, std::begin(a), filter_delay_length, 0);
        // Pad a with zeros at the end to hold the epilogue of the correlation
        TempMeta::insert(a, std::end(a), filter_delay_length, 0);
        break;

    case OverlapMode::Valid:
        // We have nothing to do here
        break;

    case OverlapMode::Same:
        // Pad a with zeros at the end to hold the epilogue of the correlation
        TempMeta::insert(a, std::begin(a), filter_delay_length, 0);
        break;

    default:
        // Return the signal unmodified if no supported mode was chosen
        return std::make_pair(a.begin(), a.end());

    }

    // Execute the correlation kernel
    overlap_kernel(a.begin(), a.end() - filter_delay_length, b.begin(), b.end(), a.begin());
    // Return iterators to the valid correlation result
    return std::make_pair(a.begin(), a.end() - filter_delay_length);

}

} // namespace Dsp

#endif /* #ifndef __DSP_OVERLAP_HPP */