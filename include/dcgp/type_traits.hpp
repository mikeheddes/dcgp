#ifndef DCGP_TYPE_TRAITS_H
#define DCGP_TYPE_TRAITS_H

/// Type is a gdual
/**
 * Checks whether T is a gdual type. Provides the member constant value which is
 * equal to true, if T is the type gdual<U> for any U.
 *
 * \tparam T a type to check
 */

template <typename T>
struct is_gdual : std::false_type {
};

#ifndef __EMSCRIPTEN__
template <typename T>
struct is_gdual<audi::gdual<T>> : std::true_type {
};
#endif // __EMSCRIPTEN__

#endif // DCGP_TYPE_TRAITS_H
