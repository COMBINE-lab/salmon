/*
This code is written by kerukuro and released into public domain.
*/

#ifndef DIGESTPP_DETAIL_TRAITS_HPP
#define DIGESTPP_DETAIL_TRAITS_HPP

namespace digestpp
{
namespace detail
{

template <typename T>
struct is_xof
{
	static const bool value = T::is_xof;
};

template <typename T>
struct is_byte
{
	static const bool value = std::is_same<T, char>::value ||
			std::is_same<T, signed char>::value ||
			std::is_same<T, unsigned char>::value;
};

} // namespace detail
} // namespace digestpp

#endif // DIGESTPP_DETAIL_TRAITS_HPP

