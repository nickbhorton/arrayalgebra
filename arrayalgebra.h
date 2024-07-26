#ifndef ARRAY_ALGEBRA_HEADER_
#define ARRAY_ALGEBRA_HEADER_

#include <array>
#include <cmath>
#include <iostream>
#include <numeric>

namespace aa
{

typedef float Float;
typedef int Integer;
typedef unsigned int Uinteger;

typedef std::array<Float, 2> vec2;
typedef std::array<Float, 3> vec3;
typedef std::array<Float, 4> vec4;

typedef std::array<Integer, 2> ivec2;
typedef std::array<Integer, 3> ivec3;
typedef std::array<Integer, 4> ivec4;

typedef std::array<Uinteger, 2> uvec2;
typedef std::array<Uinteger, 3> uvec3;
typedef std::array<Uinteger, 4> uvec4;

typedef std::array<std::array<Float, 2>, 2> mat2;
typedef std::array<std::array<Float, 3>, 3> mat3;
typedef std::array<std::array<Float, 4>, 4> mat4;

typedef std::array<std::array<Integer, 2>, 2> imat2;
typedef std::array<std::array<Integer, 3>, 3> imat3;
typedef std::array<std::array<Integer, 4>, 4> imat4;

typedef std::array<std::array<Uinteger, 2>, 2> umat2;
typedef std::array<std::array<Uinteger, 3>, 3> umat3;
typedef std::array<std::array<Uinteger, 4>, 4> umat4;
} // namespace aa

// Generic size operations
template <typename T, std::size_t N>
constexpr auto
operator<<(std::ostream& os, std::array<T, N> const& a) -> std::ostream&
{
    for (size_t j = 0; j < N; j++) {
        os << a[j] << " ";
    }
    return os;
}

template <typename T, std::size_t N>
constexpr auto operator+(std::array<T, N> const& a1, std::array<T, N> const& a2)
    -> std::array<T, N>
{
    std::array<T, N> result{};
    for (size_t j = 0; j < N; j++) {
        result[j] = a1[j] + a2[j];
    }
    return result;
}

template <typename T, std::size_t N, std::size_t M>
constexpr auto operator+(
    std::array<std::array<T, M>, N> const& m1,
    std::array<std::array<T, M>, N> const& m2
) -> std::array<std::array<T, M>, N>
{
    std::array<std::array<T, M>, N> result{};
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < M; j++) {
            result[i][j] = m1[i][j] + m2[i][j];
        }
    }
    return result;
}

template <typename T, std::size_t N>
constexpr auto operator-(std::array<T, N> const& a1, std::array<T, N> const& a2)
    -> std::array<T, N>
{
    std::array<T, N> result{};
    for (size_t j = 0; j < N; j++) {
        result[j] = a1[j] - a2[j];
    }
    return result;
}

template <typename T, std::size_t N, std::size_t M>
constexpr auto operator-(
    std::array<std::array<T, M>, N> const& m1,
    std::array<std::array<T, M>, N> const& m2
) -> std::array<std::array<T, M>, N>
{
    std::array<std::array<T, M>, N> result{};
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < M; j++) {
            result[i][j] = m1[i][j] - m2[i][j];
        }
    }
    return result;
}

template <typename T, std::size_t N>
constexpr auto operator-(std::array<T, N> const& a1) -> std::array<T, N>
{
    std::array<T, N> result{};
    for (size_t j = 0; j < N; j++) {
        result[j] = -a1[j];
    }
    return result;
}

template <typename T, std::size_t N, std::size_t M>
constexpr auto operator-(std::array<std::array<T, M>, N> const& m
) -> std::array<std::array<T, M>, N>
{
    std::array<std::array<T, M>, N> result{};
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < M; j++) {
            result[i][j] = -m[i][j];
        }
    }
    return result;
}

template <typename T, std::size_t N>
constexpr auto
operator*(std::array<T, N> const& a, T const v) -> std::array<T, N>
{
    std::array<T, N> result{};
    for (size_t j = 0; j < N; j++) {
        result[j] = v * a[j];
    }
    return result;
}

template <typename T, std::size_t N, std::size_t M>
constexpr auto operator*(std::array<std::array<T, M>, N> const& m, T const& v)
    -> std::array<std::array<T, M>, N>
{
    std::array<std::array<T, M>, N> result{};
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < M; j++) {
            result[i][j] = v * m[i][j];
        }
    }
    return result;
}

template <typename T, std::size_t N>
constexpr auto
operator*(T const v, std::array<T, N> const& a) -> std::array<T, N>
{
    std::array<T, N> result{};
    for (size_t j = 0; j < N; j++) {
        result[j] = v * a[j];
    }
    return result;
}

template <typename T, std::size_t N, std::size_t M>
constexpr auto operator*(T const& v, std::array<std::array<T, M>, N> const& m)
    -> std::array<std::array<T, M>, N>
{
    std::array<std::array<T, M>, N> result{};
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < M; j++) {
            result[i][j] = v * m[i][j];
        }
    }
    return result;
}

template <typename T, std::size_t N>
constexpr auto operator*(
    std::array<std::array<T, N>, N> const& mat,
    std::array<T, N> const& vec
) -> std::array<T, N>
{
    std::array<T, N> result{};
    for (size_t i = 0; i < N; i++) {
        T sum = static_cast<T>(0);
        for (size_t j = 0; j < N; j++) {
            sum += mat[i][j] * vec[j];
        }
        result[i] = sum;
    }
    return result;
}

template <typename T, std::size_t N, std::size_t M>
constexpr auto operator*(
    std::array<std::array<T, N>, M> const& m1,
    std::array<std::array<T, M>, N> const& m2
) -> std::array<std::array<T, M>, M>
{
    std::array<std::array<T, M>, M> result{};
    for (size_t i = 0; i < M; i++) {
        for (size_t j = 0; j < M; j++) {
            for (size_t k = 0; k < N; k++) {
                result[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }
    return result;
}

template <typename T, std::size_t N>
constexpr auto dot(std::array<T, N> const& a1, std::array<T, N> const& a2) -> T
{
    return std::inner_product(
        a1.begin(),
        a1.end(),
        a2.begin(),
        static_cast<T>(0)
    );
}

template <typename T, std::size_t N>
constexpr auto transpose(std::array<std::array<T, N>, N> const& a
) -> std::array<std::array<T, N>, N>
{
    std::array<std::array<T, N>, N> result;
    for (std::size_t i = 0; i < N; i++) {
        for (std::size_t j = 0; j < N; j++) {
            result[i][j] = a[j][i];
        }
    }
    return result;
}

template <typename T>
constexpr auto cross(std::array<T, 3> const& v1, std::array<T, 3> const& v2)
    -> std::array<T, 3>
{
    std::array<T, 3> result{};
    result[0] = v1[1] * v2[2] - v1[2] * v2[1];
    result[1] = v1[2] * v2[0] - v1[0] * v2[2];
    result[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return result;
}

template <typename T, std::size_t N>
constexpr auto magnitude(std::array<T, N> const& v) -> T
{
    T acc{};
    for (size_t i = 0; i < N; i++) {
        acc += v[i] * v[i];
    }
    return std::sqrt(acc);
}

template <typename T, std::size_t N>
constexpr auto normalize(std::array<T, N> const& v) -> std::array<T, N>
{
    std::array<T, N> result{};
    double const length = magnitude(v);
    for (size_t i = 0; i < N; i++) {
        result[i] = v / length;
    }
    return result;
}

#endif
