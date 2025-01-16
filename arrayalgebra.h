#ifndef ARRAY_ALGEBRA_HEADER_
#define ARRAY_ALGEBRA_HEADER_

#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <numeric>
#include <tuple>

namespace aa
{

typedef float Float;
typedef int Integer;
typedef unsigned int Uinteger;

typedef std::array<Float, 2> vec2;
typedef std::array<Float, 3> vec3;
typedef std::array<Float, 4> vec4;

typedef std::array<double, 2> dvec2;
typedef std::array<double, 3> dvec3;
typedef std::array<double, 4> dvec4;

typedef std::array<Integer, 2> ivec2;
typedef std::array<Integer, 3> ivec3;
typedef std::array<Integer, 4> ivec4;

typedef std::array<Uinteger, 2> uvec2;
typedef std::array<Uinteger, 3> uvec3;
typedef std::array<Uinteger, 4> uvec4;

typedef std::array<uint8_t, 2> u8vec2;
typedef std::array<uint8_t, 3> u8vec3;
typedef std::array<uint8_t, 4> u8vec4;

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

// Printing function
template <typename T, std::size_t N>
constexpr auto
operator<<(std::ostream& os, std::array<T, N> const& a) -> std::ostream&
{
    for (size_t j = 0; j < N; j++) {
        os << a[j] << " ";
    }
    return os;
}

template <typename T, std::size_t N, size_t M>
constexpr auto
operator<<(std::ostream& os, std::array<std::array<T, N>, M> const& a)
    -> std::ostream&
{
    for (size_t i = 0; i < M; i++) {
        for (size_t j = 0; j < N; j++) {
            os << a[i][j] << " ";
        }
        std::cout << "\n";
    }
    return os;
}

// v1 + v2
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

// m1 + m2
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

// v1 - v2
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

// m1 - m2
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

// -v
template <typename T, std::size_t N>
constexpr auto operator-(std::array<T, N> const& a1) -> std::array<T, N>
{
    std::array<T, N> result{};
    for (size_t j = 0; j < N; j++) {
        result[j] = -a1[j];
    }
    return result;
}

// -m
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

// v * c
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

// m * c
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

// c * v
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

// c * m
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

// m * v
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

// m1 * m2
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
constexpr auto operator+=(std::array<T, N>& v1, std::array<T, N> const& v2)
    -> std::array<T, N>&
{
    for (size_t i = 0; i < N; i++) {
        v1[i] += v2[i];
    }
    return v1;
}

template <typename T, std::size_t N>
constexpr auto operator-=(std::array<T, N>& v1, std::array<T, N> const& v2)
    -> std::array<T, N>&
{
    for (size_t i = 0; i < N; i++) {
        v1[i] -= v2[i];
    }
    return v1;
}

namespace aa
{
// v1 dot v2
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
constexpr T distance(std::array<T, N> const& a1, std::array<T, N> const& a2)
{
    T sum{};
    for (size_t i = 0; i < N; i++) {
        sum += std::pow(
            static_cast<aa::Float>(a1[i] - a2[i]),
            static_cast<aa::Float>(2.0)
        );
    }
    return std::sqrt(static_cast<aa::Float>(sum));
}

// transpose m
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

template <typename T, std::size_t N>
constexpr auto identity() -> std::array<std::array<T, N>, N>
{
    std::array<std::array<T, N>, N> result{};
    for (std::size_t i = 0; i < N; i++) {
        result[i][i] = static_cast<T>(1);
    }
    return result;
}

// flatten m
template <typename T, std::size_t N>
constexpr auto flatten(std::array<std::array<T, N>, N> const& a
) -> std::array<T, N * N>
{
    std::array<T, N * N> result;
    for (std::size_t i = 0; i < N; i++) {
        for (std::size_t j = 0; j < N; j++) {
            result[i * N + j] = a[j][i];
        }
    }
    return result;
}

// v1 cross v2
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

// magnitude v
template <typename T, std::size_t N>
constexpr auto magnitude(std::array<T, N> const& v) -> T
{
    T acc{};
    for (size_t i = 0; i < N; i++) {
        acc += v[i] * v[i];
    }
    return std::sqrt(acc);
}

// normalize v
template <typename T, std::size_t N>
constexpr auto normalize(std::array<T, N> const& v) -> std::array<T, N>
{
    std::array<T, N> result{};
    double const length = magnitude(v);
    for (size_t i = 0; i < N; i++) {
        result[i] = v[i] * (static_cast<T>(1) / length);
    }
    return result;
}

// Tranformations needed for graphics
template <typename T>
auto constexpr scale(T x, T y, T z) -> std::array<std::array<T, 4>, 4>
{
    return std::array<std::array<T, 4>, 4>{
        {{x, static_cast<T>(0), static_cast<T>(0), static_cast<T>(0)},
         {static_cast<T>(0), y, static_cast<T>(0), static_cast<T>(0)},
         {static_cast<T>(0), static_cast<T>(0), z, static_cast<T>(0)},
         {static_cast<T>(0),
          static_cast<T>(0),
          static_cast<T>(0),
          static_cast<T>(1)}}
    };
}

template <typename T>
auto constexpr translate(std::array<T, 3> translation_vector
) -> std::array<std::array<T, 4>, 4>
{
    return std::array<std::array<T, 4>, 4>{
        {{static_cast<T>(1),
          static_cast<T>(0),
          static_cast<T>(0),
          translation_vector[0]},
         {static_cast<T>(0),
          static_cast<T>(1),
          static_cast<T>(0),
          translation_vector[1]},
         {static_cast<T>(0),
          static_cast<T>(0),
          static_cast<T>(1),
          translation_vector[2]},
         {static_cast<T>(0),
          static_cast<T>(0),
          static_cast<T>(0),
          static_cast<T>(1)}}
    };
}

// converts degrees to radians
template <typename T> static T radians(T degrees)
{
    return degrees * (static_cast<T>(M_PI) / static_cast<T>(180));
}

template <typename T>
auto rotate_z(T degrees) -> std::array<std::array<T, 4>, 4>
{
    std::array<std::array<T, 4>, 4> result{identity<T, 4>()};
    T rad = radians(degrees);
    result[0][0] = std::cos(rad);
    result[0][1] = -std::sin(rad);
    result[1][0] = std::sin(rad);
    result[1][1] = std::cos(rad);
    return result;
}

template <typename T>
auto rotate_y(T degrees) -> std::array<std::array<T, 4>, 4>
{
    std::array<std::array<T, 4>, 4> result{identity<T, 4>()};
    T rad = radians(degrees);
    result[0][0] = std::cos(rad);
    result[0][2] = std::sin(rad);
    result[2][0] = -std::sin(rad);
    result[2][2] = std::cos(rad);
    return result;
}

template <typename T>
auto rotate_x(T degrees) -> std::array<std::array<T, 4>, 4>
{
    std::array<std::array<T, 4>, 4> result{identity<T, 4>()};
    T rad = radians(degrees);
    result[1][1] = std::cos(rad);
    result[1][2] = -std::sin(rad);
    result[2][1] = std::sin(rad);
    result[2][2] = std::cos(rad);
    return result;
}

template <typename T>
auto perspective(T fov_degrees, T aspect_ratio, T near, T far)
    -> std::array<std::array<T, 4>, 4>
{
    // perspective projection
    std::array<std::array<T, 4>, 4> proj{identity<T, 4>()};
    proj[3][3] = 0;
    T fov = radians<T>(fov_degrees);
    T const n = near;
    T const f = far;

    T h = static_cast<T>(2) * n * std::tan(fov / static_cast<T>(2));
    T w = h * aspect_ratio;

    T r = (w / 2.0f);
    T l = -(w / 2.0f);
    T b = -(h / 2.0f);
    T t = (h / 2.0f);

    proj[0][0] = (static_cast<T>(2) * n) / (r - l);
    proj[1][1] = (static_cast<T>(2) * n) / (t - b);
    proj[2][2] = -(f + n) / (f - n);
    proj[2][1] = (t + b) / (t - b);
    proj[2][0] = (r + l) / (r - l);
    proj[2][3] = static_cast<T>(-2) * (f * n) / (f - n);
    proj[3][2] = static_cast<T>(-1);
    return proj;
}

template <typename T>
auto orthogonal(T fov_degrees, T aspect_ratio, T near, T far)
    -> std::array<std::array<T, 4>, 4>
{
    // orthogonal projection
    std::array<std::array<T, 4>, 4> proj{identity<T, 4>()};
    T fov = radians<T>(fov_degrees);
    T const n = near;
    T const f = far;

    T h = static_cast<T>(2) * n * std::tan(fov / static_cast<T>(2));
    T w = h * aspect_ratio;

    T r = (w / 2.0f);
    T l = -(w / 2.0f);
    T b = -(h / 2.0f);
    T t = (h / 2.0f);

    proj[0][0] = (static_cast<T>(2)) / (r - l);
    proj[1][1] = (static_cast<T>(2)) / (t - b);
    proj[2][2] = -2.0f / (f - n);
    proj[3][3] = 1;
    proj[3][0] = -(r + l) / (r - l);
    proj[3][1] = -(t + b) / (t - b);
    proj[3][2] = -(f + n) / (f - n);
    return proj;
}

// TODO: make this not inline idk its scary
inline std::tuple<mat4, vec3, vec3> view(
    aa::Float yaw,
    aa::Float pitch,
    vec3 const& world_up,
    vec3 const& camera_position
)
{
    vec3 rev_target_direction{0, 0, 0};
    rev_target_direction[0] = std::cos(radians(yaw)) * std::cos(radians(pitch));
    rev_target_direction[1] = std::sin(radians(pitch));
    rev_target_direction[2] = std::sin(radians(yaw)) * std::cos(radians(pitch));

    vec3 camera_front = normalize(rev_target_direction);

    vec3 camera_right = normalize(cross(world_up, rev_target_direction));

    vec3 camera_up = normalize(cross(rev_target_direction, camera_right));

    mat4 view_r{identity<aa::Float, 4>()};
    view_r[0][3] = -camera_position[0];
    view_r[1][3] = -camera_position[1];
    view_r[2][3] = -camera_position[2];

    mat4 view_l{identity<aa::Float, 4>()};
    view_l[0][0] = camera_right[0];
    view_l[0][1] = camera_right[1];
    view_l[0][2] = camera_right[2];

    view_l[1][0] = camera_up[0];
    view_l[1][1] = camera_up[1];
    view_l[1][2] = camera_up[2];

    view_l[2][0] = rev_target_direction[0];
    view_l[2][1] = rev_target_direction[1];
    view_l[2][2] = rev_target_direction[2];

    mat4 view_t = view_l * view_r;
    return {view_t, camera_front, camera_up};
}

} // namespace aa

#endif
