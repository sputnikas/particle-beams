#pragma once

#include <cmath>
#include <iostream>
#include <complex>

#if __cplusplus > 199711L
#define __VEC3T(a, b, c) {(a), (b), (c)}
#else
#define __VEC3T(a, b, c) Vec3<T>((a), (b), (c))
#endif // __cplusplus

template <typename T>
class Vec3 {
	public:
		T x;
		T y;
		T z;

		Vec3();
		Vec3(const T &, const T &, const T &);
		Vec3(const Vec3<T> &);
		~Vec3();
};

template <typename T> std::ostream &operator<< (std::ostream &, const Vec3<T> &);

// Операторы сложения и вычитания
template <typename T> inline Vec3<T> operator + (const Vec3<T> &, const Vec3<T> &);
template <typename T> inline Vec3<T> operator - (const Vec3<T> &, const Vec3<T> &);
template <typename T> inline Vec3<T> operator - (const Vec3<T> &);
template <typename T> inline void operator += (Vec3<T> &, const Vec3<T> &);
template <typename T> inline void operator -= (Vec3<T> &, const Vec3<T> &);
// Операторы умножения
template <typename T> inline Vec3<T> operator * (const T &, const Vec3<T> &);
template <typename T> inline Vec3<T> operator * (const Vec3<T> &, const T &);
template <typename T> inline void operator *= (Vec3<T> &, const T &);
template <typename T> inline Vec3<T> operator * (const Vec3<T> &, const Vec3<T> &);
template <typename T> inline void operator *= (Vec3<T> &, const Vec3<T> &);
template <typename T> inline Vec3<T> operator / (const Vec3<T> &, const T &);
template <typename T> inline void operator /= (Vec3<T> &, const T &);
// Скалярное произведение
template <typename T> inline T dot(const Vec3<T> &, const Vec3<T> &);
// Нормы
template <typename T> inline T norm(const Vec3<T> &);
template <typename T> inline T norm2(const Vec3<T> &);
// Операторы сравнения
template <typename T> inline bool operator < (const Vec3<T> &, const T &);
template <typename T> inline bool operator > (const Vec3<T> &, const T &);
template <typename T> inline bool operator == (const Vec3<T> &, const Vec3<T> &);

// Тестовая функция
#ifdef HAS_TEST
	void test_vec3();
#endif // HAS_TEST

////////////////////////////////////////////////////////////////////////////////////////////////////
// implementation for template functions
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
Vec3<T>::Vec3() : x(0), y(0), z(0)
{
}

template <typename T>
Vec3<T>::Vec3(const T &x, const T &y, const T &z)
	: x(x), y(y), z(z)
{
}

template <typename T>
Vec3<T>::Vec3(const Vec3<T> &a) : x(a.x), y(a.y), z(a.z)
{
}

template <typename T>
Vec3<T>::~Vec3()
{
}

template <typename T>
std::ostream &operator<< (std::ostream &os, const Vec3<T> &a)
{
	os << "{" << a.x << ", " << a.y << ", " << a.z << "}";
	return os;
};

// Операторы сложения и вычитания
template <typename T>
inline Vec3<T> operator + (const Vec3<T> &a, const Vec3<T> &b)
{
	return __VEC3T(a.x + b.x, a.y + b.y, a.z + b.z);
};

template <typename T>
inline Vec3<T> operator - (const Vec3<T> &a, const Vec3<T> &b)
{
	return __VEC3T(a.x - b.x, a.y - b.y, a.z - b.z);
};

template <typename T>
inline Vec3<T> operator - (const Vec3<T> &a)
{
	return __VEC3T(- a.x, - a.y, - a.z);
};

template <typename T>
inline void operator += (Vec3<T> &a, const Vec3<T> &b)
{
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
};

template <typename T>
inline void operator -= (Vec3<T> &a, const Vec3<T> &b)
{
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
};

// Операторы умножения
template <typename T>
inline Vec3<T> operator * (const T &f, const Vec3<T> &a)
{
	return __VEC3T(f * a.x, f * a.y, f * a.z);
};

template <typename T>
inline Vec3<T> operator * (const Vec3<T> &a, const T &f)
{
	return __VEC3T(f * a.x, f * a.y, f * a.z);
};

template <typename T>
inline void operator *= (Vec3<T> &a, const T &f)
{
	a.x *= f;
	a.y *= f;
	a.z *= f;
};

template <typename T>
inline Vec3<T> operator / (const Vec3<T> &a, const T &f)
{
	return __VEC3T(a.x / f, a.y / f, a.z / f);
};

template <typename T>
inline void operator /= (Vec3<T> &a, const T &f)
{
	a.x /= f;
	a.y /= f;
	a.z /= f;
};

template <typename T>
inline Vec3<T> operator * (const Vec3<T> &a, const Vec3<T> &b)
{
	return __VEC3T(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
};

template <typename T>
inline void operator *= (Vec3<T> &a, const Vec3<T> &b)
{
	T tx = a.x, ty = a.y;
	a.x = a.y * b.z - a.z * b.y;
	a.y = a.z * b.x - tx * b.z;
	a.z = tx * b.y  - ty * b.z;
};

// Скалярное произведение
template <typename T>
inline T dot(const Vec3<T> &a, const Vec3<T> &b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
};

template <>
inline std::complex<float> dot<std::complex<float>> (const Vec3<std::complex<float>> &a,
		const Vec3<std::complex<float>> &b)
{
	return std::conj(a.x) * b.x + std::conj(a.y) * b.y + std::conj(a.z) * b.z;
};

template <>
inline std::complex<double> dot<std::complex<double>> (const Vec3<std::complex<double>> &a,
		const Vec3<std::complex<double>> &b)
{
	return std::conj(a.x) * b.x + std::conj(a.y) * b.y + std::conj(a.z) * b.z;
};

// Нормы
template <typename T>
inline T norm2(const Vec3<T> &a)
{
	return dot(a, a);
};

template <typename T>
inline T norm(const Vec3<T> &a)
{
	return sqrt(norm2(a));
};

// Операторы сравнения
template <typename T>
inline bool operator < (const Vec3<T> &a, const T &f)
{
	return norm2(a) < f*f;
};

template <typename T>
inline bool operator > (const Vec3<T> &a, const T &f)
{
	return norm2(a) > f*f;
};

template <typename T>
inline bool operator == (const Vec3<T> &a, const Vec3<T> &b)
{
	return norm2(a - b) == 0;
};

