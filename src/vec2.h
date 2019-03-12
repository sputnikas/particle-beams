#pragma once

#include <cmath>
#include <iostream>
#include <complex>

#if __cplusplus > 199711L
#define __VEC2T(a, b) {(a), (b)}
#else
#define __VEC2T(a, b) Vec2<T>((a), (b))
#endif // __cplusplus

template <typename T>
class Vec2 {
	public:
		T x;
		T y;

		Vec2();
		Vec2(const T &, const T &);
		Vec2(const Vec2<T> &);
		~Vec2();
};

template <typename T> std::ostream &operator<< (std::ostream &, const Vec2<T> &);

// Операторы сложения и вычитания
template <typename T> inline Vec2<T> operator + (const Vec2<T> &, const Vec2<T> &);
template <typename T> inline Vec2<T> operator - (const Vec2<T> &, const Vec2<T> &);
template <typename T> inline Vec2<T> operator - (const Vec2<T> &);
template <typename T> inline void operator += (Vec2<T> &, const Vec2<T> &);
template <typename T> inline void operator -= (Vec2<T> &, const Vec2<T> &);
// Операторы умножения
template <typename T> inline Vec2<T> operator * (const T &, const Vec2<T> &);
template <typename T> inline Vec2<T> operator * (const Vec2<T> &, const T &);
template <typename T> inline void operator *= (Vec2<T> &, const T &);
template <typename T> inline T operator * (const Vec2<T> &, const Vec2<T> &);
template <typename T> inline Vec2<T> operator / (const Vec2<T> &, const T &);
template <typename T> inline void operator /= (Vec2<T> &, const T &);
// Скалярное произведение
template <typename T> inline T dot(const Vec2<T> &, const Vec2<T> &);
// Нормы
template <typename T> inline T norm(const Vec2<T> &);
template <typename T> inline T norm2(const Vec2<T> &);
// Операторы сравнения
template <typename T> inline bool operator < (const Vec2<T> &, const T &);
template <typename T> inline bool operator > (const Vec2<T> &, const T &);
template <typename T> inline bool operator == (const Vec2<T> &, const Vec2<T> &);

// Тестовая функция
#ifdef HAS_TEST
	void test_vec2();
#endif // HAS_TEST

////////////////////////////////////////////////////////////////////////////////////////////////////
// implementation for template functions
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
Vec2<T>::Vec2() : x(0), y(0)
{
}

template <typename T>
Vec2<T>::Vec2(const T &x, const T &y)
	: x(x), y(y)
{
}

template <typename T>
Vec2<T>::Vec2(const Vec2<T> &a) : x(a.x), y(a.y)
{
}

template <typename T>
Vec2<T>::~Vec2()
{
}

template <typename T>
std::ostream &operator<< (std::ostream &os, const Vec2<T> &a)
{
	os << "{" << a.x << ", " << a.y << "}";
	return os;
};

// Операторы сложения и вычитания
template <typename T>
inline Vec2<T> operator + (const Vec2<T> &a, const Vec2<T> &b)
{
	return __VEC2T(a.x + b.x, a.y + b.y);
};

template <typename T>
inline Vec2<T> operator - (const Vec2<T> &a, const Vec2<T> &b)
{
	return __VEC2T(a.x - b.x, a.y - b.y);
};

template <typename T>
inline Vec2<T> operator - (const Vec2<T> &a)
{
	return __VEC2T(- a.x, - a.y);
};

template <typename T>
inline void operator += (Vec2<T> &a, const Vec2<T> &b)
{
	a.x += b.x;
	a.y += b.y;
};

template <typename T>
inline void operator -= (Vec2<T> &a, const Vec2<T> &b)
{
	a.x -= b.x;
	a.y -= b.y;
};

// Операторы умножения
template <typename T>
inline Vec2<T> operator * (const T &f, const Vec2<T> &a)
{
	return __VEC2T(f * a.x, f * a.y);
};

template <typename T>
inline Vec2<T> operator * (const Vec2<T> &a, const T &f)
{
	return __VEC2T(f * a.x, f * a.y);
};

template <typename T>
inline void operator *= (Vec2<T> &a, const T &f)
{
	a.x *= f;
	a.y *= f;
};

template <typename T>
inline Vec2<T> operator / (const Vec2<T> &a, const T &f)
{
	return __VEC2T(a.x / f, a.y / f);
};

template <typename T>
inline void operator /= (Vec2<T> &a, const T &f)
{
	a.x /= f;
	a.y /= f;
};

template <typename T>
inline T operator * (const Vec2<T> &a, const Vec2<T> &b)
{
	return a.x * b.y - a.y * b.x;
};

// Скалярное произведение
template <typename T>
inline T dot(const Vec2<T> &a, const Vec2<T> &b)
{
	return a.x * b.x + a.y * b.y;
};

template <>
inline std::complex<float> dot<std::complex<float>> (const Vec2<std::complex<float>> &a,
		const Vec2<std::complex<float>> &b)
{
	return std::conj(a.x) * b.x + std::conj(a.y) * b.y;
};

template <>
inline std::complex<double> dot<std::complex<double>> (const Vec2<std::complex<double>> &a,
		const Vec2<std::complex<double>> &b)
{
	return std::conj(a.x) * b.x + std::conj(a.y) * b.y;
};

// Нормы
template <typename T>
inline T norm2(const Vec2<T> &a)
{
	return dot(a, a);
};

template <typename T>
inline T norm(const Vec2<T> &a)
{
	return sqrt(norm2(a));
};

// Операторы сравнения
template <typename T>
inline bool operator < (const Vec2<T> &a, const T &f)
{
	return norm2(a) < f*f;
};

template <typename T>
inline bool operator > (const Vec2<T> &a, const T &f)
{
	return norm2(a) > f*f;
};

template <typename T>
inline bool operator == (const Vec2<T> &a, const Vec2<T> &b)
{
	return norm2(a - b) == 0;
};
