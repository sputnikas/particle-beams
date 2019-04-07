#pragma once

#include <cmath>
#include <iostream>
#include <complex>

#if __cplusplus > 199711L
#define __VEC3T(a, b, c) {(a), (b), (c)}
#else
#define __VEC3T(a, b, c) Vec3<T>((a), (b), (c))
#endif // __cplusplus

#ifndef INLINE
#define INLINE inline
#endif // INLINE

template <typename T>
class Vec3 {
	public:
		T x;
		T y;
		T z;

		#ifdef HAS_TEST_CONSTRUCTOR
		static int index_constructor;
		#endif // HAS_TEST_CONSTRUCTOR

		Vec3();
		Vec3(const T &x, const T &y, const T &z);
		Vec3(const Vec3<T> &a);
		~Vec3();

		// Операторы сложения и вычитания
		INLINE Vec3<T> operator + (const Vec3<T> &a) const;
		INLINE Vec3<T> operator - (const Vec3<T> &a) const;
		INLINE Vec3<T> operator - ();
		INLINE void operator += (const Vec3<T> &a);
		INLINE void operator -= (const Vec3<T> &a);
		// Операторы умножения
		INLINE Vec3<T> operator * (const T &f) const;
		INLINE void operator *= (const T &f);
		INLINE Vec3<T> operator / (const T &f) const;
		INLINE void operator /= (const T &f);
		// Скалярное произведение
		INLINE T dot(const Vec3<T> &a) const;
		// Векторное произведение - краеугольный камень любой шаблонной реализации размерности
		// если размерность пространства = N, то векторное произведение определяется
		// как вектор, нормальный к подпространству из N-1 векторов и по длине равный
		// N-мерному объёму параллелипепипеда построенного на N-1 векторах + направление
		// определяется способом обхода! Получается, что количество аргументов функции cross
		// равно N - 1 для внешней функции и N - 2 для функции класса, то есть зависит от параметра шаблона!
		INLINE Vec3<T> cross(const Vec3<T> &a) const;
		// Нормы
		INLINE T norm() const;
		INLINE T norm2() const;
		// Операторы сравнения
		INLINE bool operator < (const T &f) const;
		INLINE bool operator > (const T &f) const;
		INLINE bool operator == (const Vec3<T> &a) const;
};

#ifdef HAS_TEST
	void test_vec3();
#endif // HAS_TEST

////////////////////////////////////////////////////////////////////////////////////////////////////
// implementation for template functions
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
Vec3<T>::Vec3() : x(0), y(0), z(0)
{
	#ifdef HAS_TEST_CONSTRUCTOR
	index_constructor++;
	std::cout << "constructor Vec3() : " << index_constructor << std::endl;
	#endif // HAS_TEST_CONSTRUCTOR
};

template <typename T>
Vec3<T>::Vec3(const T &x, const T &y, const T &z)
	: x(x), y(y), z(z)
{
	#ifdef HAS_TEST_CONSTRUCTOR
	index_constructor++;
	std::cout << "constructor Vec3(x,y) : " << index_constructor << std::endl;
	#endif // HAS_TEST_CONSTRUCTOR
};

template <typename T>
Vec3<T>::Vec3(const Vec3<T> &a) : x(a.x), y(a.y), z(a.z)
{
	#ifdef HAS_TEST_CONSTRUCTOR
	index_constructor++;
	std::cout << "copy constructor Vec3 : " << index_constructor << std::endl;
	#endif // HAS_TEST_CONSTRUCTOR
};

template <typename T>
Vec3<T>::~Vec3()
{
	#ifdef HAS_TEST_CONSTRUCTOR
	index_constructor--;
	std::cout << "destructor Vec3(const Vec3<T> &) : " << index_constructor << std::endl;
	#endif // HAS_TEST_CONSTRUCTOR
};

template <typename T>
std::ostream &operator << (std::ostream &os, const Vec3<T> &a)
{
	os << "{" << a.x << ", " << a.y << ", " << a.z << "}";
	return os;
};

// Операторы сложения и вычитания
template <typename T>
INLINE Vec3<T> Vec3<T>::operator + (const Vec3<T> &a) const
{
	return __VEC3T(x + a.x, y + a.y, z + a.z);
};

template <typename T>
INLINE Vec3<T> Vec3<T>::operator - (const Vec3<T> &a) const
{
	return __VEC3T(x - a.x, y - a.y, z - a.z);
};

template <typename T>
INLINE Vec3<T> Vec3<T>::operator - ()
{
	return __VEC3T(- x, - y, - z);
};

template <typename T>
INLINE void Vec3<T>::operator += (const Vec3<T> &a)
{
	x += a.x;
	y += a.y;
	z += a.z;
};

template <typename T>
INLINE void Vec3<T>::operator -= (const Vec3<T> &a)
{
	x -= a.x;
	y -= a.y;
	z -= a.z;
};

// Операторы умножения
template <typename T>
INLINE Vec3<T> operator * (const T &f, const Vec3<T> &a)
{
	return __VEC3T(f * a.x, f * a.y, f * a.z);
};

template <typename T>
INLINE Vec3<T> Vec3<T>::operator * (const T &f) const
{
	return __VEC3T(f * x, f * y, f * z);
};

template <typename T>
INLINE void Vec3<T>::operator *= (const T &f)
{
	x *= f;
	y *= f;
	z *= f;
};

template <typename T>
INLINE Vec3<T> Vec3<T>::operator / (const T &f) const
{
	return __VEC3T(x / f, y / f, z / f);
};

template <typename T>
INLINE void Vec3<T>::operator /= (const T &f)
{
	x /= f;
	y /= f;
	z /= f;
};

// Скалярное произведение
template <typename T>
INLINE T Vec3<T>::dot(const Vec3<T> &a) const
{
	return a.x * x + a.y * y + a.z * z;
};

// Векторное произведение
template <typename T>
INLINE Vec3<T> Vec3<T>::cross(const Vec3<T> &a) const
{
	return __VEC3T(y*a.z - z*a.y, z*a.x - x*a.z, x*a.y - y*a.x);
};

// Нормы
template <typename T>
INLINE T Vec3<T>::norm2() const
{
	return x*x + y*y + z*z;
};

template <typename T>
INLINE T Vec3<T>::norm() const
{
	return std::sqrt(x*x + y*y + z*z);
};

// Операторы сравнения
template <typename T>
INLINE bool Vec3<T>::operator < (const T &f) const
{
	return norm2() < f*f;
};

template <typename T>
INLINE bool Vec3<T>::operator > (const T &f) const
{
	return norm2() > f*f;
};

template <typename T>
INLINE bool Vec3<T>::operator == (const Vec3<T> &a) const
{
	return ((x - a.x) == 0) && ((y - a.y) == 0) && ((z - a.z) == 0);
};
