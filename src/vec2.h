#pragma once

#include <cmath>
#include <iostream>
#include <complex>

#if __cplusplus > 199711L
#define __VEC2T(a, b) {(a), (b)}
#else
#define __VEC2T(a, b) Vec2<T>((a), (b))
#endif // __cplusplus

#ifndef INLINE
#define INLINE inline
#endif // INLINE

template <typename T>
class Vec2 {
	public:
		T x;
		T y;

		#ifdef HAS_TEST_CONSTRUCTOR
		static int index_constructor;
		#endif // HAS_TEST_CONSTRUCTOR

		Vec2();
		Vec2(const T &x, const T &y);
		Vec2(const Vec2<T> &a);
		~Vec2();

		// Операторы сложения и вычитания
		INLINE Vec2<T> operator + (const Vec2<T> &a) const;
		INLINE Vec2<T> operator - (const Vec2<T> &a) const;
		INLINE Vec2<T> operator - ();
		INLINE void operator += (const Vec2<T> &a);
		INLINE void operator -= (const Vec2<T> &a);
		// Операторы умножения
		INLINE Vec2<T> operator * (const T &f) const;
		INLINE void operator *= (const T &f);
		INLINE Vec2<T> operator / (const T &f) const;
		INLINE void operator /= (const T &f);
		// Скалярное произведение
		INLINE T dot(const Vec2<T> &a) const;
		// Векторное произведение - краеугольный камень любой шаблонной реализации
		// если размерность пространства = N, то векторное произведение определяется
		// как вектор, нормальный к подпространству из N-1 векторов и по длине равный
		// N-мерному объёму параллелипепипеда построенного на N-1 векторах + направление
		// определяется способом обхода! Получается, что количество аргументов функции cross
		// равно N - 1 для внешней функции и N-2 для функции класса, то есть зависит от параметра шаблона!
		INLINE Vec2<T> cross() const;
		// Нормы
		INLINE T norm() const;
		INLINE T norm2() const;
		// Операторы сравнения
		INLINE bool operator < (const T &f) const;
		INLINE bool operator > (const T &f) const;
		INLINE bool operator == (const Vec2<T> &a) const;
};

#ifdef HAS_TEST
	void test_vec2();
#endif // HAS_TEST

////////////////////////////////////////////////////////////////////////////////////////////////////
// implementation for template functions
////////////////////////////////////////////////////////////////////////////////////////////////////

template <typename T>
Vec2<T>::Vec2() : x(0), y(0)
{
	#ifdef HAS_TEST_CONSTRUCTOR
	index_constructor++;
	std::cout << "constructor Vec2() : " << index_constructor << std::endl;
	#endif // HAS_TEST_CONSTRUCTOR
};

template <typename T>
Vec2<T>::Vec2(const T &x, const T &y)
	: x(x), y(y)
{
	#ifdef HAS_TEST_CONSTRUCTOR
	index_constructor++;
	std::cout << "constructor Vec2(x,y) : " << index_constructor << std::endl;
	#endif // HAS_TEST_CONSTRUCTOR
};

template <typename T>
Vec2<T>::Vec2(const Vec2<T> &a) : x(a.x), y(a.y) 
{
	#ifdef HAS_TEST_CONSTRUCTOR
	index_constructor++;
	std::cout << "copy constructor Vec2 : " << index_constructor << std::endl;
	#endif // HAS_TEST_CONSTRUCTOR
};

template <typename T>
Vec2<T>::~Vec2()
{
	#ifdef HAS_TEST_CONSTRUCTOR
	index_constructor--;
	std::cout << "destructor Vec2(const Vec2<T> &) : " << index_constructor << std::endl;
	#endif // HAS_TEST_CONSTRUCTOR
};

template <typename T>
std::ostream &operator << (std::ostream &os, const Vec2<T> &a)
{
	os << "{" << a.x << ", " << a.y << "}";
	return os;
};

// Операторы сложения и вычитания
template <typename T>
INLINE Vec2<T> Vec2<T>::operator + (const Vec2<T> &a) const
{
	return __VEC2T(x + a.x, y + a.y);
};

template <typename T>
INLINE Vec2<T> Vec2<T>::operator - (const Vec2<T> &a) const
{
	return __VEC2T(x - a.x, y - a.y);
};

template <typename T>
INLINE Vec2<T> Vec2<T>::operator - ()
{
	return __VEC2T(- x, - y);
};

template <typename T>
INLINE void Vec2<T>::operator += (const Vec2<T> &a)
{
	x += a.x;
	y += a.y;
};

template <typename T>
INLINE void Vec2<T>::operator -= (const Vec2<T> &a)
{
	x -= a.x;
	y -= a.y;
};

// Операторы умножения
template <typename T>
INLINE Vec2<T> operator * (const T &f, const Vec2<T> &a)
{
	return __VEC2T(f * a.x, f * a.y);
};

template <typename T>
INLINE Vec2<T> Vec2<T>::operator * (const T &f) const
{
	return __VEC2T(f * x, f * y);
};

template <typename T>
INLINE void Vec2<T>::operator *= (const T &f)
{
	x *= f;
	y *= f;
};

template <typename T>
INLINE Vec2<T> Vec2<T>::operator / (const T &f) const
{
	return __VEC2T(x / f, y / f);
};

template <typename T>
INLINE void Vec2<T>::operator /= (const T &f)
{
	x /= f;
	y /= f;
};

// Скалярное произведение
template <typename T>
INLINE T Vec2<T>::dot(const Vec2<T> &a) const
{
	return a.x * x + a.y * y;
};

// Векторное произведение
template <typename T>
INLINE Vec2<T> Vec2<T>::cross() const
{
	return __VEC2T(y, -x);
};

// Нормы
template <typename T>
INLINE T Vec2<T>::norm2() const
{
	return x*x + y*y;
};

template <typename T>
INLINE T Vec2<T>::norm() const
{
	return std::sqrt(x*x + y*y);
};

// Операторы сравнения
template <typename T>
INLINE bool Vec2<T>::operator < (const T &f) const
{
	return norm2() < f*f;
};

template <typename T>
INLINE bool Vec2<T>::operator > (const T &f) const
{
	return norm2() > f*f;
};

template <typename T>
INLINE bool Vec2<T>::operator == (const Vec2<T> &a) const
{
	return ((x - a.x) == 0) && ((y - a.y) == 0);
};
