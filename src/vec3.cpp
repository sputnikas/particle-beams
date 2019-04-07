#include "vec3.h"

#ifdef HAS_TEST_CONSTRUCTOR
template<typename T>
int Vec3<T>::index_constructor = 0;
#endif // HAS_TEST_CONSTRUCTOR

template <>
INLINE std::complex<float> Vec3<std::complex<float>>::norm2() const
{
	return std::conj(x)*x + std::conj(y)*y + std::conj(z)*z;
};

template <>
INLINE std::complex<double> Vec3<std::complex<double>>::norm2() const
{
	return std::conj(x)*x + std::conj(y)*y + std::conj(z)*z;
};

template <>
INLINE std::complex<float> Vec3<std::complex<float>>::norm() const
{
	return std::sqrt(std::conj(x)*x + std::conj(y)*y) + std::conj(z)*z;
};

template <>
INLINE std::complex<double> Vec3<std::complex<double>>::norm() const
{
	return std::sqrt(std::conj(x)*x + std::conj(y)*y) + std::conj(z)*z;
};

template <>
INLINE std::complex<float> Vec3<std::complex<float>>::dot(const Vec3<std::complex<float>> &a) const
{
	return std::conj(a.x) * x + std::conj(a.y) * y + std::conj(z)*z;
};

template <>
INLINE std::complex<double> Vec3<std::complex<double>>::dot(const Vec3<std::complex<double>> &a) const
{
	return std::conj(a.x) * x + std::conj(a.y) * y + std::conj(z)*z;
};

#ifdef HAS_TEST

void test_vec3()
{
	Vec3<double> a = Vec3<double> ( 1, 2, 3 );
	Vec3<double> b = Vec3<double> ( 3, 2, 1 );
	Vec3<double> c;
	std::cout << "a = " << a << std::endl;
	std::cout << "b = " << b << std::endl;
	std::cout << "c = " << c << std::endl;
	std::cout << "a + b = " << a + b << std::endl;
	c += a;
	std::cout << "c += a \nc = " << c << std::endl;
	c -= b;
	std::cout << "c -= b" << std::endl;
	std::cout << "c = a - b = " << a - b << " = " << c << std::endl;
	std::cout << "- a = " << - a << std::endl;
	c = a;
	c *= 2.;
	std::cout << "c = a" << std::endl;
	std::cout << "c *= 2" << std::endl;
	std::cout << "c = 2*a = " << a * 2. << " = " << 2.*a << " = " << c << std::endl;
	std::cout << "6*a = " << 3.*a * 2. << std::endl;
	c = a;
	c /= 2.;
	std::cout << "c = a" << std::endl;
	std::cout << "c /= 2" << std::endl;
	std::cout << "a/2 = " << a / 2. << " = " << c << std::endl;
	std::cout << "a.cross(b) = " << a.cross(b) << std::endl;
	std::cout << "a.dot(b) = " << a.dot(b) << std::endl;
	std::cout << "a.norm() = " << a.norm() << std::endl;
	std::cout << "a.norm2() = " << a.norm2() << std::endl;
	std::cout << "a < 3 : " << ( a < 3. ) << std::endl;
	std::cout << "a > 3 : " << ( a > 3. ) << std::endl;
	std::cout << "a == b : " << ( a == b ) << std::endl;
	std::cout << "a == a : " << ( a == a ) << std::endl;
};

#endif // HAS_TEST
