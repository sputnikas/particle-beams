#include "vec2.h"

#ifdef HAS_TEST

void test_vec2()
{
	Vec2<double> a = Vec2<double> ( 1, 2 );
	Vec2<double> b = Vec2<double> ( 3, 2 );
	Vec2<double> c;
	std::cout << "a = " << a << std::endl;
	std::cout << "b = " << b << std::endl;
	std::cout << "c = " << c << std::endl;
	std::cout << "a + b = " << a + b << std::endl;
	c += a;
	std::cout << "c = a" << c << std::endl;
	c -= b;
	std::cout << "a - b = " << a - b << " " << c << std::endl;
	std::cout << "- a = " << - a << std::endl;
	c = a;
	c *= 2.;
	std::cout << "2*a = " << a * 2. << " " << 2.*a << " " << c << std::endl;
	c = a;
	c /= 2.;
	std::cout << "a/2 = " << a / 2. << " " << c << std::endl;
	std::cout << "a*b = " << a * b << std::endl;
	std::cout << "dot(a, b) = " << dot ( a, b ) << std::endl;
	std::cout << "norm(a) = " << norm ( a ) << std::endl;
	std::cout << "norm2(a) = " << norm2 ( a ) << std::endl;
	std::cout << "a < 3 " << ( a < 3. ) << std::endl;
	std::cout << "a > 3 " << ( a > 3. ) << std::endl;
	std::cout << "a == b " << ( a == b ) << std::endl;
	std::cout << "a == a " << ( a == a ) << std::endl;
};

#endif // HAS_TEST
