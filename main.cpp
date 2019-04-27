#include "main.h"

int main()
{
	std::cout << "pbeams - program for calculation particle beams in diferent systems" << std::endl;
	#ifdef HAS_TEST
		// test_vec3();
		// test_vec2();
		test_sphere_method();
		//test_expression();
		//test_number_method();
	#endif // HAS_TEST

	return 0;
}
