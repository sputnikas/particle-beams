#include "main.h"

int main(int argc, char *argv[])
{
	std::cout << "pbeams - program for calculation particle beams in diferent systems" << std::endl;
	#ifdef HAS_TEST
		// test_vec3();
		// test_vec2();
		// test_sphere_method();
		// test_expression();
		// test_number_method();
	#endif // HAS_TEST
	std::string filename("test.dat");
    std::string filename2("test2.dat");
    std::cout << filename << std::endl;

    double I = -3000;
    double dt = 1e-12;
    double By = CL*1;
    double vz = 0.8;
    double Ex = vz*By;
    size_t n = 100;
    double wx = 1e-2;
    double wy = 1e-2;
    size_t N = 30;

    ParticleType *ptype = new ParticleTypePoint(I*dt/N, fabs(I*dt/N/EQ*EM));
    std::cout << "rho = " << I/vz/CL/wx/wy << std::endl;
    Injector *injector = new InjectorRectangleZ(wx, wy, vz, N, ptype, I, I/vz/CL/wx/wy, fabs(I*dt/N/EQ));
    Space3d *space3d = new Space3dHalf(Vec3<double>(0, 0, 1));
    ExternField *extfield = new ExternFieldConst(Vec3<double>(Ex, 0, 0), Vec3<double>(0, By, 0), space3d);

    PPSolver pp;
    pp.extern_fields.push_back(extfield);
    pp.spaces.push_back(space3d);
    pp.particle_types.push_back(ptype);
    pp.injectors.push_back(injector);
    pp.number_method = 3;
    pp.field_interaction = &pp.field_interaction1;
    pp.interaction_type = 0;
    pp.sphere_method = 0;
    pp.nt = 0;
    pp.t = 0.0;
	pp.dt = dt;
	pp.cdt = CL*dt;
	pp.ct = 0.0;

    for (size_t i = 0; i<n; i++) {
        //std::cout << i << std::endl;
        
        //getchar();
        std::cout << i << ": " << pp.calc() << std::endl;
    }
	
    //std::cout << "We are in function" << std::endl;
    pp.save(filename);
    pp.save(filename2, n);
    pp.free();
	return 0;
}
