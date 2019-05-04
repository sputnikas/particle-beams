#include "main.h"

int main() {
    std::string filename("test.dat");
    std::string filename2("test2.dat");
    std::cout << filename << std::endl;

    double I = 10;
    double dt = 1e-12;
    double By = CL*1;
    double vz = 0.8;
    double Ex = vz*By;
    size_t n = 100;
    double wx = 1e-2;
    double wy = 1e-2;
    size_t N = 30;

    ParticleTypePoint *ptype = new ParticleTypePoint(I*dt/N, I*dt/N/EQ*EM);
    InjectorRectangleZ *injector = new InjectorRectangleZ(wx, wy, vz, N, ptype, I, I*dt/vz/CL/wx/wy/dt, fabs(I*dt/N/EQ));
    Space3dHalf *space3d = new Space3dHalf(Vec3<double>(0, 0, 1));
    ExternFieldConst *extfield = new ExternFieldConst(Vec3<double>(Ex, 0, 0), Vec3<double>(0, By, 0), space3d);

    PPSolver pp;
    pp.extern_fields.push_back(extfield);
    pp.spaces.push_back(space3d);
    pp.particle_types.push_back(ptype);
    pp.injectors.push_back(injector);
    pp.number_method = 1;
    pp.field_interaction = &pp.field_interaction0;
    pp.nt = 0;
    pp.t = 0.0;

    for (size_t i = 0; i<n; i++) {
        std::cout << i << std::endl;
        pp.calc();
    }

    pp.save(filename);

    return 0;
}