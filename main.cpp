#include "main.h"

// PPSolver load_ini(const std::string &filename) {
//     const int NBUF = 512;
//     const int NINI = 256;
//     char line[NBUF];
//     char uline[NBUF];
//     char key[NINI];
//     char val[NINI];

//     PPSolver ps;
//     FILE *fd = fopen(filename.c_str(), "r");
//     if (fd != NULL) {
//     if (fgets(line, sizeof(line), fd) != NULL) {
//         while (feof(fd) == 0) {
//             if (strcmp(line, "[Main]") == 0) {
//                 while (fgets(uline, sizeof(uline), fd) != NULL) {
//                     if (uline[0] == '[') {
//                         strcpy(line, uline);
//                         break;
//                     }
//                     if (strspn(uline, "; ") == strspn(uline, " ")) {
//                         if (sscanf(uline, "%s = %[^;\n]", key, val) == 2) {
//                             if (strcmp(key, "dt")) sscanf(val, "%le", &ps.dt);
//                             if (strcmp(key, "ncounts")) sscanf(val, "%llu", &ps.ncounts);
//                             if (strcmp(key, "number_method")) sscanf(val, "%llu", &ps.number_method);
//                             if (strcmp(key, "interaction_relativistic")) sscanf(val, "%llu", &ps.interaction_type);
//                             size_t interaction_type;
//                             if (strcmp(key, "interaction_type")) sscanf(val, "%llu", &interaction_type);
//                             switch (interaction_type) {
//                             case 1: 
//                                 ps.field_interaction = PPSolver::field_interaction1;
//                                 break;
//                             default:
//                                 ps.field_interaction = PPSolver::field_interaction0;
//                             }
//                             if (strcmp(key, "sphere_method")) sscanf(val, "%llu", &ps.sphere_method);
//                         }
//                     }
//                 }
//             }
//             if (strcmp(line, "[Injector]") == 0) {
//                 while (fgets(uline, sizeof(uline), fd) != NULL) {
//                     if (uline[0] == '[') {
//                         strcpy(line, uline);
//                         break;
//                     }
//                     if (strspn(uline, "; ") == strspn(uline, " ")) {
//                         if (sscanf(uline, "%s = %[^;\n]", key, val) == 2) {
//                             size_t injector_type;
//                             if (strcmp(key, "type")) sscanf(val, "%llu", &injector_type);
//                             if (injector_type == 0) {
//                                 InjectorRectangleZ *injector = new InjectorRectangleZ();
//                                 if (strcmp(key, "centerx")) sscanf(val, "%le", &injector->center.x);
//                                 if (strcmp(key, "centery")) sscanf(val, "%le", &injector->center.y);
//                                 if (strcmp(key, "centerz")) sscanf(val, "%le", &injector->center.z);
//                                 if (strcmp(key, "wx")) sscanf(val, "%le", &injector->wx);
//                                 if (strcmp(key, "wy")) sscanf(val, "%le", &injector->wy);
//                                 if (strcmp(key, "vz")) sscanf(val, "%le", &injector->vz);
//                                 if (strcmp(key, "I")) sscanf(val, "%le", &injector->I);
//                                 if (strcmp(key, "nparticles")) sscanf(val, "%le", &injector->N);
//                                 size_t particle_type;
//                                 size_t particle_type_form;
//                                 if (strcmp(key, "particle_type")) sscanf(val, "%llu", &particle_type);
//                                 if (strcmp(key, "particle_type_form")) sscanf(val, "%llu", &particle_type);
//                                 if (particle_type == 0) {
                                    
//                                 }
//                             }
//                             if (strcmp(key, "ncounts")) sscanf(val, "%llu", &ps.ncounts);
//                             if (strcmp(key, "number_method")) sscanf(val, "%llu", &ps.number_method);
//                             if (strcmp(key, "interaction_relativistic")) sscanf(val, "%llu", &ps.interaction_type);
//                             if (strcmp(key, "interaction_type")) sscanf(val, "%llu", &interaction_type);
//                             if (strcmp(key, "sphere_method")) sscanf(val, "%llu", &ps.sphere_method);
//                         }
//                     }
//                 }
//             }
//         }
//     }
//     fclose(fd);
// }

int main(int argc, char *argv[])
{
    unsigned NUM_THREADS = 2;
    omp_set_num_threads(NUM_THREADS);
    srand((unsigned) time(NULL));
	std::cout << "pbeams - program for calculation particle beams in diferent systems" << std::endl;
	#ifdef HAS_TEST
		// test_vec3();
		// test_vec2();
		// test_sphere_method();
		// test_expression();
		// test_number_method();
        // test_ppsolver1();
        // test_ppsolver2("test2.dat", 600, 10000, 1e-12, 0.1, 0.001, 0);
	#endif // HAS_TEST

    PPSolver ps;
    ps.load("test2.dat");
    // ps.save("test2.025.dat",  25-1);
    // ps.save("test2.050.dat",  50-1);
    //ps.save("test2.100.dat", 100-1);
    ps.save("test2.200.dat", 200-1);
    ps.save("test2.300.dat", 300-1);
    // ps.save("test2.400.dat", 400-1);
    // ps.save("test2.500.dat", 500-1);
    // ps.save("test2.600.dat", 600-1);
    ps.free();

	return 0;
}
