#include "ISOP2P1.h"

void ISOP2P1::config(std::string _config_file)
{
	std::string trash;
	std::ifstream input(_config_file.c_str());
	input >> trash >> trash >> mesh_file;
	input >> trash >> trash >> l_tol;
	input >> trash >> trash >> l_Euler_tol;
	input >> trash >> trash >> n_tol;
	input >> trash >> trash >> record;
	input >> trash >> trash >> check;
	input >> trash >> trash >> error_check;
	input >> trash >> trash >> viscosity;
	input >> trash >> trash >> t0;
	input >> trash >> trash >> t1;
	input >> trash >> trash >> t;
	input >> trash >> trash >> dt;
	input >> trash >> trash >> n_method;
	input >> trash >> trash >> time_step_control;
	input >> trash >> trash >> Stokes_init;
	input >> trash >> trash >> NS_init;
	input >> trash >> trash >> scheme;
	input >> trash >> trash >> isMoving;
	input >> trash >> trash >> scale;
	input >> trash >> trash >> max_step;
	input >> trash >> trash >> alpha;
	input >> trash >> trash >> beta;
	input >> trash >> trash >> output_vorticity;
	input >> trash >> trash >> output_divergence;
	input >> trash >> trash >> body_force;
	input >> trash >> trash >> angle;
	input.close();
};

void ISOP2P1::config_debug()
{
	std::cout << "mesh_file = " << mesh_file << std::endl;
	std::cout << "l_tol"<< l_tol << std::endl;
	std::cout << "l_Euler_tol =" << l_Euler_tol << std::endl;
	std::cout << "n_tol = " << n_tol << std::endl;
	std::cout << "record = " << record << std::endl;
	std::cout << "check =" << check << std::endl;
	std::cout << "error_check =" << error_check << std::endl;
	std::cout << "viscosity =" << viscosity << std::endl;
	std::cout << "t0 =" << t0 << std::endl;
	std::cout << "t1 =" << t1 << std::endl;
	std::cout << "t =" << t << std::endl;
	std::cout << "dt =" << dt << std::endl;
	std::cout << "n_method =" << n_method << std::endl;
	std::cout << "time_step_control =" << time_step_control << std::endl;
	std::cout << "Stokes_init =" << Stokes_init << std::endl;
	std::cout << "NS_init =" << NS_init << std::endl;
	std::cout << "scheme =" << scheme << std::endl;
	std::cout << "isMoving =" << isMoving << std::endl;
	std::cout << "scale =" << scale << std::endl;
	std::cout << "max_step =" << max_step << std::endl;
	std::cout << "alpha =" << alpha << std::endl;
	std::cout << "beta =" << beta << std::endl;
	std::cout << "output_vorticity =" << output_vorticity << std::endl;
	std::cout << "output_divergence =" << output_divergence << std::endl;
	std::cout << "body_force =" << body_force << std::endl;
	std::cout << "angle = " << angle << std::endl;

}
