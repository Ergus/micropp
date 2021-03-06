/*
 *  This is a test example for MicroPP: a finite element library
 *  to solve microstructural problems for composite materials.
 *
 *  Copyright (C) - 2018 - Jimmy Aguilar Mena <kratsbinovish@gmail.com>
 *                         Guido Giuntoli <gagiuntoli@gmail.com>
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <iomanip>

#include <ctime>
#include <cassert>

#include "micro.hpp"

using namespace std;

#define D_EPS 0.01

int main (int argc, char *argv[])
{
	const int dim = 3;
	if (argc < 5) {
		cerr << "Usage: " << argv[0] << " nx ny nz dir [steps]" << endl;
		return(1);
	}

	const int nx = atoi(argv[1]);
	const int ny = atoi(argv[2]);
	const int nz = atoi(argv[3]);
	const int dir = atoi(argv[4]);
	const int time_steps = (argc > 5 ? atoi(argv[5]) : 10);  // Optional value
	int size[3];
	size[0] = nx;
	size[1] = ny;
	size[2] = nz;

	int micro_type = 1; // 2 materiales matriz y fibra (3D esfera en matriz)
	double micro_params[5];
	micro_params[0] = 1.0; // lx
	micro_params[1] = 1.0; // ly
	micro_params[2] = 1.0; // lz
	micro_params[3] = 0.1; // grosor capa de abajo
	micro_params[4] = 0.0; // inv_tol

	int mat_types[2] = { 1, 0 }; // dos materiales lineales (type = 0)

	material_t mat_params[2];
	mat_params[0].set(1.0e6, 0.25, 1.0e2, 2.0e5, 1);
	mat_params[1].set(1.0e6, 0.25, 1.0e2, 2.0e5, 1);

	micropp<3> micro(1, size, micro_type, micro_params, mat_params);

	double sig[6], ctan[36];
	double eps[6] = { 0. };

	for (int t = 0; t < time_steps; ++t) {

		cout << "time step = " << t << endl;

		if (t < 20)
			eps[dir] += D_EPS;
		else if (t < 40)
			eps[dir] -= D_EPS;
		else
			eps[dir] += D_EPS;

		micro.set_macro_strain(0, eps);
		micro.homogenize();
		micro.get_macro_stress(0, sig);
		micro.get_macro_ctan(0, ctan);

		micro.update_vars();
		micro.write_info_files ();

		cout << "eps =\t";
		for (int i = 0; i < 6; ++i)
			cout << setw(14) << eps[i] << "\t";
		cout << endl;

		cout << "sig =\t";
		for (int i = 0; i < 6; ++i)
			cout << setw(14) << sig[i] << "\t";
		cout << endl;

		cout << "ctan =\t";
		for (int i = 0; i < 6; ++i)
			cout << setw(14) << ctan[i] << "\t";
		cout << endl;

		cout << endl;
		micro.output (t, 0);
	}
	return 0;
}
