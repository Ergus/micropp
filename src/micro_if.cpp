/*
 *  This source code is part of MicroPP: a finite element library
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

#include <cmath>

#include "instrument.hpp"
#include "micro.hpp"

using namespace std;

bool micropp::plastic_law_if(const material_t *material, const double eps[6],
                             const double eps_p_old[6], double alpha_old,
                             double *_dl, double _normal[6], double _s_trial[6]) const
{
	if (dim == 3) {
		/* Calculates _dl, _normal and _s_trial */

		double eps_dev[6], eps_p_dev_1[6];

		get_dev_tensor_3D(eps_p_old, eps_p_dev_1);
		get_dev_tensor_3D(eps, eps_dev);

		for (int i = 0; i < 3; ++i)
			_s_trial[i] = 2 * material->mu * (eps_dev[i] - eps_p_dev_1[i]);

		for (int i = 3; i < 6; ++i)
			_s_trial[i] = material->mu * (eps_dev[i] - eps_p_dev_1[i]);

		double tmp = 0.0;
		for (int i = 0; i < 6; ++i)
			tmp += _s_trial[i] * _s_trial[i];
		double s_norm = sqrt(tmp);

		double f_trial = s_norm -
			sqrt(2. / 3.) * (material->Sy + material->Ka * alpha_old);

		if (f_trial > 0 && material->plasticity) {

			for (int i = 0; i < 6; ++i)
				_normal[i] = _s_trial[i] / s_norm;
			*_dl = f_trial / (2. * material->mu *
			                  (1. + (0. * material->Ka) / (3. * material->mu)));
			return true;

		} else {

			memset(_normal, 0, 6 * sizeof(double));
			*_dl = 0;
			return false;
		}
	} else if (dim == 2){
		return false;
	}

	perror("Invalid dimension value\n");
	abort();
}


bool micropp::plastic_evolute_if(const material_t *material,
                                 const double eps[6],
								 const double eps_p_old[6],
								 double alpha_old, 
								 double *eps_p_new,
								 double *alpha_new) const
{
	if (dim == 3) {
		double dl, normal[6], s_trial[6];
		bool nl_flag = plastic_law_if(material, eps, eps_p_old, alpha_old,
		                              &dl, normal, s_trial);

		for (int i = 0; i < 6; ++i)
			eps_p_new[i] = eps_p_old[i] + dl * normal[i];

		*alpha_new = alpha_old + sqrt(2. / 3.) * dl;

		return nl_flag;
	} else if (dim == 2) {
		return false;
	}

	perror("Invalid dimension value\n");
	abort();
}


void micropp::isolin_get_stress_if(const material_t *material,
                                      const double eps[6],
                                      double stress[6]) const
{
	if (dim == 3) {
		for (int i = 0; i < 3; ++i)
			stress[i] = material->lambda * (eps[0] + eps[1] + eps[2]) \
				+ 2 * material->mu * eps[i];

		for (int i = 3; i < 6; ++i)
			stress[i] = material->mu * eps[i];

	} else if (dim == 2) {
		const double E = material->E;
		const double nu = material->nu;

		double ctan[3][3] = { { (1 - nu),     	nu,                0 },
		                      {       nu, (1 - nu),                0 },
		                      {        0,        0, (1 - 2 * nu) / 2 } };

		for (int i = 0; i < nvoi; i++)
			for (int j = 0; j < nvoi; j++)
				ctan[i][j] *= E / ((1 + nu) * (1 - 2 * nu));

		memset(stress, 0, nvoi * sizeof(double));
		for (int i = 0; i < nvoi; i++)
			for (int j = 0; j < nvoi; j++)
				stress[i] += ctan[i][j] * eps[j];
	} else {

	perror("Invalid dimension value\n");
	abort();
	}

}
