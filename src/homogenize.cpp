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

#include <cmath>
#include <cassert>

#include "instrument.hpp"
#include "micro.hpp"

template <int tdim>
void micropp<tdim>::set_macro_strain(const int gp_id,
									 const double *macro_strain)
{
	assert(gp_id < ngp);
	assert(ngp > 0);
	memcpy(gp_list[gp_id].macro_strain, macro_strain, nvoi * sizeof(double));
}


template <int tdim>
void micropp<tdim>::get_macro_stress(const int gp_id,
									 double *macro_stress) const
{
	assert(gp_id < ngp);
	assert(ngp > 0);
	memcpy(macro_stress, gp_list[gp_id].macro_stress, nvoi * sizeof(double));
}


template <int tdim>
void micropp<tdim>::get_macro_ctan(const int gp_id, double *macro_ctan) const
{
	assert(gp_id < ngp);
	assert(ngp > 0);
	memcpy(macro_ctan, gp_list[gp_id].macro_ctan, nvoi * nvoi * sizeof(double));
}

template <int tdim>
void micropp<tdim>::homogenize()
{
	INST_START;

	for (int gp = 0; gp < ngp; ++gp) {
		gp_t<tdim> * const gp_ptr = &gp_list[gp];

		if (gp_ptr->is_linear(ctan_lin, inv_tol, -1.0e10) && (!gp_ptr->allocated)) {

			/* This is a risky optimization and should be used with extreme
			 * caution setting the variable <inv_tol> to a right and small value
			 */

			for (int i = 0; i < nvoi; ++i) {
				gp_ptr->macro_stress[i] = 0.0;
				for (int j = 0; j < nvoi; ++j)
					gp_ptr->macro_stress[i] += ctan_lin[i * nvoi + j]
						* gp_ptr->macro_strain[j];
			}
			memcpy(gp_ptr->macro_ctan, ctan_lin, nvoi * nvoi * sizeof(double));
			memset(gp_ptr->nr_its, 0, (1 + nvoi) * sizeof(double));
			memset(gp_ptr->nr_err, 0, (1 + nvoi) * sizeof(double));

		} else {

			double *u_aux = (double *) malloc(nndim * sizeof(double));
			double *du_aux = (double *) malloc(nndim * sizeof(double));
			double *b = (double *) malloc(nndim * sizeof(double));

			const int ns[3] = { nx, ny, nz };
			ell_matrix A;
			ell_init(&A, ell_cols, dim, dim, ns, CG_MIN_ERR, CG_MAX_ITS);

			double *vold, *vnew, *aux_old = nullptr, *aux_new = nullptr;
			if (!gp_ptr->allocated) {

				aux_old = (double *) calloc(num_int_vars, sizeof(double));
				aux_new = (double *) malloc(num_int_vars * sizeof(double));

				vold = aux_old;
				vnew = aux_new;

			} else {
				vold = gp_ptr->int_vars_n;
				vnew = gp_ptr->int_vars_k;
			}

			// SIGMA (1 Newton-Raphson)
			memcpy(gp_ptr->u_k, gp_ptr->u_n, nndim * sizeof(double));

			double nr_err;
			int nr_its = newton_raphson(gp_ptr->macro_strain, &A, gp_ptr->u_k,
			                            b, du_aux, vold, &nr_err);

			gp_ptr->nr_its[0] = nr_its;
			gp_ptr->nr_err[0] = nr_err;

			calc_ave_stress(gp_ptr->u_k, vold, gp_ptr->macro_stress);

			bool nl_flag = calc_vars_new(gp_ptr->u_k, vold, vnew);

			if (nl_flag) {
				if (!gp_ptr->allocated) {
					gp_ptr->allocate(num_int_vars);
					memcpy(gp_ptr->int_vars_k, vnew, num_int_vars * sizeof(double));
				}
			}

			// CTAN (6 Newton-Raphson's)
			memcpy(u_aux, gp_ptr->u_k, nndim * sizeof(double));
			double eps_1[6], sig_0[6], sig_1[6];
			memcpy(sig_0, gp_ptr->macro_stress, nvoi * sizeof(double));

			for (int i = 0; i < nvoi; ++i) {

				memcpy(eps_1, gp_ptr->macro_strain, nvoi * sizeof(double));
				eps_1[i] += D_EPS_CTAN_AVE;

				nr_its = newton_raphson(eps_1, &A, u_aux, b, du_aux, vold, &nr_err);
				gp_ptr->nr_its[i + 1] = nr_its;
				gp_ptr->nr_err[i + 1] = nr_err;

				calc_ave_stress(u_aux, vold, sig_1);

				for (int v = 0; v < nvoi; ++v)
					gp_ptr->macro_ctan[v * nvoi + i] =
						(sig_1[v] - sig_0[v]) / D_EPS_CTAN_AVE;

			}
			free(aux_old);
			free(aux_new);

			ell_free(&A);
			free(u_aux);
			free(du_aux);
			free(b);
		}
	}
}

template <int tdim>
void micropp<tdim>::update_vars()
{
	for (int igp = 0; igp < ngp; ++igp)
		gp_list[igp].update_vars();
}

// Explicit instantiation
template class micropp<2>;
template class micropp<3>;
