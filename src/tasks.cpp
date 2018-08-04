/*
 * This source code is part of MicroPP: a finite element library
 * to solve microstructural problems for composite materials.
 *
 * Copyright (C) - 2018 - Jimmy Aguilar Mena <kratsbinovish@gmail.com>
 *                        Guido Giuntoli <gagiuntoli@gmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "micro.hpp"
#include "tasks.hpp"

template <int tdim>
void micropp<tdim>::homogenize_conditional(data in_data,
                                           int *ell_cols, int ell_cols_size,
                                           material_t *material_list, int numMaterials,
                                           int *elem_type, int nelem,
                                           gp_t<tdim> *gp_ptr,
                                           double *u_k, double *u_n, int nndim,
                                           bool allocated, double *vars_n_old,
                                           double *vars_k_new, int num_int_vars)
{
	double *vold = vars_n_old, *vnew = vars_n_old;
	double *aux_old = nullptr, *aux_new = nullptr;

	if (!allocated) {
		aux_old = (double *) calloc(num_int_vars, sizeof(double));
		aux_new = (double *) malloc(num_int_vars * sizeof(double));

		vold = aux_old;
		vnew = aux_new;
	}

	double *u_aux = (double *) malloc(nndim * sizeof(double));
	double *du_aux = (double *) malloc(nndim * sizeof(double));
	double *b = (double *) malloc(nndim * sizeof(double));
	const int ns[3] = { nx, ny, nz };
	ell_matrix A;
	ell_init(&A, ell_cols, dim, dim, ns, CG_MIN_ERR, CG_MAX_ITS);

	// SIGMA (1 Newton-Raphson)
	memcpy(u_k, u_n, nndim * sizeof(double));

	double nr_err;
	int nr_its = newton_raphson(gp_ptr->macro_strain, &A, gp_ptr->u_k,
	                            b, du_aux, vold, &nr_err);
	gp_ptr->nr_its[0] = nr_its;
	gp_ptr->nr_err[0] = nr_err;

	calc_ave_stress(u_k, vold, gp_ptr->macro_stress);

	bool nl_flag = calc_vars_new(gp_ptr->u_k, vold, vnew);

	if (nl_flag && !allocated) {
			allocate_gp(gp_ptr, gp_ptr->int_vars_n, num_int_vars);
			memcpy(gp_ptr->int_vars_k, vnew,
			       num_int_vars * sizeof(double));
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
			gp_ptr->macro_ctan[v * nvoi + i] = (sig_1[v] - sig_0[v]) / D_EPS_CTAN_AVE;

	}

	ell_free(&A);
	free(u_aux);
	free(du_aux);
	free(b);
	free(aux_old);
	free(aux_new);
}


template <int tdim>
void micropp<tdim>::homogenize_task(data in_data,
                                    int *ell_cols, int ell_cols_size,
                                    material_t *material_list, int numMaterials,
                                    int *elem_type, int nelem,
                                    gp_t<tdim> *gp_ptr,
                                    double *u_k, double *u_n, int nndim,
                                    double *vars_n_old, double *vars_k_new,
                                    int num_int_vars)
{
	micropp<tdim> self(&in_data, gp_ptr);

	if (gp_ptr->is_linear(self.ctan_lin, self.inv_tol) && (!gp_ptr->allocated)) {

		for (int i = 0; i < nvoi; ++i) {
			gp_ptr->macro_stress[i] = 0.0;
			for (int j = 0; j < nvoi; ++j)
				gp_ptr->macro_stress[i] += self.ctan_lin[i * nvoi + j]
					* gp_ptr->macro_strain[j];
		}
		memcpy(gp_ptr->macro_ctan, self.ctan_lin, nvoi * nvoi * sizeof(double));
		memset(gp_ptr->nr_its, 0, (1 + nvoi) * sizeof(double));
		memset(gp_ptr->nr_err, 0, (1 + nvoi) * sizeof(double));

	} else {

		bool allocated = gp_ptr->allocated;

		if (allocated) {
			#pragma oss task in(ell_cols[0; ell_cols_size]) \
				in(material_list[0; numMaterials]) \
				in(elem_type[0; nelem]) \
				firstprivate(in_data) \
				inout(gp_ptr[0]) \
				out(u_k[0; nndim]) \
				inout(u_n[0; nndim]) \
				inout(vars_n_old[0; num_int_vars]) \
				inout(vars_k_new[0; num_int_vars])
			homogenize_conditional(in_data,
				ell_cols, ell_cols_size,
				material_list, numMaterials,
				elem_type, nelem,
				gp_ptr,
				u_k, u_n, nndim,
				true, vars_n_old, vars_k_new, num_int_vars);

		} else {
			#pragma oss task in(ell_cols[0; ell_cols_size]) \
				in(material_list[0; numMaterials]) \
				in(elem_type[0; nelem]) \
				 \
				inout(gp_ptr[0]) \
				out(u_k[0; nndim]) \
				inout(u_n[0; nndim]) \
				out(vars_n_old[0; num_int_vars]) \
				out(vars_k_new[0; num_int_vars])
			homogenize_conditional(in_data,
			                       ell_cols, ell_cols_size,
			                       material_list, numMaterials,
			                       elem_type, nelem,
			                       gp_ptr,
			                       u_k, u_n, nndim,
			                       false, vars_n_old, vars_k_new, num_int_vars);

		}
	}
}


// Explicit instantiation
template class micropp<2>;
template class micropp<3>;