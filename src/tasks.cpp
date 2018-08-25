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

#include "tasks.hpp"

void homogenize_conditional_task(micropp self, int nvoi,
                                 int *ell_cols, const int ell_cols_size,
                                 const material_t *material_list, const int numMaterials,
                                 int *elem_type, int nelem,
                                 gp_t *gp_ptr,
                                 int nndim, int num_int_vars,
                                 const bool allocated)
{
	dprintf("Nanos cluster %d/%d \n", get_node_id(), get_nodes_nr());

	double *u_n = (double *) malloc(nndim * sizeof(double));
	double *vnew = (double *) malloc(num_int_vars * sizeof(double));

	double *vold = nullptr, *aux_old = nullptr;

	if (!allocated) {
		aux_old = (double *) calloc(num_int_vars, sizeof(double));
		vold = aux_old;
	} else {
		vold = gp_ptr->int_vars_k;
	}

	double *u_aux = (double *) malloc(nndim * sizeof(double));
	double *du_aux = (double *) malloc(nndim * sizeof(double));
	double *b = (double *) malloc(nndim * sizeof(double));

	const int ns[3] = { self.nx, self.ny, self.nz };
	ell_matrix A;
	ell_init(&A, ell_cols, self.dim, self.dim, ns, CG_MIN_ERR, CG_MAX_ITS);

	// SIGMA (1 Newton-Raphson)
	memcpy(u_n, gp_ptr->u_k, nndim * sizeof(double));

	double nr_err;
	int nr_its = self.newton_raphson(gp_ptr->macro_strain, &A, u_n,
	                                 b, du_aux, vold, &nr_err);
	gp_ptr->nr_its[0] = nr_its;
	gp_ptr->nr_err[0] = nr_err;

	self.calc_ave_stress(u_n, vold, gp_ptr->macro_stress);

	bool nl_flag = self.calc_vars_new(u_n, vold, vnew);

	if (nl_flag) {
		if (!allocated) {
			gp_ptr->allocate(num_int_vars);
			memcpy(gp_ptr->int_vars_k, vnew, num_int_vars * sizeof(double));
		}
	}

	// CTAN (6 Newton-Raphson's)
	memcpy(u_aux, u_n, nndim * sizeof(double));
	double eps_1[6], sig_0[6], sig_1[6];
	memcpy(sig_0, gp_ptr->macro_stress, self.nvoi * sizeof(double));

	for (int i = 0; i < self.nvoi; ++i) {

		memcpy(eps_1, gp_ptr->macro_strain, self.nvoi * sizeof(double));
		eps_1[i] += D_EPS_CTAN_AVE;

		nr_its = self.newton_raphson(eps_1, &A, u_aux, b, du_aux, vold, &nr_err);
		gp_ptr->nr_its[i + 1] = nr_its;
		gp_ptr->nr_err[i + 1] = nr_err;

		self.calc_ave_stress(u_aux, vold, sig_1);

		for (int v = 0; v < self.nvoi; ++v)
			gp_ptr->macro_ctan[v * self.nvoi + i] = (sig_1[v] - sig_0[v]) / D_EPS_CTAN_AVE;

	}

	memcpy(gp_ptr->u_k, u_n, nndim * sizeof(double));

	if (gp_ptr->allocated)
		memcpy(gp_ptr->int_vars_k, vnew, num_int_vars * sizeof(double));

	free(aux_old);

	ell_free(&A);
	free(u_aux);
	free(du_aux);
	free(b);
}


void homogenize_weak_task(micropp self, int nvoi,
                          int *ell_cols, const int ell_cols_size,
                          const material_t *material_list, const int numMaterials,
                          int *elem_type, int nelem,
                          gp_t *gp_ptr, int nndim, int num_int_vars)
{

	if (gp_ptr->is_linear(self.ctan_lin, self.inv_tol, -1.0e10)
	    && (!gp_ptr->allocated)) {

		/* This is a risky optimization and should be used with extreme
		 * caution setting the variable <inv_tol> to a right and small value
		 */

		for (int i = 0; i < nvoi; ++i) {
			gp_ptr->macro_stress[i] = 0.0;
			for (int j = 0; j < nvoi; ++j)
				gp_ptr->macro_stress[i] += self.ctan_lin[i * nvoi + j]
					*gp_ptr->macro_strain[j];
		}

		memcpy(gp_ptr->macro_ctan, self.ctan_lin, nvoi * nvoi * sizeof(double));
		memset(gp_ptr->nr_its, 0, (1 + nvoi) * sizeof(double));
		memset(gp_ptr->nr_err, 0, (1 + nvoi) * sizeof(double));

	} else {

		double *tv_k = gp_ptr->int_vars_k;
		double *tu_k = gp_ptr->u_k;

		if (gp_ptr->allocated) {
			#pragma oss task in(ell_cols[0; ell_cols_size]) \
				in(material_list[0; numMaterials]) \
				in(elem_type[0; nelem]) \
                                                                        \
				inout(gp_ptr) \
				inout(tu_k[0; nndim]) \
				inout(tv_k[0; num_int_vars])
			{
				printf("CALLER: %p (%p) ", gp_ptr, gp_ptr->u_k);
				gp_ptr->print();
			homogenize_conditional_task(self, nvoi,
			                            ell_cols, ell_cols_size,
			                            material_list, numMaterials,
			                            elem_type, nelem,
			                            gp_ptr, nndim, num_int_vars,
			                            true);


			}
		} else {

			#pragma oss task in(ell_cols[0; ell_cols_size]) \
				in(material_list[0; numMaterials]) \
				in(elem_type[0; nelem]) \
                                                                        \
				inout(gp_ptr) \
				inout(tu_k[0; nndim]) \
				inout(tv_k[0; num_int_vars])
			{
				printf("CALLER: %p (%p) ", gp_ptr, gp_ptr->u_k);
				gp_ptr->print();
			homogenize_conditional_task(self, nvoi,
			                            ell_cols, ell_cols_size,
			                            material_list, numMaterials,
			                            elem_type, nelem,
			                            gp_ptr, nndim, num_int_vars,
			                            false);
			}
		}
	}
}

