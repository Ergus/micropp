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

template <int tdim>
void homogenize_conditional_task(struct data self_data, int nvoi,
                            int *ell_cols, const int ell_cols_size,
                            const material_t *material_list, const int numMaterials,
                            int *elem_type, int nelem,
                            gp_t<tdim> *gp_ptr,
                            double *u_n, double *u_k, int nndim,
                            const bool allocated,
                            double *vars_n_old, double *vars_k_new, int num_int_vars)
{
	dprintf("Nanos cluster %d/%d \n", get_node_id(), get_nodes_nr());

	printf("gp: %p %p != %p\n", gp_ptr, gp_ptr->u_k, u_k);

	micropp<tdim> self(&self_data, gp_ptr);

	assert(gp_ptr->u_k == u_k);
	assert(gp_ptr->u_n == u_n);

	double *vold = vars_n_old, *vnew = vars_k_new;
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

	const int ns[3] = { self.nx, self.ny, self.nz };
	ell_matrix A;
	ell_init(&A, ell_cols, self.dim, self.dim, ns, CG_MIN_ERR, CG_MAX_ITS);

	// SIGMA (1 Newton-Raphson)
	memcpy(gp_ptr->u_k, gp_ptr->u_n, nndim * sizeof(double));

	double nr_err;
	int nr_its = self.newton_raphson(gp_ptr->macro_strain, &A, gp_ptr->u_k,
	                                 b, du_aux, vold, &nr_err);
	gp_ptr->nr_its[0] = nr_its;
	gp_ptr->nr_err[0] = nr_err;

	self.calc_ave_stress(gp_ptr->u_k, vold, gp_ptr->macro_stress);

	bool nl_flag = self.calc_vars_new(gp_ptr->u_k, vold, vnew);

	if (nl_flag) {
		if (!allocated) {
			gp_ptr->allocate(num_int_vars);
			memcpy(gp_ptr->int_vars_k, vnew, num_int_vars * sizeof(double));
		}
	}

	// CTAN (6 Newton-Raphson's)
	memcpy(u_aux, gp_ptr->u_k, nndim * sizeof(double));
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

	free(aux_old);
	free(aux_new);

	ell_free(&A);
	free(u_aux);
	free(du_aux);
	free(b);
}


template <int tdim>
void homogenize_weak_task(struct data self, int nvoi,
                     int *ell_cols, const int ell_cols_size,
                     const material_t *material_list, const int numMaterials,
                     int *elem_type, int nelem,
                     gp_t<tdim> *gp_ptr,
                     double *u_n, double *u_k, int nndim,
                     double *vars_n_old, double *vars_k_new, int num_int_vars)
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

		if (gp_ptr->allocated) {
			#pragma oss task in(ell_cols[0; ell_cols_size]) \
				in(material_list[0; numMaterials]) \
				in(elem_type[0; nelem]) \
                                                                        \
				inout(gp_ptr[0]) \
				inout(u_n[0; nndim]) \
				inout(u_k[0; nndim]) \
				inout(vars_n_old[0; num_int_vars]) \
				inout(vars_k_new[0; num_int_vars])
			{
				printf("CALLER: %p (%p) ", gp_ptr, u_k);
				gp_ptr->print();
			homogenize_conditional_task<tdim>(self, nvoi,
			                                  ell_cols, ell_cols_size,
			                                  material_list, numMaterials,
			                                  elem_type, nelem,
			                                  gp_ptr,
			                                  u_n, u_k, nndim,
			                                  true,
			                                  vars_n_old, vars_k_new, num_int_vars);


			}
		} else {

			#pragma oss task in(ell_cols[0; ell_cols_size]) \
				in(material_list[0; numMaterials]) \
				in(elem_type[0; nelem]) \
                                                                        \
				inout(gp_ptr[0]) \
				inout(u_n[0; nndim]) \
				inout(u_k[0; nndim]) \
				inout(vars_n_old[0; num_int_vars]) \
				inout(vars_k_new[0; num_int_vars])
			{
				printf("CALLER: %p (%p) ", gp_ptr, u_k);
				gp_ptr->print();
			homogenize_conditional_task<tdim>(self, nvoi,
			                                  ell_cols, ell_cols_size,
			                                  material_list, numMaterials,
			                                  elem_type, nelem,
			                                  gp_ptr,
			                                  u_n, u_k, nndim,
			                                  false,
			                                  vars_n_old, vars_k_new, num_int_vars);
			}
		}
	}
}

// Explicit instantiation
template class micropp<2>;
template class micropp<3>;

template
void homogenize_conditional_task<2>(struct data self, int nvoi,
                                 int *ell_cols, const int ell_cols_size,
                                 const material_t *material_list, const int numMaterials,
                                 int *elem_type, int nelem,
                                 gp_t<2> *gp_ptr,
                                 double *u_k, double *u_n, int nndim,
                                 const bool allocated, double *vars_n_old,
                                 double *vars_k_new, int num_int_vars);

template
void homogenize_conditional_task<3>(struct data self, int nvoi,
                                 int *ell_cols, const int ell_cols_size,
                                 const material_t *material_list, const int numMaterials,
                                 int *elem_type, int nelem,
                                 gp_t<3> *gp_ptr,
                                 double *u_k, double *u_n, int nndim,
                                 const bool allocated, double *vars_n_old,
                                 double *vars_k_new, int num_int_vars);

template
void homogenize_weak_task<2>(data self,int nvoi,
                          int *ell_cols, const int ell_cols_size,
                          const material_t *material_list, const int numMaterials,
                          int *elem_type, int nelem,
                          gp_t<2> *gp_ptr,
                          double *u_k, double *u_n, int nndim,
                          double *vars_n_old, double *vars_k_new, int num_int_vars);

template
void homogenize_weak_task<3>(data self,int nvoi,
                          int *ell_cols, const int ell_cols_size,
                          const material_t *material_list, const int numMaterials,
                          int *elem_type, int nelem,
                          gp_t<3> *gp_ptr,
                          double *u_k, double *u_n, int nndim,
                          double *vars_n_old, double *vars_k_new, int num_int_vars);
