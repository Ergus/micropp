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
#include "gp.hpp"


template <int tdim>
void micropp<tdim>::set_macro_strain(const int gp_id,
				     const double *macro_strain)
{
	assert(gp_id >= 0);
	assert(gp_id < ngp);

	double *tpgp_macro_strain = gp_list[gp_id].macro_strain;
	const int tnvoi = nvoi;

	#pragma oss task out(tpgp_macro_strain[0; tnvoi]) in(macro_strain[0; tnvoi])
	memcpy(tpgp_macro_strain, macro_strain, tnvoi * sizeof(double));

}


template <int tdim>
void micropp<tdim>::get_macro_stress(const int gp_id, double *macro_stress) const
{
	INST_START;

	assert(gp_id >= 0);
	assert(gp_id < ngp);

	double *tpgp_macro_stress = gp_list[gp_id].macro_stress;
	const int tnvoi = nvoi;

	#pragma oss task in(tpgp_macro_stress[0; tnvoi]) out(macro_stress[0; tnvoi])
	memcpy(macro_stress, tpgp_macro_stress, tnvoi * sizeof(double));
}


template <int tdim>
void micropp<tdim>::get_macro_ctan(const int gp_id, double *macro_ctan) const
{
	INST_START;

	assert(gp_id >= 0);
	assert(gp_id < ngp);

	double *tpgp_macro_ctan = gp_list[gp_id].macro_ctan;
	const int tnvoi2 = nvoi * nvoi;

	#pragma oss task in(tpgp_macro_ctan[0; tnvoi2]) out(macro_ctan[0; tnvoi2])
	memcpy(macro_ctan, tpgp_macro_ctan, tnvoi2 * sizeof(double));
}


template <int tdim>
void micropp<tdim>::homogenize()
{
	INST_START;

	#pragma omp parallel for schedule(dynamic,1)
	for (int igp = 0; igp < ngp; ++igp) {

		const int lnvoi = nvoi;
		int *tpell_cols = ell_cols;
		const int tell_cols_size = ell_cols_size;

		material_t *tpmaterial = material_list;
		const int tnumMaterials = numMaterials;

		int *tpelem_type = elem_type;
		const int tnelem = nelem;

		gp_t<tdim> *gp_ptr = &gp_list[igp];
		const int tnndim = nndim;
		const int tnum_int_vars = num_int_vars;

		double *tpu_k = gp_ptr->u_k;
		double *tpint_vars_n = gp_ptr->int_vars_k;

		#pragma oss task in(this[0])				\
			weakin(tpell_cols[0; tell_cols_size])		\
			weakin(tpmaterial[0; tnumMaterials])		\
			weakin(tpelem_type[0; tnelem])			\
									\
			inout(gp_ptr[0])				\
			weakinout(tpu_k[0; tnndim])			\
			weakinout(tpint_vars_n[0; tnum_int_vars])
		{

			if (gp_ptr->allocated) {
				#pragma oss task in(this[0])		\
					in(tpell_cols[0; tell_cols_size]) \
					in(tpmaterial[0; tnumMaterials]) \
					in(tpelem_type[0; tnelem])	\
									\
					inout(gp_ptr[0])		\
					inout(tpu_k[0; tnndim])		\
					inout(tpint_vars_n[0; tnum_int_vars])
				homogenize_conditional_task(lnvoi,
							    tpell_cols, tell_cols_size,
							    tpmaterial, tnumMaterials,
							    tpelem_type, tnelem,
							    gp_ptr, tnndim, tnum_int_vars);
			} else {
				#pragma oss task in(this[0])		\
					in(tpell_cols[0; tell_cols_size]) \
					in(tpmaterial[0; tnumMaterials]) \
					in(tpelem_type[0; tnelem])	\
									\
					inout(gp_ptr[0])		\
					out(tpu_k[0; tnndim])		\
					out(tpint_vars_n[0; tnum_int_vars])
				homogenize_conditional_task(lnvoi,
							    tpell_cols, tell_cols_size,
							    tpmaterial, tnumMaterials,
							    tpelem_type, tnelem,
							    gp_ptr, tnndim, tnum_int_vars);
			}
		}
	}
}


template <int tdim>
void micropp<tdim>::homogenize_conditional_task(const int nvoi,
						int *ell_cols,
						const int ell_cols_size,
						const material_t *material_list,
						const int numMaterials,
						int *elem_type, int nelem,
						gp_t<tdim> *gp_ptr,
						int nndim, int num_int_vars)
{
	newton_t newton;
	newton.max_its = NR_MAX_ITS;
	newton.max_tol = NR_MAX_TOL;
	newton.rel_tol = NR_REL_TOL;

	int ierr;

	double *u = (double *) malloc(nndim * sizeof(double));
	double *du = (double *) malloc(nndim * sizeof(double));
	double *b = (double *) malloc(nndim * sizeof(double));
	const int ns[3] = { nx, ny, nz };

	ell_matrix A;
	ell_init(&A, ell_cols, dim, dim, ns, CG_MIN_ERR, CG_REL_ERR, CG_MAX_ITS);

	int newton_its, solver_its[NR_MAX_ITS];
	double newton_err[NR_MAX_ITS], solver_err[NR_MAX_ITS];

	gp_ptr->cost = 0;

	double *vars_new,  *vars_new_aux = nullptr;

	if (gp_ptr->allocated)
		vars_new = gp_ptr->int_vars_k;
	else
		vars_new = vars_new_aux = (double *) malloc(num_int_vars * sizeof(double));

	// SIGMA 1 Newton-Raphson
	memcpy(u, gp_ptr->u_n, nndim * sizeof(double));

	newton_raphson(&A, b, u, du,
		       gp_ptr->allocated, gp_ptr->macro_strain,
		       gp_ptr->int_vars_n, &newton);

	memcpy(gp_ptr->u_k, u, nndim * sizeof(double));
	gp_ptr->newton = newton;

	for (int i = 0; i < newton.its; ++i)
		gp_ptr->cost += newton.solver_its[i];

	if (coupling == ONE_WAY) {

		double *stress = gp_ptr->macro_stress;
		double *strain = gp_ptr->macro_strain;
		memset (stress, 0.0, nvoi * sizeof(double));
		for (int i = 0; i < nvoi; ++i) {
			for (int j = 0; j < nvoi; ++j)
				stress[i] += ctan_lin[i * nvoi + j] * strain[j];
		}

	} else if (coupling == FULL || coupling == NO_COUPLING) {

		calc_ave_stress(gp_ptr->u_k, gp_ptr->int_vars_n, gp_ptr->macro_stress);
		filter(gp_ptr->macro_stress, nvoi, FILTER_REL_TOL);
	}

	bool nl_flag = calc_vars_new(gp_ptr->u_k, gp_ptr->int_vars_n, vars_new, &f_trial_max);

	if (nl_flag && (gp_ptr->allocated == false)) {
		gp_ptr->allocate();
		memcpy(gp_ptr->int_vars_k, vars_new, num_int_vars * sizeof(double));
	}

	if (gp_ptr->allocated && coupling == FULL) {
		// CTAN 3/6 Newton-Raphsons in 2D/3D
		double eps_1[6], sig_0[6], sig_1[6];

		memcpy(u, gp_ptr->u_k, nndim * sizeof(double));
		memcpy(sig_0, gp_ptr->macro_stress, nvoi * sizeof(double));

		for (int i = 0; i < nvoi; ++i) {

			memcpy(eps_1, gp_ptr->macro_strain, nvoi * sizeof(double));
			eps_1[i] += D_EPS_CTAN_AVE;

			newton_raphson(&A, b, u, du,
				       true, eps_1, gp_ptr->int_vars_n, &newton);

			calc_ave_stress(u, gp_ptr->int_vars_n, sig_1);

			for (int v = 0; v < nvoi; ++v)
				gp_ptr->macro_ctan[v * nvoi + i] =
					(sig_1[v] - sig_0[v]) / D_EPS_CTAN_AVE;

		}

		filter(gp_ptr->macro_ctan, nvoi * nvoi, FILTER_REL_TOL);
	}

	free(u);
	free(du);
	free(b);

	if (vars_new_aux)
		free(vars_new_aux);
}



template <int tdim>
void micropp<tdim>::update_vars()
{
	INST_START;

	gp_t<tdim> *tpgp = gp_list;

	for (int igp = 0; igp < ngp; ++igp){
		#pragma oss task inout(tpgp[igp])
		tpgp[igp].update_vars();
	}
}




template class micropp<2>;
template class micropp<3>;
