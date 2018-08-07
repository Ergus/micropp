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

	const int tnvoi = nvoi;

	for (int gp = 0; gp < ngp; ++gp) {
		gp_t<tdim> * const gp_ptr = &gp_list[gp];

		double *tv_n = &(dint_vars_n[num_int_vars * gp]);
		double *tv_k = &(dint_vars_k[num_int_vars * gp]);
		double *tu_n = &(du_n[nndim *gp]);
		double *tu_k = &(du_k[nndim *gp]);

		printf("%d %p %p\n", gp, gp_ptr, gp_ptr->u_k);

		homogenize_weak_task(*this, tnvoi,
		                     ell_cols, ell_cols_size,
		                     material_list, numMaterials,
		                     elem_type, nelem,
		                     gp_ptr,
		                     tu_n, tu_k, nndim,
		                     tv_n, tv_k, num_int_vars);
	}
	#pragma oss taskwait
}

template <int tdim>
void micropp<tdim>::update_vars()
{
	for (int gp = 0; gp < ngp; ++gp) {
		gp_list[gp].update_vars();
	}
}

// Explicit instantiation
template class micropp<2>;
template class micropp<3>;
