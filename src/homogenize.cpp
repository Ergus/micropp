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

void micropp::set_macro_strain(const int gp_id,
                               const double *macro_strain)
{
	assert(gp_id < ngp);
	assert(ngp > 0);
	memcpy(gp_list[gp_id].macro_strain, macro_strain, nvoi * sizeof(double));
}


void micropp::get_macro_stress(const int gp_id,
                               double *macro_stress) const
{
	assert(gp_id < ngp);
	assert(ngp > 0);
	memcpy(macro_stress, gp_list[gp_id].macro_stress, nvoi * sizeof(double));
}


void micropp::get_macro_ctan(const int gp_id, double *macro_ctan) const
{
	assert(gp_id < ngp);
	assert(ngp > 0);
	memcpy(macro_ctan, gp_list[gp_id].macro_ctan, nvoi * nvoi * sizeof(double));
}

void micropp::homogenize()
{
	const int tnvoi = nvoi;

	for (int gp = 0; gp < ngp; ++gp) {
		gp_t * const gp_ptr = &gp_list[gp];

		int *ell_cols_ptr = ell_cols;
		int ell_cols_size_tmp = ell_cols_size;

		material_t *material_ptr = material_list;
		int numMaterials_tmp = numMaterials;

		int *elem_type_ptr = elem_type;
		int nelem_tmp = nelem;

		double *tv_k = &(dint_vars_k[num_int_vars * gp]);
		double *tu_k = &(du_k[nndim *gp]);
		int nndim_tmp = nndim;
		int num_int_vars_tmp = num_int_vars;

		printf("%d %p %p\n", gp, gp_ptr, gp_ptr->u_k);

		#pragma oss task in(ell_cols_ptr[0; ell_cols_size_tmp]) \
			in(material_ptr[0; numMaterials_tmp]) \
			in(elem_type_ptr[0; nelem_tmp]) \
			 \
			inout(gp_ptr[0]) \
			weakinout(tu_k[0; nndim_tmp]) \
			weakinout(tv_k[0; num_int_vars_tmp])
		homogenize_weak_task(*this, tnvoi,
		                     ell_cols_ptr, ell_cols_size_tmp,
		                     material_ptr, numMaterials_tmp,
		                     elem_type_ptr, nelem_tmp,
		                     gp_ptr, nndim_tmp, num_int_vars_tmp);
	}
	#pragma oss taskwait
}

