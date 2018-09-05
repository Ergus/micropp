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

	const int tnvoi = nvoi;

	double *lstrain = (double *) rrl_malloc(tnvoi * sizeof(double));
	memcpy(lstrain, macro_strain, tnvoi * sizeof(double));

	double *tm_strain = gp_list[gp_id].macro_strain;

	#pragma oss task out(tm_strain[0;tnvoi]) in(lstrain[0;tnvoi]) //firstprivate(macro_strain[0;tnvoi])
	memcpy(tm_strain, lstrain, tnvoi * sizeof(double));

	#pragma oss task in(gp_list[gp_id])
	gp_list[gp_id].print_strain();

	#pragma oss taskwait
	rrl_free(lstrain);
}


template <int tdim>
void micropp<tdim>::get_macro_stress(const int gp_id,
									 double *macro_stress) const
{
	assert(gp_id < ngp);
	assert(ngp > 0);

	const int tnvoi = nvoi;
	const double *tm_stress = gp_list[gp_id].macro_stress;

	#pragma oss task in(gp_list[gp_id])
	gp_list[gp_id].print_stress();

	#pragma oss task in(tm_stress[0;tnvoi]) if(0)
	memcpy(macro_stress, tm_stress, tnvoi * sizeof(double));
	#pragma oss taskwait
}


template <int tdim>
void micropp<tdim>::get_macro_ctan(const int gp_id, double *macro_ctan) const
{
	assert(gp_id < ngp);
	assert(ngp > 0);

	const int tnvoi2 = nvoi * nvoi;
	const double *tm_ctan = gp_list[gp_id].macro_ctan;

	#pragma oss task in(gp_list[gp_id])
	gp_list[gp_id].print_ctan();

	#pragma oss task in(tm_ctan[0; tnvoi2]) if(0)
	memcpy(macro_ctan, tm_ctan, tnvoi2 * sizeof(double));
	#pragma oss taskwait
}

template <int tdim>
void micropp<tdim>::homogenize()
{
	INST_START;

	const int tnvoi = nvoi;

	for (int gp = 0; gp < ngp; ++gp) {
		gp_t<tdim> * const gp_ptr = &gp_list[gp];

		int *ell_cols_ptr = ell_cols;
		const int ell_cols_size_tmp = ell_cols_size;

		material_t *material_ptr = material_list;
		const int numMaterials_tmp = numMaterials;

		int *elem_type_ptr = elem_type;
		const int nelem_tmp = nelem;

		double *tv_k = &(dint_vars_k[num_int_vars * gp]);
		double *tu_k = &(du_k[nndim *gp]);
		const int nndim_tmp = nndim;
		const int num_int_vars_tmp = num_int_vars;

		printf("%d gp_ptr = %p tv_k = %p tu_k = %p\n",
		       gp, gp_ptr, tv_k, tu_k);

		#pragma oss task weakin(ell_cols_ptr[0; ell_cols_size_tmp]) \
			weakin(material_ptr[0; numMaterials_tmp]) \
			weakin(elem_type_ptr[0; nelem_tmp]) \
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

// Explicit instantiation
template class micropp<2>;
template class micropp<3>;
