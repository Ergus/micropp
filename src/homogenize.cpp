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
#include "tasks.hpp"

template <int tdim>
void micropp<tdim>::set_macro_strain(const int gp_id, const double *const macro_strain)
{
	assert(gp_id < ngp);
	assert(ngp > 0);
	memcpy(gp_list[gp_id].macro_strain, macro_strain, nvoi * sizeof(double));
}


template <int tdim>
void micropp<tdim>::get_macro_stress(const int gp_id, double *macro_stress) const
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

	for (int igp = 0; igp < ngp; ++igp) {
		gp_t<tdim> * const gp_ptr = &gp_list[igp];

		homogenize_task(*this,
		                ell_cols, ell_cols_size,
		                material_list, numMaterials,
		                elem_type, nelem,
		                &gp_list[igp],
		                gp_ptr->u_k, gp_ptr->u_n, nndim,
		                gp_ptr->int_vars_n, gp_ptr->int_vars_k, num_int_vars);
	}
}

template <int tdim>
void micropp<tdim>::update_vars()
{
	for (int igp = 0; igp < ngp; ++igp) {
		gp_list[igp].update_vars();
	}
}

// Explicit instantiation
template class micropp<2>;
template class micropp<3>;


