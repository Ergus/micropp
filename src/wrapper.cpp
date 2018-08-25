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

#include "micro.hpp"

extern "C" {

	void material_set_(material_t *in, const double *E,
	                   const double *nu, const double *Ka,
	                   const double *Sy, const int *type)
	{
		in->set(*E, *nu, *Ka, *Sy, *type);
	}

	//  Micro type functions
	micropp *init_(int *dim, int *ngp, const int size[3], const int *micro_type,
	                   const double *micro_params, const material_t *materials)
	{
		return new micropp(*dim ,*ngp, size, *micro_type, micro_params, materials);
	}

	void free_(micropp **in)
	{
		delete (*in);
	}

	void get_nl_flag_(const micropp **self, int *gp_id, int *nl_flag)
	{
		(*nl_flag) = (*self)->get_nl_flag(*gp_id);
	}

	void set_macro_strain_(micropp **self, const int *gp_id,
	                        const double *macro_strain)
	{
		(*self)->set_macro_strain(*gp_id, macro_strain);
	}

	void get_macro_stress_(const micropp **self,
	                        const int *gp_id, double *macro_stress)
	{
		(*self)->get_macro_stress(*gp_id, macro_stress);
	}

	void get_macro_ctan_(const micropp **self, const int *gp_id,
	                     double *macro_ctan)
	{
		(*self)->get_macro_ctan(*gp_id, macro_ctan);
	}

	void homogenize3_(micropp **self)
	{
		(*self)->homogenize();
	}

	void update_vars3_(micropp **self)
	{
		(*self)->update_vars();
	}

	void output3_(micropp **self, int *tstep, int *gp_id)
	{
		(*self)->output(*tstep, *gp_id);
	}

	void write_info_files3_(micropp **self)
	{
		(*self)->write_info_files();
	}

	void print_info3_(micropp **self)
	{
		printf("ptr2 %p\n", self);
		(*self)->print_info();
	}
}
