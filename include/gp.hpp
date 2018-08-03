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

#ifndef GP_HPP
#define GP_HPP

#include <cassert>
#include <cstdlib>
#include <cmath>

template <int dim>
class gp_t {
		static constexpr int nvoi = dim * (dim + 1) / 2;  // 3, 6
	public:
		double macro_strain[nvoi];
		double macro_stress[nvoi];
		double macro_ctan[nvoi * nvoi];

		bool allocated; // flag for memory optimization

		double *int_vars_n; // vectors for calculations
		double *int_vars_k;
		double *u_n;
		double *u_k;

		int nr_its[nvoi + 1]; // measurements
		double nr_err[nvoi + 1];
		double inv_max;

		gp_t() = delete;

		void init(double *_int_vars_n, double *_int_vars_k,
		          double *_u_n, double *_u_k, int nndim)
		{
			assert(nndim > 0);
			allocated = false;

			int_vars_n = _int_vars_n;
			int_vars_k = _int_vars_k;

			u_n = _u_n;
			u_k = _u_k;

			inv_max = -1.0;

			memset(u_n, 0, nndim * sizeof(double));
		}

		~gp_t()
		{}

		void allocate(const int num_int_vars)
		{
			assert(!allocated);

			allocated = true;
			memset(int_vars_n, 0, num_int_vars * sizeof(double));
		}


		void update_vars()
		{
			double *tmp = int_vars_n;
			int_vars_n = int_vars_k;
			int_vars_k = tmp;

			tmp = u_n;
			u_n = u_k;
			u_k = tmp;
		}


		double get_inv_1(const double *tensor) const
		{
			const double ret = tensor[0] + tensor[1];

			if (dim == 2)
				return ret;

			return ret + tensor[2];
		}


		bool is_linear(double *ctan_lin, double inv_tol)
		{
			double macro_stress[6] = { 0.0 };

			inv_max = -1.0e10;

			for (int i = 0; i < nvoi; ++i)
				for (int j = 0; j < nvoi; ++j)
					macro_stress[i] += ctan_lin[i * nvoi + j] * macro_strain[j];

			const double inv = get_inv_1(macro_stress);
			if (fabs(inv) > inv_max)
				inv_max = fabs(inv);

			return (fabs(inv) < inv_tol);
		}
};

#endif // GP_HPP
