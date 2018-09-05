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

#include <iostream>
#include <iomanip>

using namespace std;

template <int dim>
class gp_t {
		static constexpr int nvoi = dim * (dim + 1) / 2;  // 3, 6
	public:
		double macro_strain[nvoi];
		double macro_stress[nvoi];
		double macro_ctan[nvoi * nvoi];

		bool allocated; // flag for memory optimization

		double *int_vars_k;
		double *u_k;

		int nr_its[nvoi + 1]; // measurements
		double nr_err[nvoi + 1];
		double inv_max;

		gp_t() = delete;
		~gp_t()	{}

		void init(double *_int_vars_k, double *_u_k, int nndim)
		{
			assert(nndim > 0);
			allocated = false;
			inv_max = -1.0;

			int_vars_k = _int_vars_k;
			u_k = _u_k;
		}


		void allocate(const int num_int_vars)
		{
			assert(!allocated);
			allocated = true;
		}

		bool is_linear(const double *ctan_lin, const double _inv_tol,
		               double _inv_max)
		{
			double macro_stress[nvoi] = { 0.0 };
			for (int i = 0; i < nvoi; ++i)
				for (int j = 0; j < nvoi; ++j)
					macro_stress[i] += ctan_lin[i * nvoi + j] * macro_strain[j];

			double inv = macro_stress[0];
			for(int i = 1; i < dim; ++i)
				inv += macro_stress[i];

			if (fabs(inv) > _inv_max)
				inv_max = fabs(inv);

			return (fabs(inv) < _inv_tol);
		}

		void print_stress() const
		{
			cout << "\t\t";
			for(int i = 0; i < nvoi; ++i)
				cout << setw(14) << macro_stress[i] << "\t";
			cout << endl;
		}

		void print_strain() const
		{
			cout << "\t\t";
			for(int i = 0; i < nvoi; ++i)
				cout << setw(14) << macro_strain[i] << "\t";
			cout << endl;
		}

		void print_ctan() const
		{
			for(int i = 0; i < nvoi; ++i) {
				cout << "\t\t";
				for(int j = 0; j < nvoi; ++j)
					cout << setw(14) << macro_ctan[i * nvoi + j] << "\t";
				cout << endl;
			}
		}

		void print()
		{
			printf("gp: %p u_k: %p v_k: %p\n", this, u_k, int_vars_k);
		}
};

#endif
