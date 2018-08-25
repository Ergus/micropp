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

void micropp::set_displ_bc_2D(const double *eps, double *u) const
{
	const double eps_t[] = {       eps[0], 0.5 * eps[2],
	                         0.5 * eps[2],       eps[1] };


	for (int i = 0; i < nx; ++i) {
		const int n = nod_index(i, 0, 0); // y = 0
		const double coor[2] = { i * dx, 0 };
		mvp(dim, eps_t, coor, &u[n * dim]);
	}

	for (int i = 0; i < nx; ++i) {
		const int n = nod_index(i, ny - 1, 0); // y = ly
		const double coor[2] = { i * dx, ly };
		mvp(dim, eps_t, coor, &u[n * dim]);
	}

	for (int j = 1; j < ny - 1; ++j) {
		const int n = nod_index(0, j, 0); // x = 0
		const double coor[2] = { 0, j * dy };
		mvp(dim, eps_t, coor, &u[n * dim]);
	}

	for (int j = 1; j < ny - 1; ++j) {
		const int n = nod_index(nx - 1, j, 0); // x = lx
		const double coor[2] = { lx, j * dy };
		mvp(dim, eps_t, coor, &u[n * dim]);
	}
}


void micropp::calc_bmat_2D(int gp, double *bmat)	const
{
	const int npedim = npe * dim;

	constexpr double xg[] = { -CONSTXG, -CONSTXG,
	                          +CONSTXG, -CONSTXG,
	                          +CONSTXG, +CONSTXG,
	                          -CONSTXG, +CONSTXG };

	const double dsh[] = {
		-(1 - xg[gp * dim + 1]) / 4 * 2 / dx, -(1 - xg[gp * dim + 0]) / 4 * 2 / dy,
		+(1 - xg[gp * dim + 1]) / 4 * 2 / dx, -(1 + xg[gp * dim + 0]) / 4 * 2 / dy,
		+(1 + xg[gp * dim + 1]) / 4 * 2 / dx, +(1 + xg[gp * dim + 0]) / 4 * 2 / dy,
		-(1 + xg[gp * dim + 1]) / 4 * 2 / dx, +(1 - xg[gp * dim + 0]) / 4 * 2 / dy };

	for (int i = 0; i < npe; ++i) {
		bmat[0 * npedim + i * dim    ] = dsh[i * dim + 0];
		bmat[0 * npedim + i * dim + 1] = 0;
		bmat[1 * npedim + i * dim    ] = 0;
		bmat[1 * npedim + i * dim + 1] = dsh[i * dim + 1];
		bmat[2 * npedim + i * dim    ] = dsh[i * dim + 1];
		bmat[2 * npedim + i * dim + 1] = dsh[i * dim + 0];
	}
}


double micropp::assembly_rhs_2D(const double *u, const double *old,
                                double *b) const
{
	memset(b, 0, nndim * sizeof(double));

	double be[dim * npe];
	int index[dim * npe];

	for (int ex = 0; ex < nex; ++ex) {
		for (int ey = 0; ey < ney; ++ey) {

			int n[npe];
			get_elem_nodes(n, ex, ey);

			for (int j = 0; j < npe; ++j)
				for (int d = 0; d < dim; ++d)
					index[j * dim + d] = n[j] * dim + d;

			get_elem_rhs(u, old, be, ex, ey);

			for (int i = 0; i < npe * dim; ++i)
				b[index[i]] += be[i];
		}
	}

	// boundary conditions
	for (int i = 0; i < nx; i++) {
		const int n = nod_index2D(i, 0); // y = 0
		memset(&b[n * dim], 0, dim * sizeof(double));
	}
	for (int i = 0; i < nx; i++) {
		const int n = nod_index2D(i, ny - 1); // y = ly
		memset(&b[n * dim], 0, dim * sizeof(double));
	}
	for (int j = 1; j < ny - 1; j++) {
		const int n = nod_index2D(0, j); // x = 0
		memset(&b[n * dim], 0, dim * sizeof(double));
	}
	for (int j = 1; j < ny - 1; j++) {
		const int n = nod_index2D(nx - 1, j); // x = lx
		memset(&b[n * dim], 0, dim * sizeof(double));
	}

	// Common part
	for (int i = 0; i < nndim; ++i)
		b[i] = -b[i];

	double norm = 0.0;
	for (int i = 0; i < nndim; ++i)
		norm += b[i] * b[i];
	norm = sqrt(norm);

	return norm;
}


void micropp::get_elem_mat_2D(const double *u, const double *old,
                           double *Ae, int ex, int ey) const
{
	const int e = glo_elem(ex, ey, 0);
	const material_t material = get_material(e);

	const int npedim = npe * dim;
	const int npedim2 = npedim * npedim;

	double *TAe = (double *) alloca(npedim2 * sizeof(double));
	memset(TAe, 0, npedim2 * sizeof(double));

	const double E = material.E;
	const double nu = material.nu;
	const bool plasticity = material.plasticity;

	double ctan[3][3] = {{(1 - nu),       nu,                0 },
	                     {      nu, (1 - nu),                0 },
	                     {       0,        0, (1 - 2 * nu) / 2 }};

	for (int i = 0; i < nvoi; i++)
		for (int j = 0; j < nvoi; j++)
			ctan[i][j] *= E / ((1 + nu) * (1 - 2 * nu));

	for (int gp = 0; gp < npe; ++gp) {

		double bmat[nvoi][npedim], cxb[nvoi][npedim];

		calc_bmat_2D(gp, (double *) bmat);

		for (int i = 0; i < nvoi; ++i) {
			for (int j = 0; j < npedim; ++j) {
				double tmp = 0.0;
				for (int k = 0; k < nvoi; ++k)
					tmp += ctan[i][k] * bmat[k][j];
				cxb[i][j] = tmp;
			}
		}

		for (int m = 0; m < nvoi; ++m) {
			for (int i = 0; i < npedim; ++i) {
				const int inpedim = i * npedim;
				const double bmatmi = bmat[m][i];
				for (int j = 0; j < npedim; ++j)
					TAe[inpedim + j] += bmatmi * cxb[m][j] * wg;
			}
		}
		memcpy(Ae, TAe, npedim2 * sizeof(double));
	}
}


void micropp::assembly_mat_2D(const double *u, const double *old,
                              ell_matrix *A) const
{
	ell_set_zero_mat(A);

	double Ae[npe * dim * npe * dim];
	for (int ex = 0; ex < nex; ++ex) {
		for (int ey = 0; ey < ney; ++ey) {
			get_elem_mat_2D(u, old, Ae, ex, ey);
			ell_add_2D(A, ex, ey, Ae);
		}
	}
	ell_set_bc_2D(A);
}

