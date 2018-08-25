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

void micropp::set_displ_bc_3D(const double *eps, double *u) const
{
	const double eps_t[] = {       eps[0], 0.5 * eps[3], 0.5 * eps[4] ,
	                         0.5 * eps[3],       eps[1], 0.5 * eps[5] ,
	                         0.5 * eps[4], 0.5 * eps[5],       eps[2] };


	for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			const int n = nod_index3D(i, j, 0); // z = 0
			const double coor[3] = { i * dx, j * dy, 0 };
			mvp(dim, eps_t, coor, &u[n * dim]);
		}
	}

	for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			const int n = nod_index3D(i, j, nz - 1); // z = lz
			const double coor[3] = { i * dx, j * dy, lz };
			mvp(dim, eps_t, coor, &u[n * dim]);
		}
	}

	for (int i = 0; i < nx; ++i) {
		for (int k = 1; k < nz - 1; ++k) {
			const int n = nod_index3D(i, 0, k); // y = 0
			const double coor[3] = { i * dx, 0, k * dz };
			mvp(dim, eps_t, coor, &u[n * dim]);
		}
	}

	for (int i = 0; i < nx; ++i) {
		for (int k = 1; k < nz - 1; ++k) {
			const int n = nod_index3D(i, ny - 1, k); // y = ly
			const double coor[3] = { i * dx, ly , k * dz };
			mvp(dim, eps_t, coor, &u[n * dim]);
		}
	}

	for (int j = 1; j < ny - 1; ++j) {
		for (int k = 1; k < nz - 1; ++k) {
			const int n = nod_index3D(0, j, k); // x = 0
			const double coor[3] = { 0, j * dy , k * dz };
			mvp(dim, eps_t, coor, &u[n * dim]);
		}
	}

	for (int j = 1; j < ny - 1; ++j) {
		for (int k = 1; k < nz - 1; ++k) {
			const int n = nod_index3D(nx - 1, j, k); // x = lx
			const double coor[3] = { ly, j * dy , k * dz };
			mvp(dim, eps_t, coor, &u[n * dim]);
		}
	}
}


void micropp::calc_bmat_3D(int gp, double *bmat) const
{
	const int npedim = npe * dim;

	constexpr double xg[] = { -CONSTXG, -CONSTXG, -CONSTXG ,
	                          +CONSTXG, -CONSTXG, -CONSTXG ,
	                          +CONSTXG, +CONSTXG, -CONSTXG ,
	                          -CONSTXG, +CONSTXG, -CONSTXG ,
	                          -CONSTXG, -CONSTXG, +CONSTXG ,
	                          +CONSTXG, -CONSTXG, +CONSTXG ,
	                          +CONSTXG, +CONSTXG, +CONSTXG ,
	                          -CONSTXG, +CONSTXG, +CONSTXG  };

	const double dsh[] = {
		  -(1 - xg[gp * dim + 1]) * (1 - xg[gp * dim + 2]) / 8. * 2. / dx,
		  -(1 - xg[gp * dim + 0]) * (1 - xg[gp * dim + 2]) / 8. * 2. / dy,
		  -(1 - xg[gp * dim + 0]) * (1 - xg[gp * dim + 1]) / 8. * 2. / dz,
		  +(1 - xg[gp * dim + 1]) * (1 - xg[gp * dim + 2]) / 8. * 2. / dx,
		  -(1 + xg[gp * dim + 0]) * (1 - xg[gp * dim + 2]) / 8. * 2. / dy,
		  -(1 + xg[gp * dim + 0]) * (1 - xg[gp * dim + 1]) / 8. * 2. / dz,
		  +(1 + xg[gp * dim + 1]) * (1 - xg[gp * dim + 2]) / 8. * 2. / dx,
		  +(1 + xg[gp * dim + 0]) * (1 - xg[gp * dim + 2]) / 8. * 2. / dy,
		  -(1 + xg[gp * dim + 0]) * (1 + xg[gp * dim + 1]) / 8. * 2. / dz,
		  -(1 + xg[gp * dim + 1]) * (1 - xg[gp * dim + 2]) / 8. * 2. / dx,
		  +(1 - xg[gp * dim + 0]) * (1 - xg[gp * dim + 2]) / 8. * 2. / dy,
		  -(1 - xg[gp * dim + 0]) * (1 + xg[gp * dim + 1]) / 8. * 2. / dz,
		  -(1 - xg[gp * dim + 1]) * (1 + xg[gp * dim + 2]) / 8. * 2. / dx,
		  -(1 - xg[gp * dim + 0]) * (1 + xg[gp * dim + 2]) / 8. * 2. / dy,
		  +(1 - xg[gp * dim + 0]) * (1 - xg[gp * dim + 1]) / 8. * 2. / dz,
		  +(1 - xg[gp * dim + 1]) * (1 + xg[gp * dim + 2]) / 8. * 2. / dx,
		  -(1 + xg[gp * dim + 0]) * (1 + xg[gp * dim + 2]) / 8. * 2. / dy,
		  +(1 + xg[gp * dim + 0]) * (1 - xg[gp * dim + 1]) / 8. * 2. / dz,
		  +(1 + xg[gp * dim + 1]) * (1 + xg[gp * dim + 2]) / 8. * 2. / dx,
		  +(1 + xg[gp * dim + 0]) * (1 + xg[gp * dim + 2]) / 8. * 2. / dy,
		  +(1 + xg[gp * dim + 0]) * (1 + xg[gp * dim + 1]) / 8. * 2. / dz,
		  -(1 + xg[gp * dim + 1]) * (1 + xg[gp * dim + 2]) / 8. * 2. / dx,
		  +(1 - xg[gp * dim + 0]) * (1 + xg[gp * dim + 2]) / 8. * 2. / dy,
		  +(1 - xg[gp * dim + 0]) * (1 + xg[gp * dim + 1]) / 8. * 2. / dz };

	for (int i = 0; i < npe; ++i) {
		bmat[0 * npedim + i * dim    ] = dsh[i * dim + 0];
		bmat[0 * npedim + i * dim + 1] = 0;
		bmat[0 * npedim + i * dim + 2] = 0;
		bmat[1 * npedim + i * dim    ] = 0;
		bmat[1 * npedim + i * dim + 1] = dsh[i * dim + 1];
		bmat[1 * npedim + i * dim + 2] = 0;
		bmat[2 * npedim + i * dim    ] = 0;
		bmat[2 * npedim + i * dim + 1] = 0;
		bmat[2 * npedim + i * dim + 2] = dsh[i * dim + 2];
		bmat[3 * npedim + i * dim    ] = dsh[i * dim + 1];
		bmat[3 * npedim + i * dim + 1] = dsh[i * dim + 0];
		bmat[3 * npedim + i * dim + 2] = 0;
		bmat[4 * npedim + i * dim    ] = dsh[i * dim + 2];
		bmat[4 * npedim + i * dim + 1] = 0;
		bmat[4 * npedim + i * dim + 2] = dsh[i * dim + 0];
		bmat[5 * npedim + i * dim    ] = 0;
		bmat[5 * npedim + i * dim + 1] = dsh[i * dim + 2];
		bmat[5 * npedim + i * dim + 2] = dsh[i * dim + 1];
	}
}


double micropp::assembly_rhs_3D(const double *u, const double *old,
                                double *b) const
{
	memset(b, 0, nndim * sizeof(double));

	double be[dim * npe];
	int index[dim * npe];

	for (int ez = 0; ez < nez; ++ez) {
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex) {

				int n[npe];
				get_elem_nodes(n, ex, ey, ez); // final call, no array access

				for (int j = 0; j < npe; ++j)
					for (int d = 0; d < dim; ++d)
						index[j * dim + d] = n[j] * dim + d;

				get_elem_rhs(u, old, be, ex, ey, ez);

				for (int i = 0; i < npe * dim; ++i)
					b[index[i]] += be[i];
			}
		}
	}

	// boundary conditions
	for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			const int n = nod_index3D(i, j, 0); // z = 0
			memset(&b[n * dim], 0, dim * sizeof(double));
		}
	}

	for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			const int n = nod_index3D(i, j, nz - 1); // z = lx
			memset(&b[n * dim], 0, dim * sizeof(double));
		}
	}

	for (int i = 0; i < nx; ++i) {
		for (int k = 1; k < nz - 1; ++k) {
			const int n = nod_index3D(i, 0, k); // y = 0
			memset(&b[n * dim], 0, dim * sizeof(double));
		}
	}

	for (int i = 0; i < nx; ++i) {
		for (int k = 1; k < nz - 1; ++k) {
			const int n = nod_index3D(i, ny - 1, k); // y = ly
			memset(&b[n * dim], 0, dim * sizeof(double));
		}
	}

	for (int j = 1; j < ny - 1; ++j) {
		for (int k = 1; k < nz - 1; ++k) {
			const int n = nod_index3D(0, j, k); // x = 0
			memset(&b[n * dim], 0, dim * sizeof(double));
		}
	}

	for (int j = 1; j < ny - 1; j++) {
		for (int k = 1; k < nz - 1; ++k) {
			const int n = nod_index3D(nx - 1, j, k); // x = lx
			memset(&b[n * dim], 0, dim * sizeof(double));
		}
	}

	// Common
	for (int i = 0; i < nndim; ++i)
		b[i] = -b[i];

	double norm = 0.0;
	for (int i = 0; i < nndim; ++i)
		norm += b[i] * b[i];
	norm = sqrt(norm);

	return norm;
}


void micropp::get_elem_mat_3D(const double *u, const double *old,
                              double *Ae, int ex, int ey, int ez) const
{
	const int e = glo_elem(ex, ey, ez);
	const material_t material = get_material(e);

	double ctan[6][6];
	const int npedim = npe * dim;
	const int npedim2 = npedim * npedim;

	double *TAe = (double *) alloca(npedim2 * sizeof(double));
	memset(TAe, 0, npedim2 * sizeof(double));

	// memset(Ae, 0, npedim2 * sizeof(double));

	for (int gp = 0; gp < npe; ++gp) {

		double eps[6];
		get_strain(u, gp, eps, ex, ey, ez);
		const double *eps_p_old = &old[intvar_ix(e, gp, 0)];
		double alpha_old = old[intvar_ix(e, gp, 6)];

		if (material.plasticity)
			plastic_get_ctan_3D(&material, eps, eps_p_old, alpha_old, ctan);
		else
			isolin_get_ctan_3D(&material, ctan);

		double bmat[nvoi][npedim], cxb[nvoi][npedim];
		calc_bmat_3D(gp, (double *) bmat);

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


void micropp::assembly_mat_3D(const double *u, const double *old,
                              ell_matrix *A) const
{
	ell_set_zero_mat(A);

	double Ae[npe * dim * npe * dim];
	for (int ex = 0; ex < nex; ++ex) {
		for (int ey = 0; ey < ney; ++ey) {
			for (int ez = 0; ez < nez; ++ez) {
				get_elem_mat_3D(u, old, Ae, ex, ey, ez);
				ell_add_3D(A, ex, ey, ez, Ae);
			}
		}
	}
	ell_set_bc_3D(A);
}


void micropp::plastic_get_stress_3D(const material_t *material,
                                    const double eps[6],
                                    const double eps_p_old[6],
                                    double alpha_old,
                                    double stress[6]) const
{
	double dl, normal[6], s_trial[6];
	bool nl_flag = plastic_law_if(material, eps, eps_p_old, alpha_old,
	                              &dl, normal, s_trial);

	//sig_2 = s_trial + K * tr(eps) * 1 - 2 * mu * dl * normal;
	memcpy(stress, s_trial, 6 * sizeof(double));

	for (int i = 0; i < 3; ++i)
		stress[i] += material->k * (eps[0] + eps[1] + eps[2]);

	for (int i = 0; i < 6; ++i)
		stress[i] -= 2 * material->mu * dl * normal[i];
}

void micropp::plastic_get_ctan_3D(const material_t *material,
								  const double eps[6],
								  const double eps_p_old[6],
								  double alpha_old,
								  double ctan[6][6]) const
{
	double stress_0[6];
	plastic_get_stress_3D(material, eps, eps_p_old, alpha_old, stress_0);

	for (int i = 0; i < nvoi; ++i) {

		double eps_1[6];
		memcpy(eps_1, eps, nvoi * sizeof(double));
		eps_1[i] += D_EPS_CTAN;

		double stress_1[6];
		plastic_get_stress_3D(material, eps_1, eps_p_old, alpha_old, stress_1);

		for (int j = 0; j < nvoi; ++j)
			ctan[j][i] = (stress_1[j] - stress_0[j]) / D_EPS_CTAN;
	}
}


void micropp::get_dev_tensor_3D(const double tensor[6],
                                double tensor_dev[6]) const
{
	memcpy(tensor_dev, tensor, nvoi * sizeof(double));
	for (int i = 0; i < 3; i++)
		tensor_dev[i] -= (1 / 3.0) * (tensor[0] + tensor[1] + tensor[2]);
}


void micropp::isolin_get_ctan_3D(const material_t *material,
								 double ctan[6][6]) const
{
	// C = lambda * (1x1) + 2 mu I
	memset(ctan, 0, nvoi * nvoi * sizeof(double));

	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			ctan[i][j] += material->lambda;

	for (int i = 0; i < 3; ++i)
		ctan[i][i] += 2 * material->mu;

	for (int i = 3; i < 6; ++i)
		ctan[i][i] = material->mu;
}
