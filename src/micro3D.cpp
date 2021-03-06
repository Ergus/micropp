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


template<>
micropp<3>::micropp(const int _ngp, const int size[3], const int _micro_type,
                    const double *_micro_params, const material_t *_materials):
	ngp(_ngp),
	nx(size[0]), ny(size[1]), nz(size[2]),
	nn(nx * ny * nz),
	nex(nx - 1), ney(ny - 1), nez(nz - 1),
	nelem(nex * ney * nez),
	lx(_micro_params[0]), ly(_micro_params[1]),	lz(_micro_params[2]),
	dx(lx / nex), dy(ly / ney),	dz(lz / nez),

	width(_micro_params[3]), inv_tol(_micro_params[4]),
	micro_type(_micro_type), num_int_vars(nelem * 8 * NUM_VAR_GP)
{
	initialize(_micro_params, _materials);
}


template<>
int micropp<3>::get_elem_type(int ex, int ey, int ez)
{
	assert(micro_type == 0 || micro_type == 1);

	if (micro_type == 0) { // sphere in the center

		const double coor[3] = {
			ex * dx + dx / 2,
			ey * dy + dy / 2,
			ez * dz + dz / 2 };
		const double center[3] = { lx / 2, ly / 2, lz / 2 };
		const double rad = micro_params[3];
		double tmp = 0.;
		for (int i = 0; i < 3; ++i)
			tmp += (center[i] - coor[i]) * (center[i] - coor[i]);

		return (tmp < rad * rad);

	} else if (micro_type == 1) { // 2 flat layers in y dir

		const double y = ey * dy + dy / 2;
		return (y < width);
	}

	cerr << "Invalid micro_type = " << micro_type << endl;
	return -1;
}


template <>
void micropp<3>::set_displ_bc(const double *eps)
{
	const double eps_t[3][3] = {
		{ eps[0], 0.5 * eps[3], 0.5 * eps[4] },
		{ 0.5 * eps[3], eps[1], 0.5 * eps[5] },
		{ 0.5 * eps[4], 0.5 * eps[5], eps[2] } };


	for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			const int n = nod_index3D(i, j, 0); // z = 0
			const double coor[3] = { i * dx, j * dy, 0 };
			mvp_3(eps_t, coor, &u[n * dim]);
		}
	}

	for (int i = 0; i < nx; ++i) {
		for (int j = 0; j < ny; ++j) {
			const int n = nod_index3D(i, j, nz - 1); // z = lz
			const double coor[3] = { i * dx, j * dy, lz };
			mvp_3(eps_t, coor, &u[n * dim]);
		}
	}

	for (int i = 0; i < nx; ++i) {
		for (int k = 1; k < nz - 1; ++k) {
			const int n = nod_index3D(i, 0, k); // y = 0
			const double coor[3] = { i * dx, 0, k * dz };
			mvp_3(eps_t, coor, &u[n * dim]);
		}
	}

	for (int i = 0; i < nx; ++i) {
		for (int k = 1; k < nz - 1; ++k) {
			const int n = nod_index3D(i, ny - 1, k); // y = ly
			const double coor[3] = { i * dx, ly , k * dz };
			mvp_3(eps_t, coor, &u[n * dim]);
		}
	}

	for (int j = 1; j < ny - 1; ++j) {
		for (int k = 1; k < nz - 1; ++k) {
			const int n = nod_index3D(0, j, k); // x = 0
			const double coor[3] = { 0, j * dy , k * dz };
			mvp_3(eps_t, coor, &u[n * dim]);
		}
	}

	for (int j = 1; j < ny - 1; ++j) {
		for (int k = 1; k < nz - 1; ++k) {
			const int n = nod_index3D(nx - 1, j, k); // x = lx
			const double coor[3] = { ly, j * dy , k * dz };
			mvp_3(eps_t, coor, &u[n * dim]);
		}
	}
}


template <>
void micropp<3>::calc_bmat(int gp, double bmat[nvoi][npe * dim]) const
{
	INST_START;

	const double dsh[8][3] = {
		{ -(1 - xg[gp][1]) * (1 - xg[gp][2]) / 8. * 2. / dx,
		  -(1 - xg[gp][0]) * (1 - xg[gp][2]) / 8. * 2. / dy,
		  -(1 - xg[gp][0]) * (1 - xg[gp][1]) / 8. * 2. / dz },
		{ +(1 - xg[gp][1]) * (1 - xg[gp][2]) / 8. * 2. / dx,
		  -(1 + xg[gp][0]) * (1 - xg[gp][2]) / 8. * 2. / dy,
		  -(1 + xg[gp][0]) * (1 - xg[gp][1]) / 8. * 2. / dz },
		{ +(1 + xg[gp][1]) * (1 - xg[gp][2]) / 8. * 2. / dx,
		  +(1 + xg[gp][0]) * (1 - xg[gp][2]) / 8. * 2. / dy,
		  -(1 + xg[gp][0]) * (1 + xg[gp][1]) / 8. * 2. / dz },
		{ -(1 + xg[gp][1]) * (1 - xg[gp][2]) / 8. * 2. / dx,
		  +(1 - xg[gp][0]) * (1 - xg[gp][2]) / 8. * 2. / dy,
		  -(1 - xg[gp][0]) * (1 + xg[gp][1]) / 8. * 2. / dz },
		{ -(1 - xg[gp][1]) * (1 + xg[gp][2]) / 8. * 2. / dx,
		  -(1 - xg[gp][0]) * (1 + xg[gp][2]) / 8. * 2. / dy,
		  +(1 - xg[gp][0]) * (1 - xg[gp][1]) / 8. * 2. / dz },
		{ +(1 - xg[gp][1]) * (1 + xg[gp][2]) / 8. * 2. / dx,
		  -(1 + xg[gp][0]) * (1 + xg[gp][2]) / 8. * 2. / dy,
		  +(1 + xg[gp][0]) * (1 - xg[gp][1]) / 8. * 2. / dz },
		{ +(1 + xg[gp][1]) * (1 + xg[gp][2]) / 8. * 2. / dx,
		  +(1 + xg[gp][0]) * (1 + xg[gp][2]) / 8. * 2. / dy,
		  +(1 + xg[gp][0]) * (1 + xg[gp][1]) / 8. * 2. / dz },
		{ -(1 + xg[gp][1]) * (1 + xg[gp][2]) / 8. * 2. / dx,
		  +(1 - xg[gp][0]) * (1 + xg[gp][2]) / 8. * 2. / dy,
		  +(1 - xg[gp][0]) * (1 + xg[gp][1]) / 8. * 2. / dz } };

	for (int i = 0; i < 8; ++i) {
		bmat[0][i * dim    ] = dsh[i][0];
		bmat[0][i * dim + 1] = 0;
		bmat[0][i * dim + 2] = 0;
		bmat[1][i * dim    ] = 0;
		bmat[1][i * dim + 1] = dsh[i][1];
		bmat[1][i * dim + 2] = 0;
		bmat[2][i * dim    ] = 0;
		bmat[2][i * dim + 1] = 0;
		bmat[2][i * dim + 2] = dsh[i][2];
		bmat[3][i * dim    ] = dsh[i][1];
		bmat[3][i * dim + 1] = dsh[i][0];
		bmat[3][i * dim + 2] = 0;
		bmat[4][i * dim    ] = dsh[i][2];
		bmat[4][i * dim + 1] = 0;
		bmat[4][i * dim + 2] = dsh[i][0];
		bmat[5][i * dim    ] = 0;
		bmat[5][i * dim + 1] = dsh[i][2];
		bmat[5][i * dim + 2] = dsh[i][1];
	}
}


template <>
template <>
void micropp<3>::get_elem_rhs(double be[npe * dim],
                              int ex, int ey, int ez) const
{
	constexpr int npedim = npe * dim;
	double bmat[nvoi][npe * dim], stress_gp[nvoi], strain_gp[nvoi];
	const double wg = dx * dy * dz / npe;

	memset(be, 0, npedim * sizeof(double));

	for (int gp = 0; gp < npe; ++gp) {

		calc_bmat(gp, bmat);

		get_strain(gp, strain_gp, ex, ey, ez);
		get_stress(gp, strain_gp, stress_gp, ex, ey, ez);

		for (int i = 0; i < npe * dim; i++)
			for (int j = 0; j < nvoi; j++)
				be[i] += bmat[j][i] * stress_gp[j] * wg;
	}
}


template<>
double micropp<3>::assembly_rhs()
{
	INST_START;

	memset(b, 0.0, nn * dim * sizeof(double));

	double be[3 * 8];
	int index[3 * 8];

	for (int ex = 0; ex < nex; ex++) {
		for (int ey = 0; ey < ney; ey++) {
			for (int ez = 0; ez < nez; ez++) {

				int n[8];
				get_elem_nodes(n, ex, ey, ez);

				for (int j = 0; j < npe; ++j)
					for (int d = 0; d < dim; ++d)
						index[j * dim + d] = n[j] * dim + d;

				get_elem_rhs(be, ex, ey, ez);

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
	for (int i = 0; i < nn * dim; ++i)
		b[i] = -b[i];

	double norm = 0.0;
	for (int i = 0; i < nn * dim; ++i)
		norm += b[i] * b[i];
	norm = sqrt(norm);

	return norm;
}


template <>
template <>
void micropp<3>::get_elem_mat(double Ae[npe * dim * npe * dim],
		int ex, int ey, int ez) const
{
	INST_START;
	const int e = glo_elem(ex, ey, ez);
	const material_t material = get_material(e);

	double ctan[nvoi][nvoi];
	constexpr int npedim = npe * dim;
	constexpr int npedim2 = npedim * npedim;
	const double wg = (1. / 8.) * dx * dy * dz;

	double TAe[npedim2] = { 0.0 };

	// memset(Ae, 0, npedim2 * sizeof(double));

	for (int gp = 0; gp < npe; ++gp) {

		double eps[6];
		get_strain(gp, eps, ex, ey, ez);
		const double *eps_p_old = &vars_old[intvar_ix(e, gp, 0)];
		double alpha_old = vars_old[intvar_ix(e, gp, 6)];

		if (material.plasticity)
			plastic_get_ctan(&material, eps, eps_p_old, alpha_old, ctan);
		else
			isolin_get_ctan(&material, ctan);

		double bmat[nvoi][npedim], cxb[nvoi][npedim];
		calc_bmat(gp, bmat);

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


template <>
void micropp<3>::assembly_mat()
{
	INST_START;

	ell_set_zero_mat(&A);

	double Ae[npe * dim * npe * dim];
	for (int ex = 0; ex < nex; ++ex) {
		for (int ey = 0; ey < ney; ++ey) {
			for (int ez = 0; ez < nez; ++ez) {
				get_elem_mat(Ae, ex, ey, ez);
				ell_add_3D(&A, ex, ey, ez, Ae);
			}
		}
	}
	ell_set_bc_3D(&A);
}


template <>
void micropp<3>::calc_ave_stress(double stress_ave[6]) const
{
	const double wg = (1. / 8.) * dx * dy * dz;

	memset(stress_ave, 0, nvoi * sizeof(double));

	for (int ex = 0; ex < nex; ++ex) {
		for (int ey = 0; ey < ney; ++ey) {
			for (int ez = 0; ez < nez; ++ez) {

				double stress_aux[6] = { 0.0 };

				for (int gp = 0; gp < npe; ++gp) {

					double stress_gp[6];
					double strain_gp[6];

					get_strain(gp, strain_gp, ex, ey, ez);
					get_stress(gp, strain_gp, stress_gp, ex, ey, ez);
					for (int v = 0; v < nvoi; ++v)
						stress_aux[v] += stress_gp[v] * wg;

				}
				for (int v = 0; v < nvoi; ++v)
					stress_ave[v] += stress_aux[v];
			}
		}
	}

	for (int v = 0; v < nvoi; ++v)
		stress_ave[v] /= (lx * ly);
}


template <>
void micropp<3>::calc_ave_strain(double strain_ave[6]) const
{
	const double wg = (1. / 8.) * dx * dy * dz;
	memset(strain_ave, 0, nvoi * sizeof(double));

	for (int ex = 0; ex < nex; ++ex) {
		for (int ey = 0; ey < ney; ++ey) {
			for (int ez = 0; ez < nez; ++ez) {

				double strain_aux[6] = { 0.0 };

				for (int gp = 0; gp < npe; ++gp) {
					double strain_gp[6];

					get_strain(gp, strain_gp, ex, ey, ez);
					for (int v = 0; v < nvoi; ++v)
						strain_aux[v] += strain_gp[v] * wg;
				}

				for (int v = 0; v < nvoi; v++)
					strain_ave[v] += strain_aux[v];
			}
		}
	}

	for (int v = 0; v < nvoi; v++)
		strain_ave[v] /= (lx * ly);
}


template<>
void micropp<3>::calc_fields()
{
	const double ivol = 1. / (dx * dy * dz);
	const double wg = (1. / 8.) * dx * dy * dz;

	for (int ex = 0; ex < nex; ++ex) {
		for (int ey = 0; ey < ney; ++ey) {
			for (int ez = 0; ez < nez; ++ez) {

				double strain_aux[6] = { 0.0 };
				double stress_aux[6] = { 0.0 };

				for (int gp = 0; gp < npe; ++gp) {

					double stress_gp[6], strain_gp[6];

					get_strain(gp, strain_gp, ex, ey, ez);
					get_stress(gp, strain_gp, stress_gp, ex, ey, ez);
					for (int v = 0; v < nvoi; v++) {
						strain_aux[v] += strain_gp[v] * wg;
						stress_aux[v] += stress_gp[v] * wg;
					}
				}

				const int e = glo_elem(ex, ey, ez);
				for (int v = 0; v < nvoi; v++) {
					elem_strain[e * nvoi + v] = strain_aux[v] * ivol;
					elem_stress[e * nvoi + v] = stress_aux[v] * ivol;
				}
			}
		}
	}
}


template<>
bool micropp<3>::calc_vars_new()
{
	INST_START;

    bool nl_flag = false;

	for (int ex = 0; ex < nex; ex++) {
		for (int ey = 0; ey < ney; ey++) {
			for (int ez = 0; ez < nez; ez++) {

				for (int gp = 0; gp < 8; gp++) {

					int e = glo_elem(ex, ey, ez);
					const material_t material = get_material(e);

					const double *eps_p_old = &vars_old[intvar_ix(e, gp, 0)];
					double alpha_old = vars_old[intvar_ix(e, gp, 6)];
					double *eps_p_new = &vars_new[intvar_ix(e, gp, 0)];
					double *alpha_new = &vars_new[intvar_ix(e, gp, 6)];
					double eps[6];
					get_strain(gp, eps, ex, ey, ez);

					nl_flag |= plastic_evolute(
							&material, eps, eps_p_old, alpha_old, 
							eps_p_new, alpha_new);
				}
			}
		}
	}

	return nl_flag;
}

// Explicit instantiation
template class micropp<3>;
