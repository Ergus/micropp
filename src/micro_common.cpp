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

template<int tdim>
micropp<tdim>::micropp(const int _ngp, const int size[3],
                       const int _micro_type,
                       const double _micro_params[5],
                       const material_t *_materials):
	ngp(_ngp),
	nx(size[0]), ny(size[1]),
	nz((tdim == 3) ? size[2] : 1),

	nn(nx * ny * nz),
	nndim(nn * dim),

	nex(nx - 1), ney(ny - 1),
	nez((tdim == 3) ? (nz - 1) : 1),

	nelem(nex * ney * nez),
	lx(_micro_params[0]), ly(_micro_params[1]),
	lz((tdim == 3) ? _micro_params[2] : 0.0),
	dx(lx / nex), dy(ly / ney),
	dz((tdim == 3) ? lz / nez : 0.0),

	special_param(_micro_params[3]), inv_tol(_micro_params[4]),

	wg(((tdim == 3) ? dx * dy * dz : dx * dy) / npe),
	vol_tot((tdim == 3) ? lx * ly * lz : lx * ly),
	ivol(1.0 / (wg * npe)),
	micro_type(_micro_type), num_int_vars(nelem * 8 * NUM_VAR_GP), 
	ell_cols_size(mypow(3, dim) * dim * nn * dim)
{
	INST_CONSTRUCT; // Initialize the Intrumentation
	// GP list here
	gp_list = new gp_t<tdim>[ngp]();

	for (int gp = 0; gp < ngp; gp++) {
		gp_list[gp].u_n = (double *) calloc(nndim, sizeof(double));
		gp_list[gp].u_k = (double *) malloc(nndim * sizeof(double));
	}

	// Material list Here
	int nParams;
	if (micro_type == 0) {
		// mat 1 = matrix
		// mat 2 = sphere
		numMaterials = 2;
		nParams = 5;
	} else if (micro_type == 1) {
		// mat 1 = layer 1
		// mat 2 = layer 2
		numMaterials = 2;
		nParams = 5;
	}

	material_list = (material_t *) rrd_malloc(numMaterials * sizeof(material_t));

	{
		ofstream file;
		file.open("micropp_materials.dat");
		file << scientific;
		for (int i = 0; i < numMaterials; ++i) {
			set(_materials[i], &material_list[i]);

			file << setw(14)
			     << _materials[i].E << " " << material_list[i].nu << " "
			     << _materials[i].Sy << " " << material_list[i].Ka << endl;
		}
	}

	elem_type = (int *) rrd_malloc(nelem * sizeof(int));
	for (int ez = 0; ez < nez; ++ez) {
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex) {
				const int e_i = glo_elem(ex, ey, ez);
				int type = get_elem_type(ex, ey, ez);
				set(type, &elem_type[e_i]);
			}
		}
	}

	elem_stress = (double *) rrd_malloc(nelem * nvoi * sizeof(double));
	elem_strain = (double *) rrd_malloc(nelem * nvoi * sizeof(double));
	// vars_old_aux = (double *) calloc(num_int_vars, sizeof(double));
	// vars_new_aux = (double *) malloc(num_int_vars * sizeof(double));

	assert(elem_stress && elem_strain && elem_type);

	output_files_header = false;

	const int ns[3] = { nx, ny, nz };
	// This needs to be released manually!!
	ell_cols = ell_init_cols(dim, dim, ns);


	calc_ctan_lin();

	ofstream file;
	file.close();

	file.open("micropp_convergence.dat");
	file.close();
	file.open("micropp_eps_sig_ctan.dat");
	file.close();
	file.open("micropp_int_vars_n.dat");
	file.close();
}

template <int tdim>
micropp<tdim>::~micropp()
{
	INST_DESTRUCT;

	free(elem_stress);
	free(elem_strain);
	free(elem_type);

	delete [] gp_list;
}


template <int tdim>
int micropp<tdim>::get_nl_flag(int gp_id) const
{
	assert(gp_id < ngp);
	assert(gp_id >= 0);
	return gp_list[gp_id].allocated;
}


template <int tdim>
bool micropp<tdim>::is_linear(const double *macro_strain, double *inv_max) const
{
	double macro_stress[6] = { 0.0 };
	for (int i = 0; i < nvoi; ++i)
		for (int j = 0; j < nvoi; ++j)
			macro_stress[i] += ctan_lin[i * nvoi + j] * macro_strain[j];

	const double inv = get_inv_1(macro_stress);
	if (fabs(inv) > *inv_max)
		*inv_max = fabs(inv);

	return (fabs(inv) < inv_tol);
}


template <int tdim>
void micropp<tdim>::calc_ctan_lin()
{
	double *u_aux = (double *) malloc(nndim * sizeof(double));
	double *b = (double *) malloc(nndim * sizeof(double));
	double *du = (double *) malloc(nndim * sizeof(double));
	double *old = (double *) calloc(num_int_vars, sizeof(double));

	const int ns[3] = { nx, ny, nz };
	ell_matrix A;
	ell_init(&A, ell_cols, dim, dim, ns, CG_MIN_ERR, CG_MAX_ITS);

	for (int i = 0; i < nvoi; ++i) {

		double eps_1[nvoi] = { 0.0 };
		eps_1[i] += D_EPS_CTAN_AVE;

		double nr_err;
		newton_raphson(eps_1, &A, u_aux, b, du, old, &nr_err);

		double sig_1[6];
		calc_ave_stress(u_aux, old, sig_1);

		for (int v = 0; v < nvoi; ++v)
			ctan_lin[v * nvoi + i] = sig_1[v] / D_EPS_CTAN_AVE;
	}

	ell_free(&A);
	free(u_aux);
	free(b);
	free(du);
}


template <int tdim>
double micropp<tdim>::get_inv_1(const double *tensor) const
{
	assert(dim == 2 || dim == 3);

	const double ret = tensor[0] + tensor[1];

	if (dim == 2)
		return ret;

	return ret + tensor[2];
}


template <int tdim>
material_t micropp<tdim>::get_material(const int e) const
{
	int mat_num;
	if (micro_type == 0)
		if (elem_type[e] == 0)
			mat_num = 0;
		else
			mat_num = 1;

	else if (micro_type == 1)
		if (elem_type[e] == 0)
			mat_num = 0;
		else
			mat_num = 1;

	return material_list[mat_num];
}


template <int tdim>
void micropp<tdim>::get_elem_rhs(const double *u, const double *old,
                                 double be[npe * dim],
                                 int ex, int ey, int ez) const
{
	constexpr int npedim = npe * dim;
	double bmat[nvoi][npedim], stress_gp[nvoi], strain_gp[nvoi];

	memset(be, 0, npedim * sizeof(double));

	for (int gp = 0; gp < npe; ++gp) {

		calc_bmat(gp, bmat);

		get_strain(u, gp, strain_gp, ex, ey, ez);
		get_stress(gp, strain_gp, old, stress_gp, ex, ey, ez);

		for (int i = 0; i < npedim; ++i)
			for (int j = 0; j < nvoi; ++j)
				be[i] += bmat[j][i] * stress_gp[j] * wg;
	}
}


template <int tdim>
void micropp<tdim>::get_elem_nodes(int n[npe], int ex, int ey, int ez) const
{
	const int nxny = nx * ny;
	const int n0 = ez * nxny + ey * nx + ex;
	n[0] = n0;
	n[1] = n0 + 1;
	n[2] = n0 + nx + 1;
	n[3] = n0 + nx;

	if (dim == 3) {
		n[4] = n[0] + nxny;
		n[5] = n[1] + nxny;
		n[6] = n[2] + nxny;
		n[7] = n[3] + nxny;
	}
}


template<int tdim>
int micropp<tdim>::get_elem_type(int ex, int ey, int ez) const
{
	assert(micro_type == 0 || micro_type == 1);

	if (micro_type == 0) { // sphere in the center

		const double coor[3] = { ex * dx + dx / 2,
		                         ey * dy + dy / 2,
		                         ez * dz + dz / 2 }; // 2D -> dz = 0

		const double center[3] = { lx / 2,
		                           ly / 2,
		                           lz / 2 }; // 2D -> lz = 0

		const double rad = special_param;

		double tmp = 0.;
		for (int i = 0; i < dim; ++i)
			tmp += (center[i] - coor[i]) * (center[i] - coor[i]);

		return (tmp < rad * rad);

	} else if (micro_type == 1) { // 2 flat layers in y dir

		const double y = ey * dy + dy / 2;
		const double width = special_param;
		return (y < width);
	}

	cerr << "Invalid micro_type = " << micro_type << endl;
	return -1;
}


template <int tdim>
void micropp<tdim>::get_elem_displ(const double *u,
		double elem_disp[npe * dim],
		int ex, int ey, int ez) const
{
	int n[npe] ;
	get_elem_nodes(n, ex, ey, ez);

	for (int i = 0 ; i < npe; ++i)
		for (int d = 0; d < dim; ++d)
			elem_disp[i * dim + d] = u[n[i] * dim + d];
}


template <int tdim>
void micropp<tdim>::get_strain(const double *u, int gp, double *strain_gp,
		int ex, int ey, int ez) const
{
	double elem_disp[npe * dim];
	get_elem_displ(u, elem_disp, ex, ey, ez);

	double bmat[nvoi][npe * dim];
	calc_bmat(gp, bmat);

	memset(strain_gp, 0, nvoi * sizeof(double));
	for (int v = 0; v < nvoi; ++v)
		for (int i = 0; i < npe * dim; i++)
			strain_gp[v] += bmat[v][i] * elem_disp[i];
}


template <int tdim>
void micropp<tdim>::print_info() const
{
	printf("micropp%d\n", dim);
	printf("ngp %d n = [%d, %d, %d] => nn = %d\n", ngp, nx, ny, nz, nn);
	printf("l = [%lf, %lf, %lf]; width = %lf\n", lx, ly, lz, special_param);
}


template <int tdim>
void micropp<tdim>::get_stress(int gp, const double eps[nvoi],
                               const double *old,
                               double stress_gp[nvoi],
                               int ex, int ey, int ez) const
{
	const int e = glo_elem(ex, ey, ez);
	const material_t material = get_material(e);
	const double mu = material.mu;

	if (material.plasticity == true) {

		const double *eps_p_old = &old[intvar_ix(e, gp, 0)];
		const double alpha_old = old[intvar_ix(e, gp, 6)];

		plastic_get_stress(&material, eps, eps_p_old, alpha_old, stress_gp);

	} else {

		isolin_get_stress(&material, eps, stress_gp);
	}

}


template <int tdim>
void micropp<tdim>::calc_ave_stress(const double *u, const double *old,
                                    double *stress_ave) const
{
	memset(stress_ave, 0, nvoi * sizeof(double));

	for (int ez = 0; ez < nez; ++ez) { // 2D -> nez = 1
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex) {

				double stress_aux[nvoi] = { 0.0 };

				for (int gp = 0; gp < npe; ++gp) {

					double stress_gp[nvoi], strain_gp[nvoi];

					get_strain(u, gp, strain_gp, ex, ey, ez);
					get_stress(gp, strain_gp, old, stress_gp, ex, ey, ez);
					for (int v = 0; v < nvoi; ++v)
						stress_aux[v] += stress_gp[v] * wg;

				}
				for (int v = 0; v < nvoi; ++v)
					stress_ave[v] += stress_aux[v];
			}
		}
	}

	for (int v = 0; v < nvoi; ++v)
		stress_ave[v] /= vol_tot;
}


template <int tdim>
void micropp<tdim>::calc_ave_strain(const double *u,
		double strain_ave[nvoi]) const
{
	memset(strain_ave, 0, nvoi * sizeof(double));

	for (int ez = 0; ez < nez; ++ez) { // 2D -> nez = 1
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex) {

				double strain_aux[nvoi] = { 0.0 };

				for (int gp = 0; gp < npe; ++gp) {
					double strain_gp[nvoi];

					get_strain(u, gp, strain_gp, ex, ey, ez);
					for (int v = 0; v < nvoi; ++v)
						strain_aux[v] += strain_gp[v] * wg;
				}

				for (int v = 0; v < nvoi; v++)
					strain_ave[v] += strain_aux[v];
			}
		}
	}

	for (int v = 0; v < nvoi; v++)
		strain_ave[v] /= vol_tot;
}


template<int tdim>
void micropp<tdim>::calc_fields(const double *old, double *u) const
{
	for (int ez = 0; ez < nez; ++ez) { // 2D -> nez = 1
		for (int ey = 0; ey < ney; ++ey) {
			for (int ex = 0; ex < nex; ++ex) {

				double strain_aux[nvoi] = { 0.0 };
				double stress_aux[nvoi] = { 0.0 };

				for (int gp = 0; gp < npe; ++gp) {

					double stress_gp[nvoi], strain_gp[nvoi];

					get_strain(u, gp, strain_gp, ex, ey, ez);
					get_stress(gp, strain_gp, old, stress_gp, ex, ey, ez);

					for (int v = 0; v < nvoi; ++v) {
						strain_aux[v] += strain_gp[v] * wg;
						stress_aux[v] += stress_gp[v] * wg;
					}
				}

				const int e = glo_elem(ex, ey, ez);
				for (int v = 0; v < nvoi; ++v) {
					elem_strain[e * nvoi + v] = strain_aux[v] * ivol;
					elem_stress[e * nvoi + v] = stress_aux[v] * ivol;
				}
			}
		}
	}
}


// Explicit instantiation
template class micropp<2>;
template class micropp<3>;
