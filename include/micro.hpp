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

#ifndef MICRO_HPP
#define MICRO_HPP

#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <cmath>
#include <cassert>
#include <cstring>

#include "util.hpp"
#include "ell.hpp"
#include "material.hpp"
#include "gp.hpp"
#include "instrument.hpp"
#include "tasks.hpp"

#define MAX_DIM         3
#define MAX_MATS        10
#define NUM_VAR_GP      7  // eps_p_1 (6) , alpha_1 (1)

#define CG_MIN_ERR      1.0e-8
#define CG_MAX_ITS      2000
#define NR_MAX_TOL      1.0e-5
#define NR_MAX_ITS      40
#define D_EPS_CTAN      1.0e-8
#define D_EPS_CTAN_AVE  1.0e-8

#define CONSTXG         0.577350269189626

#define glo_elem(ex,ey,ez)   ((ez) * (nx-1) * (ny-1) + (ey) * (nx-1) + (ex))
#define intvar_ix(e,gp,var)  ((e) * 8 * NUM_VAR_GP + (gp) * NUM_VAR_GP + (var))

using namespace std;

struct data {
	protected:
		data(const int tdim, const int ngp, const int size[3],
		     const int micro_type, const double *micro_params,
		     const material_t *materials);

	public:
		// Constants only vars
		int ngp, nx, ny, nz, nn, nndim;
		int nex, ney, nez, nelem;
		double lx, ly, lz;
		double dx, dy, dz;
		double vol_tot;
		double special_param, inv_tol, wg, ivol;

		int micro_type, num_int_vars;

		// Other variables
		bool output_files_header;

		double ctan_lin[36]; // nvoi * nvoi

		int numMaterials;
		material_t *material_list;

		int *elem_type;
		double *elem_stress;
		double *elem_strain;

		// Nanos stuff
		int *ell_cols;
		int ell_cols_size;

		double *du_k, *dint_vars_k;
};

template <int tdim>
class micropp : public data {
	private:
	public:
		// Variables (static constexpr)
		static constexpr int dim = tdim;                   // 2, 3
		static constexpr int npe = mypow(2, tdim);         // 4, 8
		static constexpr int nvoi = tdim * (tdim + 1) / 2;  // 3, 6
		const bool copy;

		gp_t<tdim> *gp_list;

		// Common
		void calc_ctan_lin();

		material_t get_material(const int e) const;

		void get_strain(const double *u, int gp,
						double strain_gp[nvoi],
						int ex, int ey, int ez = 0) const;

		void get_elem_nodes(int n[npe], int ex, int ey, int ez = 0) const;

		void get_elem_displ(const double *u, double elem_disp[npe * dim],
							int ex, int ey, int ez = 0) const;

		void get_stress(int gp, const double eps[nvoi], const double *old,
		                double stress_gp[nvoi],
		                int ex, int ey, int ez = 0) const;

		int get_elem_type(int ex, int ey, int ez = 0) const;

		void get_elem_rhs(const double *u, const double *old,
		                  double be[npe * dim],
		                  int ex, int ey, int ez = 0) const;

		void calc_ave_stress(const double *u, const double *out,
		                     double stress_ave[nvoi]) const;

		void calc_ave_strain(const double *u, double strain_ave[nvoi]) const;


		void calc_fields(const double *old, double *u) const;

		int newton_raphson(const double strain[nvoi], ell_matrix *A, double *u,
		                   double *b, double *du, double *old,
		                   double *_err) const;

		// Specialized
		template <typename... Rest>
		void get_elem_mat(const double *u, const double *old,
		                  double Ae[npe * dim * npe * dim],
		                  int ex, int ey, Rest...) const;

		void set_displ_bc(const double strain[nvoi], double *u) const;

		double assembly_rhs(const double *u, const double *old, double *b) const;
		void assembly_mat(const double *u, const double *old, ell_matrix *A) const;

		void calc_bmat(int gp, double bmat[nvoi][npe * dim]) const;

		bool calc_vars_new(const double *u, const double *_old, double *_new) const;

		void write_vtu(const double *old, double *u, int tstep, int gp_id);

		// Functions Only for 3D
		void get_dev_tensor(const double tensor[6], double tensor_dev[6]) const;

		bool plastic_law(
			const material_t *material, const double eps[6],
			const double eps_p_old[6], double alpha_old,
			double *_dl, double _normal[6], double _s_trial[6]) const;

		void plastic_get_ctan(
			const material_t *material, const double eps[6],
			const double eps_p_old[6], double alpha_old,
			double ctan[6][6]) const;

		bool plastic_evolute(const material_t *material, const double eps[6],
			const double eps_p_old[6], double alpha_old,
			double eps_p_new[6], double *alpha_new) const;

		void isolin_get_ctan(const material_t *material, double ctan[6][6]) const;


		void plastic_get_stress(const material_t *material, const double eps[6],
		                        const double eps_p_old[6], double alpha_old,
		                        double stress[6]) const;

		void isolin_get_stress(const material_t *material, const double eps[6],
		                       double stress[6]) const;

		micropp() = delete;

		micropp(const int ngp, const int size[3], const int micro_type,
		        const double *micro_params, const material_t *materials);

		micropp(data *in, gp_t<tdim> *list);

		~micropp();

		// common Functions
		int get_nl_flag(const int gp_id) const;
		void set_macro_strain(const int gp_id, const double *macro_strain);
		void get_macro_stress(const int gp_id, double *macro_stress) const;
		void get_macro_ctan(const int gp_id, double *macro_ctan) const;
		void homogenize();
		void output(int tstep, int gp_id);
		void write_info_files();
		void update_vars();
		void print_info() const;
};

#endif // MICRO_HPP
