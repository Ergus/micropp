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

#ifdef __cplusplus


#include <vector>
#include <iostream>
#include <fstream>
#include <iomanip>

#include <cmath>
#include <cassert>
#include <cstring>
#include <ctime>

#endif

#include "util.hpp"
#include "ell.hpp"
#include "material.hpp"
#include "instrument.hpp"
#include "newton.hpp"

#include "gp.hpp"
#include "tasks.hpp"

#define MAX_DIM         3
#define MAX_MATS        10
#define NUM_VAR_GP      7  // eps_p_1 (6) , alpha_1 (1)

#define CG_MIN_ERR      1.0e-50
#define CG_MAX_ITS      1000
#define CG_REL_ERR      1.0e-5

#define FILTER_REL_TOL  1.0e-5

#define D_EPS_CTAN      1.0e-8
#define D_EPS_CTAN_AVE  1.0e-8

#define CONSTXG         0.577350269189626
#define SQRT_2DIV3      0.816496581

#define glo_elem(ex,ey,ez)   ((ez) * (nx-1) * (ny-1) + (ey) * (nx-1) + (ex))
#define intvar_ix(e,gp,var)  ((e) * npe * NUM_VAR_GP + (gp) * NUM_VAR_GP + (var))


enum {
	MIC_SPHERE,
	MIC_LAYER_Y,
	MIC_CILI_FIB_Z,
	MIC_CILI_FIB_XZ,
	MIC_QUAD_FIB_XYZ,
	MIC_QUAD_FIB_XZ,
	MIC_QUAD_FIB_XZ_BROKEN_X,
	MIC_SPHERES
};

enum {
	NO_COUPLING,
	ONE_WAY,
	FULL
};

using namespace std;

template <int tdim>
class micropp {

	protected:
		static constexpr int dim = tdim;                  // 2, 3
		static constexpr int npe = mypow(2, dim);         // 4, 8
		static constexpr int nvoi = dim * (dim + 1) / 2;  // 3, 6

		// Constant members
		int ngp, nx, ny, nz, nn, nndim;
		int nex, ney, nez, nelem;
		double lx, ly, lz;
		double dx, dy, dz;
		double vol_tot;
		double special_param, wg, ivol;

		int micro_type, num_int_vars;
		bool copy;

		// Non constant members
		gp_t<tdim> *gp_list;

		int coupling;

		double micro_params[5];
		int numMaterials;
		material_t *material_list;
		double *ctan_lin;

		ell_matrix A0; // Linear Jacobian (constant)

		int *ell_cols;
		int ell_cols_size;

		int *elem_type;
		double *elem_stress;
		double *elem_strain;
		double *vars_new_aux;
		double f_trial_max;

		double *du_n, *du_k, *dint_vars_n, *dint_vars_k;

	protected:
		void calc_ctan_lin();
		material_t get_material(const int e) const;

		void get_elem_nodes(int n[npe], int ex, int ey, int ez = 0) const;

		void get_elem_displ(const double *u, double elem_disp[npe * dim],
				    int ex, int ey, int ez = 0) const;

		void get_strain(const double *u, int gp, double strain_gp[nvoi],
				int ex, int ey, int ez = 0) const;
		void get_stress(int gp, const double eps[nvoi], const double *vars_old,
				double stress_gp[nvoi], int ex, int ey, int ez = 0) const;

		int get_elem_type(int ex, int ey, int ez = 0) const;
		void get_elem_rhs(const double *u, const double *vars_old,
				  double be[npe * dim], int ex, int ey, int ez = 0) const;

		int newton_raphson(ell_matrix *A, double *b, double *u, double *du,
				   const bool non_linear, const double strain[nvoi],
				   const double *vars_old, newton_t *newton);

		void calc_ave_stress(double *u, double *vars_old,
				     double stress_ave[nvoi]) const;
		void calc_ave_strain(const double *u, double strain_ave[nvoi]) const;
		void calc_fields(double *u, double *vars_old);
		void calc_bmat(int gp, double bmat[nvoi][npe * dim]) const;
		bool calc_vars_new(const double *u, double *vars_old, double *vars_new,
				   double *f_trial_max);

		void get_elem_mat(const double *u, const double *vars_old,
				  double Ae[npe * dim * npe * dim],
				  int ex, int ey, int ez = 0) const;

		void set_displ_bc(const double strain[nvoi], double *u) const;

		double assembly_rhs(const double *u, const double *vars_old, double *b) const;
		void assembly_mat(ell_matrix *A, const double *u,
				  const double *vars_old) const;

		void write_vtu(double *u, const double *vars_old, const char *filename);

		void get_dev_tensor(const double tensor[6], double tensor_dev[6]) const;

		void plastic_get_stress(const material_t *material, const double eps[6],
					const double *eps_p_old, const double *alpha_old,
					double stress[6]) const;

		bool plastic_law(const material_t *material, const double eps[6],
				 const double *eps_p_old, const double *alpha_old,
				 double *_dl, double _normal[6], double _s_trial[6],
				 double *_f_trial) const;
		void plastic_get_ctan(const material_t *material, const double eps[nvoi],
				      const double *eps_p_old, const double *alpha_old,
				      double ctan[nvoi][nvoi]) const;
		bool plastic_evolute(const material_t *material, const double eps[6],
				     const double *eps_p_old, const double *alpha_old,
				     double eps_p_new[6], double *alpha_new, double *f_trial) const;


		void isolin_get_ctan(const material_t *material, double ctan[nvoi][nvoi]) const;
		void isolin_get_stress(const material_t *material, const double eps[6],
				       double stress[6]) const;

	public:

		static void* operator new(std::size_t sz);
		static void operator delete(void *p);

		micropp() = delete;

		micropp(const micropp &);

		micropp(const int ngp, const int size[3], const int micro_type,
			const double *micro_params, const material_t *materials,
			int coupling = ONE_WAY);

		~micropp();

		/* The most important functions */
		void set_macro_strain(const int gp_id, const double *macro_strain);
		void get_macro_stress(const int gp_id, double *macro_stress) const;
		void get_macro_ctan(const int gp_id, double *macro_ctan) const;
		void homogenize();

		/* Extras */
		int is_non_linear(const int gp_id) const;
		int get_non_linear_gps(void) const;
		double get_f_trial_max(void) const;

		int get_nl_flag(const int gp_id) const;

		void get_sigma_solver_its(int gp_id, int sigma_solver_err[NR_MAX_ITS]) const;
		void get_sigma_solver_err(int gp_id, double sigma_solver_err[NR_MAX_ITS]) const;
		void get_sigma_newton_err(int gp_id, double sigma_nr_err[NR_MAX_ITS]) const;

		int get_sigma_newton_its(int gp_id) const;

		int get_cost(int gp_id) const;

		void output(int gp_id, const char *filename);
		void update_vars();
		void print_info() const;


		void homogenize_conditional_task(const int nvoi,
						 int *ell_cols, const int ell_cols_size,
						 const material_t *material_list, const int numMaterials,
						 int *elem_type, int nelem,
						 gp_t<tdim> *gp_ptr,
						 int nndim, int num_int_vars);


};


#endif
