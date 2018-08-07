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

#ifndef TASKS_HPP
#define TASKS_HPP

#include <cstdlib>
#include <cassert>

#include "micro.hpp"
#include "ell.hpp"

#ifdef NANOS6

#include "nanos6.h"

static inline void *rrd_malloc(size_t size)
{
	dprintf("Using dmalloc\n");
	void *ret = nanos_dmalloc(size, DMALLOC_RR, 0, NULL);
	assert(ret != NULL);

	return ret;
}

static inline void rrd_free(void *in)
{
	dprintf("Using dfree\n");
	nanos_dfree(in);
}

#define get_node_id() nanos_get_node_id()
#define get_nodes_nr() nanos_get_nodes_nr()

#else

static inline void *rrd_malloc(size_t size)
{
	dprintf("Using libc malloc\n");
	void *ret = malloc(size);
	assert(ret != NULL);
	return ret;
}

static inline void rrd_free(void *in)
{
	dprintf("Using libc_free\n");
	free(in);
}

#define get_node_id() 0
#define get_nodes_nr() 1

#endif


template <int tdim>
void set_gp(double *const _int_vars_n, double *const _int_vars_k,
              double *const _u_n, double *const _u_k, int nndim,
              gp_t<tdim> *_out)
{
       _out[0].init(_int_vars_n, _int_vars_k, _u_n, _u_k, nndim);
}


template <int tdim>
void homogenize_conditional_task(struct data self, int nvoi,
                                 int *ell_cols, const int ell_cols_size,
                                 const material_t *material_list, const int numMaterials,
                                 int *elem_type, int nelem,
                                 gp_t<tdim> *gp_ptr,
                                 double *u_k, double *u_n, int nndim,
                                 const bool allocated, double *vars_n_old,
                                 double *vars_k_new, int num_int_vars);


#pragma oss task in(ell_cols[0; ell_cols_size]) \
	in(material_list[0; numMaterials]) \
	in(elem_type[0; nelem]) \
	 \
	inout(gp_ptr[0]) \
	inout(u_n[0; nndim]) \
	inout(u_k[0; nndim]) \
	inout(vars_n_old[0; num_int_vars]) \
	inout(vars_k_new[0; num_int_vars])
template <int tdim>
void homogenize_weak_task(data self, int nvoi,
                          int *ell_cols, const int ell_cols_size,
                          const material_t *material_list, const int numMaterials,
                          int *elem_type, int nelem,
                          gp_t<tdim> *gp_ptr,
                          double *u_k, double *u_n, int nndim,
                          double *vars_n_old, double *vars_k_new, int num_int_vars);

#endif //TASKS_HPP
