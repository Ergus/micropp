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

#pragma oss task inout(_out[0]) label(init_gp)
template <int tdim>
void set_gp(double *_int_vars_k, double *_u_k, int tnndim, gp_t<tdim> *_out)
{
	dprintf("Node: %d set_gp(%p) = {%p; %p}\n",
	        get_node_id(), _out, _int_vars_k, _u_k);
	_out[0].init(_int_vars_k, _u_k, tnndim);
}


#pragma oss task out(_out[0])
template <typename T>
void set_val(T _in, T *_out)
{
	*_out = _in;
}


template <int tdim>
void homogenize_conditional_task(struct data self, int nvoi,
                                 int *ell_cols, const int ell_cols_size,
                                 const material_t *material_list, const int numMaterials,
                                 int *elem_type, int nelem,
                                 gp_t<tdim> *gp_ptr,
                                 int nndim, int num_int_vars,
                                 const bool allocated);


template <int tdim>
void homogenize_weak_task(data self, int nvoi,
                          int *ell_cols, const int ell_cols_size,
                          const material_t *material_list, const int numMaterials,
                          int *elem_type, int nelem,
                          gp_t<tdim> *gp_ptr, int nndim, int num_int_vars);


#endif //TASKS_HPP
