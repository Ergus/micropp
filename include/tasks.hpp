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

#endif

#pragma oss task out(_out[0])
template <typename T>
void d_set(T _in, T *_out)
{
	_out[0] = _in;
}

#pragma oss task out(_out[0]) out(_u_n[0; nndim])
template <int tdim>
void d_set_gp(double *const _int_vars_n, double *const _int_vars_k,
              double *const _u_n, double *const _u_k, int nndim,
              gp_t<tdim> *_out)
{
	_out[0].init(_int_vars_n, _int_vars_k, _u_n, _u_k, nndim);
}

#endif //TASKS_HPP
