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

#include <cstdio>
#include <cstdlib>
#include <cassert>

#include "ell.hpp"
#include "util.hpp"

#ifdef NANOS6

#include "nanos6.h"

static inline void *rrd_malloc(size_t size)
{
	void *ret = nanos6_dmalloc(size, nanos6_equpart_distribution, 1, nullptr);
	assert(ret != NULL);
	dbprintf("Using nanos6_dmalloc [%p -> %p] size %d\n",
	         ret, (char*)ret + size, size);

	return ret;
}

static inline void rrd_free(void *in)
{
	dbprintf("Using nanos6_dfree(%p)\n", in);
	nanos6_dfree(in);
}

static inline void *rrl_malloc(size_t size)
{
	void *ret = nanos6_lmalloc(size);
	assert(ret != NULL);
	dbprintf("Using nanos6_lmalloc [%p -> %p] size %d\n",
		 ret, (char*)ret + size, size);

	return ret;
}

static inline void *rrl_calloc(size_t nmemb, size_t size)
{
	const int bytes = nmemb * size;

	void *ret = nanos6_lmalloc(bytes);
	assert(ret != NULL);
	dbprintf("Using nanos6_lmalloc [%p -> %p] size %d\n",
	         ret, (char*)ret + size, size);

	memset(ret, 0, bytes);

	return ret;
}


static inline void rrl_free(void *in)
{
	dbprintf("Using nanos6_lfree(%p)\n", in);
	nanos6_lfree(in, 0);
}


#define get_node_id() nanos6_get_cluster_node_id()
#define get_nodes_nr() nanos6_get_cluster_nodes()

#else

static inline void *rrd_malloc(size_t size)
{
	void *ret = malloc(size);
	assert(ret != NULL);
	dbprintf("Using libc malloc [%p -> %p] size %d\n",
	         ret, (char*)ret + size, size);

	return ret;
}

static inline void rrd_free(void *in)
{
	dbprintf("Using libc_free(%p)\n", in);
	free(in);
}

static inline void *rrl_malloc(size_t size)
{
	void *ret = malloc(size);
	assert(ret != NULL);
	dbprintf("Using libc_lmalloc [%p -> %p] size %d\n",
	         ret, (char*)ret + size, size);

	return ret;
}

static inline void *rrl_calloc(size_t nmemb, size_t size)
{
	void *ret = calloc(nmemb, size);
	assert(ret != NULL);
	dbprintf("Using libc_lcalloc [%p -> %p] size %d\n",
	         ret, (char*)ret + size, size);

	return ret;
}

static inline void rrl_free(void *in)
{
	dbprintf("Using libc_lfree(%p)\n", in);
	free(in);
}

#define get_node_id() 0
#define get_nodes_nr() 1

#endif

#endif //TASKS_HPP
