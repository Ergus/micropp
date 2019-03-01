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

#ifndef UTIL_HPP
#define UTIL_HPP

// Debug print macro.
#ifdef NDEBUG
#define dbprintf(...)
#else
#define dbprintf(...) fprintf(stderr, __VA_ARGS__)
#endif

#ifdef __cplusplus

#include <iostream>
#include <vector>
#include <fstream>
#include <numeric>

#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cmath>
#include <cstdio>
#include <ctime>
#include <cstring>
#include <cassert>


using namespace std;


inline uint64_t devest(const vector<uint64_t>  &in, const uint64_t mean)
{
	uint64_t out = 0;
	for (const auto &x : in) {
		const uint64_t tmp = (x - mean);
		out += tmp * tmp;
	}

	return sqrt(out / in.size());
}


inline void filter(double *arr, int n, double rel_tol)
{
	#ifdef FILTER
	double max = arr[0];
	for (int i = 1; i < n; ++i)
		if (arr[i] > max)
			max = arr[i];

	for (int i = 0; i < n; ++i)
		arr[i] = (fabs(arr[i]) > max * rel_tol) ? arr[i] : 0.0;
	#endif
}


constexpr int mypow(int v, int e)
{
	return (e == 0) ? 1 : v * mypow(v, e - 1);
}


inline void print_matvec(const double *vec, const int rows, const int cols,
			 FILE *out = NULL)
{
	for (int i = 0; i < rows; ++i) {
		fprintf(out, "\n[ ");
		for (int j = 0; j < cols; ++j)
			fprintf(out, "%lf, ", vec[i * cols + j]);
		fprintf(out, "]\n");
	}
}

inline void print_vec(const double *vec, int n, const char file_name[] = "")
{
	FILE *file = NULL;
	if (strlen(file_name) == 0)
		file = fopen(file_name, "w");
	else
		file = stdout;

	assert(file != NULL);

	for (int i = 0; i < n; ++i)
		fprintf(file, "[%lf]\n", vec[i]);

	if (file && file != stdout)
		fclose(file);
}

template <int tdim>
void mvp(const double m[tdim][tdim], const double x[tdim], double y[tdim])
{
	for (int i = 0; i < tdim; ++i) {
		double tmp = 0.0;
		for (int j = 0; j < tdim; ++j)
			tmp += m[i][j] * x[j];
		y[i] = tmp;
	}
}

inline void distribute(int total, int nodes, int *start, int *nelems)
{
	assert (nodes > 0);
	const int frac = total / nodes;
	const int mod = total - frac * nodes;
	int cum = 0;
	for(int i = 0; i < mod; ++i) {
		nelems[i] = frac + 1;
		start[i] = cum;
		cum += frac + 1;
	}

	for(int i = mod; i < nodes; ++i) {
		nelems[i] = frac;
		start[i] = cum;
		cum += frac;
	}
}

template <typename T>
void printarray(const int size, T *array)
{
	cerr << "[ ";
	for (int i = 0; i < size; ++i)
		cerr << array[i] << " ";
	cerr << "]\n";
}

#endif
