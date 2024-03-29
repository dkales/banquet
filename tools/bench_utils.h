/*
 *  Part of this file is is based on the optimized implementation of the Picnic
 *  signature scheme.
 *  See the accompanying documentation for complete details.
 *  The code is provided under the MIT license:
 *
 * Copyright (c) 2019-2020 Sebastian Ramacher, AIT
 * Copyright (c) 2016-2020 Graz University of Technology
 * Copyright (c) 2017 Angela Promitzer

 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the ""Software""), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED *AS IS*, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#ifndef BENCH_UTILS_H
#define BENCH_UTILS_H

#include <stdbool.h>
#include <stdint.h>

#include "../banquet_instances.h"

typedef struct {
  banquet_params_t params;
  uint32_t iter;
} bench_options_t;

typedef struct {
  uint32_t iter;
  uint32_t kappa;
  uint32_t m1;
  uint32_t m2;
  uint32_t tau;
  uint32_t N;
  uint32_t lambda;
} bench_options_free_t;

bool parse_args(bench_options_t *options, int argc, char **argv);

bool parse_args_free(bench_options_free_t *options, int argc, char **argv);

#endif
