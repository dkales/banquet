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

#include "bench_utils.h"

#include <cerrno>
#include <cstdio>
#if !defined(_MSC_VER)
#include <getopt.h>
#endif
#include <climits>
#include <cstdlib>

static bool parse_long(long *value, const char *arg) {
  errno = 0;
  const long v = strtol(arg, NULL, 10);

  if ((errno == ERANGE && (v == LONG_MAX || v == LONG_MIN)) ||
      (errno != 0 && v == 0)) {
    return false;
  }
  *value = v;

  return true;
}

static bool parse_uint32_t(uint32_t *value, const char *arg) {
  long tmp = 0;
  if (!parse_long(&tmp, arg)) {
    return false;
  }

  if (tmp < 0 || (unsigned long)tmp > UINT32_MAX) {
    return false;
  }

  *value = tmp;
  return true;
}

static void print_usage(const char *arg0) {
#if defined(_MSC_VER)
  printf("usage: %s iterations instance\n", arg0);
#else
  printf("usage: %s [-i iterations] instance\n", arg0);
#endif
}

static void print_usage_free(const char *arg0) {
#if defined(_MSC_VER)
  printf("usage: %s iterations instance\n", arg0);
#else
  printf("usage: %s [-i iterations] kappa N tau m1 m2 lambda\n", arg0);
#endif
}

bool parse_args(bench_options_t *options, int argc, char **argv) {
  if (argc <= 1) {
    print_usage(argv[0]);
    return false;
  }

  options->params = PARAMETER_SET_INVALID;
  options->iter = 10;

#if !defined(_MSC_VER)
  static const struct option long_options[] = {
      {"iter", required_argument, NULL, 'i'}, {0, 0, 0, 0}};

  int c = -1;
  int option_index = 0;

  while ((c = getopt_long(argc, argv, "i:l:", long_options, &option_index)) !=
         -1) {
    switch (c) {
    case 'i':
      if (!parse_uint32_t(&options->iter, optarg)) {
        printf("Failed to parse argument as positive base-10 number!\n");
        return false;
      }
      break;

    case '?':
    default:
      printf("usage: %s [-i iter] param\n", argv[0]);
      return false;
    }
  }

  if (optind == argc - 1) {
    uint32_t p = PARAMETER_SET_INVALID;
    if (!parse_uint32_t(&p, argv[optind])) {
      printf("Failed to parse argument as positive base-10 number!\n");
      return false;
    }

    if (p <= PARAMETER_SET_INVALID || p >= PARAMETER_SET_MAX_INDEX) {
      printf("Invalid parameter set selected!\n");
      return false;
    }
    options->params = static_cast<banquet_params_t>(p);
  } else {
    print_usage(argv[0]);
    return false;
  }
#else
  if (argc != 3) {
    print_usage(argv[0]);
    return false;
  }

  uint32_t p = PARAMETER_SET_INVALID;
  if (!parse_uint32_t(&options->iter, argv[1]) ||
      !parse_uint32_t(&p, argv[2])) {
    printf("Failed to parse argument as positive base-10 number!\n");
    return false;
  }

  if (p <= PARAMETER_SET_INVALID || p >= PARAMETER_SET_MAX_INDEX) {
    printf("Invalid parameter set selected!\n");
    return false;
  }
  options->params = p;
#endif

  return true;
}

bool parse_args_free(bench_options_free_t *options, int argc, char **argv) {
  if (argc <= 6) {
    print_usage_free(argv[0]);
    return false;
  }

  options->kappa = 0;
  options->tau = 0;
  options->N = 0;
  options->m1 = 0;
  options->m2 = 0;
  options->lambda = 0;
  options->iter = 10;

  static const struct option long_options[] = {
      {"iter", required_argument, NULL, 'i'}, {0, 0, 0, 0}};

  int c = -1;
  int option_index = 0;

  while ((c = getopt_long(argc, argv, "i:l:", long_options, &option_index)) !=
         -1) {
    switch (c) {
    case 'i':
      if (!parse_uint32_t(&options->iter, optarg)) {
        printf("Failed to parse argument as positive base-10 number!\n");
        return false;
      }
      break;

    case '?':
    default:
      print_usage_free(argv[0]);
      return false;
    }
  }

  uint32_t p = 0;
  if (!parse_uint32_t(&p, argv[optind++])) {
    printf("Failed to parse argument as positive base-10 number!\n");
    return false;
  }
  options->kappa = p;

  if (!parse_uint32_t(&p, argv[optind++])) {
    printf("Failed to parse argument as positive base-10 number!\n");
    return false;
  }
  options->N = p;

  if (!parse_uint32_t(&p, argv[optind++])) {
    printf("Failed to parse argument as positive base-10 number!\n");
    return false;
  }
  options->tau = p;

  if (!parse_uint32_t(&p, argv[optind++])) {
    printf("Failed to parse argument as positive base-10 number!\n");
    return false;
  }
  options->m1 = p;

  if (!parse_uint32_t(&p, argv[optind++])) {
    printf("Failed to parse argument as positive base-10 number!\n");
    return false;
  }
  options->m2 = p;

  if (!parse_uint32_t(&p, argv[optind++])) {
    printf("Failed to parse argument as positive base-10 number!\n");
    return false;
  }
  options->lambda = p;
  return true;
}
