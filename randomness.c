/*
 *  This file is part of the optimized implementation of the Picnic signature
 *  scheme.
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

#include "randomness.h"

#if defined(__linux__) && (__GLIBC__ > 2 || __GLIBC_MINOR__ >= 25)
#include <sys/random.h>

int rand_bytes(uint8_t *dst, size_t len) {
  const ssize_t ret = getrandom(dst, len, GRND_NONBLOCK);
  if (ret < 0 || (size_t)ret != len) {
    return 0;
  }
  return 1;
}
#elif defined(__APPLE__) && defined(HAVE_APPLE_FRAMEWORK)
#include <Security/Security.h>

int rand_bytes(uint8_t *dst, size_t len) {
  if (SecRandomCopyBytes(kSecRandomDefault, len, dst) == errSecSuccess) {
    return 1;
  }
  return 0;
}
#elif defined(__linux__) || defined(__APPLE__)
#include <errno.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#if defined(__linux__)
#include <linux/random.h>
#include <sys/ioctl.h>
#endif

#if !defined(O_NOFOLLOW)
#define O_NOFOLLOW 0
#endif
#if !defined(O_CLOEXEC)
#define O_CLOEXEC 0
#endif

int rand_bytes(uint8_t *dst, size_t len) {
  int fd;
  while ((fd = open("/dev/urandom", O_RDONLY | O_NOFOLLOW | O_CLOEXEC, 0)) ==
         -1) {
    // check if we should restart
    if (errno != EINTR) {
      return 0;
    }
  }
#if O_CLOEXEC == 0
  fcntl(fd, F_SETFD, fcntl(fd, F_GETFD) | FD_CLOEXEC);
#endif

#if defined(__linux__)
  int cnt = 0;
  if (ioctl(fd, RNDGETENTCNT, &cnt) == -1) {
    // not ready
    close(fd);
    return 0;
  }
#endif

  while (len) {
    const ssize_t ret = read(fd, dst, len);
    if (ret == -1) {
      if (errno == EAGAIN || errno == EINTR) {
        // retry
        continue;
      }
      close(fd);
      return 0;
    }

    dst += ret;
    len -= ret;
  }

  close(fd);
  return 1;
}
#elif defined(_WIN16) || defined(_WIN32) || defined(_WIN64)
#include <windows.h>

int rand_bytes(uint8_t *dst, size_t len) {
  if (len > ULONG_MAX) {
    return 0;
  }
  if (!BCRYPT_SUCCESS(BCryptGenRandom(NULL, dst, (ULONG)len,
                                      BCRYPT_USE_SYSTEM_PREFERRED_RNG))) {
    return 0;
  }
  return 1;
}
#else
#error "Unsupported OS! Please implement rand_bytes."
#endif
