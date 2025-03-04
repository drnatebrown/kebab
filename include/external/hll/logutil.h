/*
Taken from: https://github.com/mindis/hll

MIT License

Copyright (c) 2016 Daniel Baker

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#ifndef _LOG_UTIL_H__
#define _LOG_UTIL_H__
#ifndef __STDC_FORMAT_MACROS
#  define __STDC_FORMAT_MACROS
#endif
#include <cinttypes>
#include <cstdlib>
#include <cstdio>
#include <cstdarg>

#define _FUNCTION_MACRO_ __PRETTY_FUNCTION__
#define LOG_INFO(...) log_info(__func__, ##__VA_ARGS__);
#define LOG_WARNING(...) log_warning(_FUNCTION_MACRO_, ##__VA_ARGS__);
#define LOG_EXIT(...) log_error(_FUNCTION_MACRO_, __LINE__, ##__VA_ARGS__);
#if !NDEBUG
#    define LOG_DEBUG(...) log_debug(_FUNCTION_MACRO_, __LINE__, ##__VA_ARGS__);
#    define LOG_ASSERT(condition) log_assert(_FUNCTION_MACRO_, __LINE__, ((uint64_t)(condition)), (#condition))
#else
#    define LOG_DEBUG(...)
#    define LOG_ASSERT(...)
#endif

static inline void log_debug(const char *func, int line, const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    fprintf(stderr, "[D:%s:%d] ", func, line);
    vfprintf(stderr, fmt, args);
    va_end(args);
}

static inline void log_warning(const char *func, const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    fprintf(stderr, "[W:%s] ", func);
    vfprintf(stderr, fmt, args);
    va_end(args);
}

static inline void log_info(const char *func, const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    fprintf(stderr, "[%s] ", func);
    vfprintf(stderr, fmt, args);
    va_end(args);
}

static inline void log_error(const char *func, int line, const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    fprintf(stderr, "[E:%s:%d] ", func, line);
    vfprintf(stderr, fmt, args);
    va_end(args);
    exit(EXIT_FAILURE);
}

static inline void log_assert(const char *func, int line, int assertion, const char *assert_str) {
    if(!assertion) {
        fprintf(stderr, "[E:%s:%d] Assertion '%s' failed.",
                func, line, assert_str);
        exit(EXIT_FAILURE);
    }
}

#endif /* #ifndef _LOG_UTIL_H__ */
