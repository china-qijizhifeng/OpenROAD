// Run a command and count user-space hardware events for it and its threads.
#define _GNU_SOURCE
#include <errno.h>
#include <linux/perf_event.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <sys/ioctl.h>
#include <sys/syscall.h>
#include <sys/wait.h>
#include <unistd.h>

struct Event
{
  const char* name;
  uint32_t type;
  uint64_t config;
};

#define CACHE(cache, op, result) \
  ((cache) | ((op) << 8) | ((result) << 16))

static const struct Event events[] = {
    {"cycles", PERF_TYPE_HARDWARE, PERF_COUNT_HW_CPU_CYCLES},
    {"instructions", PERF_TYPE_HARDWARE, PERF_COUNT_HW_INSTRUCTIONS},
    {"branch_misses", PERF_TYPE_HARDWARE, PERF_COUNT_HW_BRANCH_MISSES},
    {"l1d_read_miss",
     PERF_TYPE_HW_CACHE,
     CACHE(PERF_COUNT_HW_CACHE_L1D,
           PERF_COUNT_HW_CACHE_OP_READ,
           PERF_COUNT_HW_CACHE_RESULT_MISS)},
    {"llc_read_access",
     PERF_TYPE_HW_CACHE,
     CACHE(PERF_COUNT_HW_CACHE_LL,
           PERF_COUNT_HW_CACHE_OP_READ,
           PERF_COUNT_HW_CACHE_RESULT_ACCESS)},
    {"llc_read_miss",
     PERF_TYPE_HW_CACHE,
     CACHE(PERF_COUNT_HW_CACHE_LL,
           PERF_COUNT_HW_CACHE_OP_READ,
           PERF_COUNT_HW_CACHE_RESULT_MISS)},
    {"dtlb_read_miss",
     PERF_TYPE_HW_CACHE,
     CACHE(PERF_COUNT_HW_CACHE_DTLB,
           PERF_COUNT_HW_CACHE_OP_READ,
           PERF_COUNT_HW_CACHE_RESULT_MISS)},
};
#define COUNT (sizeof(events) / sizeof(events[0]))

int main(int argc, char** argv)
{
  if (argc < 3) {
    fprintf(stderr, "usage: pcount OUT.json command...\n");
    return 2;
  }
  int ready[2];
  if (pipe(ready) != 0) {
    return 2;
  }
  const pid_t child = fork();
  if (child == 0) {
    char byte;
    close(ready[1]);
    if (read(ready[0], &byte, 1) != 1) {
      _exit(127);
    }
    execvp(argv[2], argv + 2);
    _exit(127);
  }
  close(ready[0]);
  int fds[COUNT];
  int errors[COUNT];
  for (size_t i = 0; i < COUNT; ++i) {
    struct perf_event_attr attr;
    memset(&attr, 0, sizeof(attr));
    attr.size = sizeof(attr);
    attr.type = events[i].type;
    attr.config = events[i].config;
    attr.disabled = 1;
    attr.enable_on_exec = 1;
    attr.inherit = 1;
    attr.exclude_kernel = 1;
    attr.exclude_hv = 1;
    attr.read_format
        = PERF_FORMAT_TOTAL_TIME_ENABLED | PERF_FORMAT_TOTAL_TIME_RUNNING;
    fds[i] = (int) syscall(SYS_perf_event_open, &attr, child, -1, -1, 0);
    errors[i] = fds[i] < 0 ? errno : 0;
  }
  if (write(ready[1], "x", 1) != 1) {
    return 2;
  }
  close(ready[1]);
  int status = 0;
  waitpid(child, &status, 0);
  FILE* out = fopen(argv[1], "w");
  fprintf(out, "{");
  for (size_t i = 0; i < COUNT; ++i) {
    uint64_t value[3] = {0, 0, 0};
    double scaled = -1;
    if (fds[i] >= 0 && read(fds[i], value, sizeof(value)) == sizeof(value)
        && value[2] > 0) {
      scaled = (double) value[0] * (double) value[1] / (double) value[2];
    }
    fprintf(out,
            "%s\"%s\": {\"value\": %.0f, \"enabled\": %llu, \"running\": %llu, \"errno\": %d}",
            i ? ", " : "",
            events[i].name,
            scaled,
            (unsigned long long) value[1],
            (unsigned long long) value[2],
            errors[i]);
  }
  fprintf(out, ", \"exit_status\": %d}\n", status);
  fclose(out);
  return WIFEXITED(status) ? WEXITSTATUS(status) : 1;
}
