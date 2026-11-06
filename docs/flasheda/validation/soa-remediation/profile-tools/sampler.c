// LD_PRELOAD wall-clock stack sampler for the main thread.
// Each sample is a fixed record: [depth, pc0, pc1, ...] written with write(2).
#define _GNU_SOURCE
#include <errno.h>
#include <execinfo.h>
#include <fcntl.h>
#include <signal.h>
#include <stdlib.h>
#include <string.h>
#include <sys/syscall.h>
#include <time.h>
#include <ucontext.h>
#include <unistd.h>

#define DEPTH 64
static int out_fd = -1;

static void handler(int sig, siginfo_t* info, void* context)
{
  (void) sig;
  (void) info;
  const int saved = errno;
  void* frames[DEPTH + 8];
  void* record[DEPTH + 1];
  const int n = backtrace(frames, DEPTH + 8);
  const ucontext_t* uc = (const ucontext_t*) context;
  void* pc = (void*) uc->uc_mcontext.gregs[REG_RIP];
  int start = -1;
  for (int i = 0; i < n && i < 8; ++i) {
    if (frames[i] == pc) {
      start = i;
      break;
    }
  }
  int depth = 0;
  if (start < 0) {
    record[1] = pc;
    depth = 1;
    start = n < 3 ? n : 3;
  }
  for (int i = start; i < n && depth < DEPTH; ++i) {
    record[1 + depth++] = frames[i];
  }
  record[0] = (void*) (long) depth;
  if (write(out_fd, record, sizeof(record)) < 0) {
    // Nothing safe to do inside a signal handler.
  }
  errno = saved;
}

static void copy_maps(const char* path)
{
  char name[4096];
  strncpy(name, path, sizeof(name) - 6);
  name[sizeof(name) - 6] = 0;
  strcat(name, ".maps");
  const int in = open("/proc/self/maps", O_RDONLY);
  const int out = open(name, O_WRONLY | O_CREAT | O_TRUNC, 0644);
  char buffer[65536];
  ssize_t got;
  while (in >= 0 && out >= 0 && (got = read(in, buffer, sizeof(buffer))) > 0) {
    if (write(out, buffer, got) != got) {
      break;
    }
  }
  close(in);
  close(out);
}

__attribute__((constructor)) static void start_sampler(void)
{
  const char* path = getenv("SAMPLER_OUT");
  if (path == NULL || *path == 0) {
    return;
  }
  char copy[4096];
  strncpy(copy, path, sizeof(copy) - 1);
  copy[sizeof(copy) - 1] = 0;
  // Child processes must not overwrite the parent's sample file.
  unsetenv("SAMPLER_OUT");
  out_fd = open(copy, O_WRONLY | O_CREAT | O_TRUNC, 0644);
  if (out_fd < 0) {
    return;
  }
  void* warm[4];
  backtrace(warm, 4);  // Load the unwinder before the first signal.
  struct sigaction action;
  memset(&action, 0, sizeof(action));
  action.sa_sigaction = handler;
  action.sa_flags = SA_SIGINFO | SA_RESTART;
  sigemptyset(&action.sa_mask);
  sigaction(SIGPROF, &action, NULL);
  struct sigevent event;
  memset(&event, 0, sizeof(event));
  event.sigev_notify = SIGEV_THREAD_ID;
  event.sigev_signo = SIGPROF;
  event._sigev_un._tid = (pid_t) syscall(SYS_gettid);
  timer_t timer;
  if (timer_create(CLOCK_MONOTONIC, &event, &timer) != 0) {
    return;
  }
  const long interval_ns = 1000000;  // 1 kHz wall clock, main thread.
  struct itimerspec spec = {{0, interval_ns}, {0, interval_ns}};
  timer_settime(timer, 0, &spec, NULL);
  // Written after the timer starts: libraries loaded later are rare here.
  copy_maps(copy);
}
