// timer.h
// Andrei Alexandru
// September 2009

#ifndef _TIMER_H
#define _TIMER_H

#include <sys/time.h>
#include <string>


struct timer
{
  bool verbose;
  std::string msg;
  timeval tm;
  double time;

  timer(bool verb = true): verbose(verb) {}

  void start(const std::string &m);
  void stop(const std::string &m = "");
  double get_time();
  void rstart();
  void rstop();
};

void time_stamp(const char*);

struct stopwatch
{
  std::string tag;
  timer tm;
  double time;

  stopwatch(const std::string &m = ""): tag(m), tm(false), time(0) {};
  double get_time() { return time; }

  void start() { tm.start(""); }
  void stop()  { tm.stop(); time += tm.get_time(); }
  void reset() { time=0; }
  void reset(const std::string &m) { tag = m; time=0; }

  void report();
};


#endif


