// timer.cpp
// Ben Gamari
// Januray 2009

#include "timer.h"
#include <cstdio>
#include <time.h>



void timer::start(const std::string &m)
{
  msg = m;
    if(verbose) printf("TM_BEGIN: %-40s ...\n", msg.c_str());
  gettimeofday(&tm, NULL);
};

void timer::stop(const std::string &m)
{
  timeval tm2;
  gettimeofday(&tm2, NULL);
  time = tm2.tv_sec - tm.tv_sec;
  time += (tm2.tv_usec - tm.tv_usec)/1000000.0;
    if(verbose) printf("TM_END:   %-20s %-20s... done in %.3e sec.\n",
      msg.c_str(), m.c_str(), time);
  fflush(stdout); 
}

void timer::rstart()
{
  gettimeofday(&tm, NULL);
}

void timer::rstop()
{ 
  timeval tm2;
  gettimeofday(&tm2, NULL);
  time = tm2.tv_sec - tm.tv_sec;
  time += (tm2.tv_usec - tm.tv_usec)/1000000.0;
}

double timer::get_time()
{
  return this->time;
}

void time_stamp(const char *msg)
{
  time_t time_stamp;

  {
    time(&time_stamp);
    fprintf(stderr, "%s: %s\n", msg, ctime(&time_stamp));
    fflush(stderr);
  }
}

void stopwatch::report()
{
    printf("SW_END:   %-20s ... done in %.3e sec.\n",
      tag.c_str(), time);
  fflush(stdout); 
}
