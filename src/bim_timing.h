#ifndef HAVE_BIM_TIMING_H
#define HAVE_BIM_TIMING_H 1

#ifdef BIM_TIMING

#include <ctime>
static clock_t c_start, c_diff;
static double c_msec;

struct event_counter
{
  double cumtime;
  int    count;
};
  
static std::map<std::string,event_counter> timing_report;

#define tic() c_start = clock ();

#define toc(X)                                                \
  c_diff = clock () - c_start;                                \
  c_msec = (double)c_diff * 1000 / (double)CLOCKS_PER_SEC;    \
  timing_report[X].cumtime += c_msec;                         \
  timing_report[X].count++;                                   \
  std::cout << X << std::endl << "Elapsed time : " << c_msec  \
    << "ms" << std::endl;                                     


#define print_timing_report()                                           \
  double time_counter = 0;                                              \
  std::cout << "Timing Report:" << std::endl;                           \
  for (std::map<std::string,event_counter>::iterator ii =               \
         timing_report.begin ();                                        \
       ii != timing_report.end (); ++ii)                                \
    time_counter += (*ii).second.cumtime;                               \
  for (std::map<std::string,event_counter>::iterator ii =               \
         timing_report.begin ();                                        \
       ii != timing_report.end (); ++ii)                                \
    std::cout << "Event: "                                              \
              << (*ii).first                                            \
              << ", total hits: "                                       \
              << (*ii).second.count                                     \
              << ", total time: "                                       \
              << (*ii).second.cumtime / 1.0e3                           \
              << " s. ("                                                \
              << 100 * (*ii).second.cumtime / time_counter              \
              << "p.c.)" << std::endl;                                             

#else 

#define tic()        \
  // timing disabled 
  
#define toc(X)       \
  // timing disabled 

#define print_timing_report() \
  // timing disabled 
  
#endif

#endif
