
#ifdef __GNUC__
#pragma implementation
#endif

#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif

#include <iostream.h>
#include <iomanip.h>

#include <util/group/pregtime.h>

ParallelRegionTimer::ParallelRegionTimer(const RefMessageGrp&msg,
                                         const char *topname,
                                         int cpu_time, int wall_time):
  RegionTimer(topname, cpu_time, wall_time),
  msg_(msg)
{
}

ParallelRegionTimer::~ParallelRegionTimer()
{
}

void
ParallelRegionTimer::print(ostream &o)
{
  int i,j;

  if (msg_->n() == 1) {
      RegionTimer::print(o);
      return;
    }

  update_top();

  int n = nregion();

  int minn = n;
  int maxn = n;
  msg_->max(maxn);
  msg_->min(minn);

  if (maxn != minn) {
      cerr << "ParallelRegionTimer::print: differing number of regions"
           << endl;
      abort();
    }

  double *min_cpu_time = 0;
  double *min_wall_time = 0;
  double *max_cpu_time = 0;
  double *max_wall_time = 0;
  double *avg_cpu_time = 0;
  double *avg_wall_time = 0;
  if (cpu_time_) {
      min_cpu_time = new double[n];
      get_cpu_times(min_cpu_time);
      max_cpu_time = new double[n];
      get_cpu_times(max_cpu_time);
      avg_cpu_time = new double[n];
      get_cpu_times(avg_cpu_time);
      msg_->max(max_cpu_time,n);
      msg_->min(min_cpu_time,n);
      msg_->sum(avg_cpu_time,n);
      for (i=0; i<n; i++) {
          avg_cpu_time[i] /= msg_->n();
        }
    }
  if (wall_time_) {
      min_wall_time = new double[n];
      get_wall_times(min_wall_time);
      max_wall_time = new double[n];
      get_wall_times(max_wall_time);
      avg_wall_time = new double[n];
      get_wall_times(avg_wall_time);
      msg_->max(max_wall_time,n);
      msg_->min(min_wall_time,n);
      msg_->sum(avg_wall_time,n);
      for (i=0; i<n; i++) {
          avg_wall_time[i] /= msg_->n();
        }
    }

  if (msg_->me() == 0) {
      const char **names = new const char*[n];
      get_region_names(names);
      int *depth = new int[n];
      get_depth(depth);

      int maxwidth = 0;
      double maxtime = 0.0;
      for (i=0; i<n; i++) {
          int width = strlen(names[i]) + 2 * depth[i] + 2;
          if (width > maxwidth) maxwidth = width;
          if (cpu_time_ && max_cpu_time[i] > maxtime)
              maxtime = max_cpu_time[i];
          if (wall_time_ && max_wall_time[i] > maxtime)
              maxtime = max_wall_time[i];
        }

      int maxtimewidth = 4;
      while (maxtime > 1.0) { maxtime/=10.0; maxtimewidth++; }

      o.setf(ios::right);

      for (i=0; i<maxwidth; i++) o << " ";
      if (cpu_time_) {
          o << setw(maxtimewidth+1) << " ";
          o << setw(maxtimewidth+1) << " CPU";
          o << setw(maxtimewidth+1) << " ";
        }
      if (wall_time_) {
          o << setw(maxtimewidth+1) << " ";
          o << setw(maxtimewidth+1) << " Wall";
          o << setw(maxtimewidth+1) << " ";
        }
      o << endl;

      for (i=0; i<maxwidth; i++) o << " ";
      if (cpu_time_) {
          o << setw(maxtimewidth+1) << " min";
          o << setw(maxtimewidth+1) << " max";
          o << setw(maxtimewidth+1) << " avg";
        }
      if (wall_time_) {
          o << setw(maxtimewidth+1) << " min";
          o << setw(maxtimewidth+1) << " max";
          o << setw(maxtimewidth+1) << " avg";
        }
      o << endl;

      o.setf(ios::fixed);
      o.precision(2);
      for (i=0; i<n; i++) {
          int width = strlen(names[i]) + 2 * depth[i] + 2;
          for (j=0; j<depth[i]; j++) o << "  ";
          o << names[i] << ": ";
          for (j=width; j<maxwidth; j++) o << " ";
          if (cpu_time_) {
              o << " " << setw(maxtimewidth) << min_cpu_time[i];
              o << " " << setw(maxtimewidth) << max_cpu_time[i];
              o << " " << setw(maxtimewidth) << avg_cpu_time[i];
            }                    
          if (wall_time_) {
              o << " " << setw(maxtimewidth) << min_wall_time[i];
              o << " " << setw(maxtimewidth) << max_wall_time[i];
              o << " " << setw(maxtimewidth) << avg_wall_time[i];
            }
          o << endl;
        }

      delete[] names;
      delete[] depth;
    }

  delete[] min_cpu_time;
  delete[] max_cpu_time;
  delete[] avg_cpu_time;
  delete[] min_wall_time;
  delete[] max_wall_time;
  delete[] avg_wall_time;
}
