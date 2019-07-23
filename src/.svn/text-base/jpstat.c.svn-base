/* 
Revision 1.1  2000/02/21 08:40:06  JussiMaki
a set of profiling statistics routines
*/


/* jpstat.c - Jussi.Maki@csc.fi a set of profiling statistics routines */
/* Created: 1999-06-07 */

#include <stdio.h>
#include <stdlib.h>
#include<string.h>
#include <time.h>
#include <sys/time.h>
#include <sys/times.h>
#include <sys/resource.h>

typedef struct jpstat {
  char opername[16];
  int oper_count;
  long long byte_count;
  double walltimesum;
  double cputimesum;
  double walltime_start;
  double cputime_start;
} jpstat_t;

#define JPSIZE 16

#ifdef CRAY
#define jpstat_init_ JPSTAT_INIT
#define jpstat_set_opername_ JPSTAT_SET_OPERNAME
#define jpstat_enable_ JPSTAT_ENABLE
#define jpstat_disable_ JPSTAT_DISABLE
#define jpstat_zero_ JPSTAT_ZERO
#define jpstat_start_ JPSTAT_START
#define jpstat_stop_ JPSTAT_STOP
#define jpstat_print_ JPSTAT_PRINT
#endif

static jpstat_t *jp;
static int jpsize=0;
static int jpuse=-1;
static int jpstat_enabled=1;

void jpstat_init_(int *count);
void jpstat_set_opername_(int *operid, char *opername);
void jpstat_start_(int *operid);
void jpstat_stop_(int *operid, int *bytes);
static double walltime();
static double cputime();

void jpstat_init_(int *count)
{
  jp = (jpstat_t *) calloc(sizeof(jpstat_t), *count);
  jpuse=0;
  jpsize=(*count);
}

void jpstat_set_opername_(int *operid, char *opername)
{
  if ((*operid)<=jpsize) {
    strncpy(jp[*operid].opername, opername,15);
    jp[*operid].opername[15]='\0';
  }
}

void jpstat_zero_()
{
  int i;
  for (i=0; i<jpsize; i++) {
    jp[i].oper_count  = 0;
    jp[i].byte_count  = 0;
    jp[i].walltimesum = 0.0;
    jp[i].cputimesum  = 0.0;
    jp[i].walltime_start = 0.0;
    jp[i].cputime_start  = 0.0;
  }
}

void jpstat_enable_()
{
  jpstat_enabled=1;
}
void jpstat_disable_()
{
  jpstat_enabled=0;
}

/* start the counting */
void jpstat_start_(int *operid)
{
  if (!jpstat_enabled) return;
  jp[*operid].walltime_start = walltime();
  jp[*operid].cputime_start  = cputime();
}

/* stop the current counting */
void jpstat_stop_(int *operid, int *bytes)
{
  if (!jpstat_enabled) return;
  jp[(*operid)].walltimesum += walltime() - jp[(*operid)].walltime_start;
  jp[(*operid)].cputimesum  += cputime()  - jp[(*operid)].cputime_start;
  jp[(*operid)].byte_count  += (*bytes);
  jp[(*operid)].oper_count ++;
}

/* Print out the current statistics */
void jpstat_print_()
{
  int i;
  
  struct timeval now;
  struct tm tm;

  if (!jpstat_enabled) return;
  fflush(stdout);
  gettimeofday(&now, NULL);
  printf("JPSTAT PROFILING STATISTICS TIMESTAMP %d.%06d\n",(int)now.tv_sec,(int)now.tv_usec);
  printf("OPERNAME          COUNT    KBYTES WALLTIME CPUTIME KBYTES/SEC KBYTES/OPER\n");
  for (i=0; i<jpsize; i++) {
    if (jp[i].oper_count) {
      if (jp[i].walltimesum==0.0) jp[i].walltimesum=0.0001;
      if (jp[i].cputimesum==0.0) jp[i].cputimesum=0.0001;
      if (jp[i].oper_count==0.0) jp[i].oper_count=1;
      printf("%-16s %6d %9.2f %7.2f %7.2f %7.2f %7.2f\n", 
	     jp[i].opername, jp[i].oper_count, jp[i].byte_count/1024.0, 
	     jp[i].walltimesum,jp[i].cputimesum, 
	     jp[i].byte_count/jp[i].walltimesum/1024.0,
	     jp[i].byte_count/(double)jp[i].oper_count/1024.0);
    }
  }
  fflush(stdout);
}

static double walltime()
{
  static struct timeval tvsave;
  struct timeval tv;
  static int initted=0;
  if (!initted) {
     initted=1;
     gettimeofday(&tvsave, NULL);
  }
  gettimeofday(&tv,NULL);
  return ((double)tv.tv_sec-1-tvsave.tv_sec + 1.0e-6*(1000000+tv.tv_usec-tvsave.tv_usec));
}

#if defined(RUSAGE_SELF)
static double cputime()
{
  struct rusage ru;
  getrusage (RUSAGE_SELF, &ru);
  return ((double)ru.ru_utime.tv_sec + ((double)ru.ru_utime.tv_usec * 1.0e-6));
}

#else
#include <sys/param.h>

static double cputime ()
{
  register long walltime;
  struct tms tm;
  walltime = times(&tm);
  return (double)tm.tms_utime / (double)CLK_TCK;
}
#endif
