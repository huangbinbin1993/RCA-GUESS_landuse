#include <sys/time.h>
#include <sys/resource.h>

void csecond_(float *cputime)
{
        struct rusage ru;
        getrusage(RUSAGE_SELF, &ru);
 
        *cputime = ru.ru_utime.tv_sec + 1.0e-6 * ru.ru_utime.tv_usec +
                   ru.ru_stime.tv_sec + 1.0e-6 * ru.ru_stime.tv_usec ;
}
