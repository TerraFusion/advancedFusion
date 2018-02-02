

#include <iostream>
#include <stdio.h>
#include <sys/time.h>
#include "AF_debug.h"


struct timeval tvdbg1, tvdbg2, tvdbg_sub, tvdbg_add;

static void PrintTimeval (struct timeval *tv)
{
	long milliseconds;

	/* Compute milliseconds from microseconds.  */
	milliseconds = (*tv).tv_usec / 1000;
	/* Print the formatted time, in seconds, followed by a decimal point
	   and the milliseconds.	*/
	printf ("%ld.%03ld sec\n", (*tv).tv_sec, milliseconds);
}


void StartElapseTime()
{
	tvdbg_add = (struct timeval){0};
	gettimeofday (&tvdbg1, NULL);
}

static void MeasureElapseTime()
{
	gettimeofday (&tvdbg2, NULL);
	timersub(&tvdbg2, &tvdbg1, &tvdbg_sub);
	timeradd(&tvdbg_sub, &tvdbg_add, &tvdbg_add);
}

void StopElapseTimeAndShow(std::string msg)
{
	MeasureElapseTime();
	std::cout << msg << "  Elapsed time: ";

	PrintTimeval(&tvdbg_add);
	fflush(stdout);	
}
