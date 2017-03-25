#ifndef NBODY_MISC_H
#define NBODY_MISC_H

#include <ctime>

#include "decl.h"
#include <signal.h>

class TimeMeasurement
{
public:
	clock_t t;
	double dif;

	void start()
	{
		t = std::clock();
	}

	void stop()
	{

		dif = (double)(std::clock() - t)/CLOCKS_PER_SEC;
	}

};

unsigned int	next_power_of_two(unsigned int x);

bool			should_print_step(int i, int m);
bool			should_save_frame(int i);



extern sig_atomic_t signaled;

void	SignalHandler(int s);
void	init_signal();

#endif

