#include "stdafx.h"

#include "Misc.h"
#include "kernel.h"



void			SignalHandler(int s)
{
	printf("Caught signal %d\n", s);
	signaled = 1;
}
void			init_signal()
{
	typedef void(*SignalHandlerPointer)(int);

	//SignalHandlerPointer previousHandler;
	signal(SIGINT, SignalHandler);
	signal(SIGTERM, SignalHandler);
	//signal(SIGABRT, SignalHandler);
}


unsigned int	next_power_of_two(unsigned int x)
{
	unsigned int ret = 1;
	while (ret < x)
	{
		ret *= 2;
	}
	return ret;
}
unsigned int	next_multiple_of(unsigned int x, unsigned int n)
{
	unsigned int ret = n;
	while (ret < x)
	{
		ret += n;
	}
	return ret;
}
bool			should_print_step(int i, int m)
{
	return true;
	
	int a = m / 100;

	if (a > 0)
	{
		if (i % a == 0) return true;
	}
	else
	{
		return true;
	}

	return false;
}
bool			should_save_frame(int i)
{
	if ((i % 10) == 0) return true;
	return false;
}

sig_atomic_t signaled = 0;