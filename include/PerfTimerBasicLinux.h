#ifndef PERFTIMERBASICLINUX_H
#define PERFTIMERBASICLINUX_H

#ifndef WIN32

#include <time.h>

class PerfTimerBasicLinux
{
public:
	PerfTimerBasicLinux();
	virtual ~PerfTimerBasicLinux();

	virtual void start();
	virtual void pause();
	virtual void reset();

	virtual double getElapsedTime();	// 1.0 == 1 second

private:
	timespec	_startTime;
	double		_elapsedTime;
	bool		_isRunning;
};

#endif	// !WIN32
#endif	// PERFTIMERBASICLINUX_H
