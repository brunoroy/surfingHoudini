#include "PerfTimerBasicLinux.h"

#ifndef WIN32

#warning IN PERFTIMERBASICLINUX

namespace
{
	double convertTimespecToSeconds(const timespec& t)
	{
		return static_cast<double>(t.tv_sec) + static_cast<double>(t.tv_nsec)*0.000000001;
	}
}

// ----------------------------------------------------------------------------
// Constructor / Destructor
PerfTimerBasicLinux::PerfTimerBasicLinux()
: _elapsedTime(0.0),
  _isRunning(false)
{
	_startTime.tv_sec = 0;
	_startTime.tv_nsec = 0;
}

PerfTimerBasicLinux::~PerfTimerBasicLinux()
{
}

// ----------------------------------------------------------------------------
// Public functions
void PerfTimerBasicLinux::start()
{
	if (!_isRunning)
	{
		if (clock_gettime(CLOCK_MONOTONIC, &_startTime) !=0)
		{
			_startTime.tv_sec = 0;
			_startTime.tv_nsec = 0;
		}
		_isRunning = true;
	}
}

void PerfTimerBasicLinux::pause()
{
	if (_isRunning)
	{
		timespec stopTime;
		timespec deltaTime;

		if (clock_gettime(CLOCK_MONOTONIC, &stopTime) != 0)
		{
			deltaTime.tv_sec = 0;
			deltaTime.tv_nsec = 0;
		}
		else
		{
			deltaTime.tv_sec = stopTime.tv_sec - _startTime.tv_sec;
			deltaTime.tv_nsec = stopTime.tv_nsec - _startTime.tv_nsec;
		}

		_elapsedTime += convertTimespecToSeconds(deltaTime);

		_isRunning = false;
	}
}

void PerfTimerBasicLinux::reset()
{
	_elapsedTime = 0.0;
	if (clock_gettime(CLOCK_MONOTONIC, &_startTime) != 0)
	{
		_startTime.tv_sec = 0;
		_startTime.tv_nsec = 0;
	}
}

double PerfTimerBasicLinux::getElapsedTime()
{
	double totalTime = _elapsedTime;
	if (_isRunning)
	{
		timespec stopTime;
		timespec deltaTime;

		if (clock_gettime(CLOCK_MONOTONIC, &stopTime) != 0)
		{
			deltaTime.tv_sec = 0;
			deltaTime.tv_nsec = 0;
		}
		else
		{
			deltaTime.tv_sec = stopTime.tv_sec - _startTime.tv_sec;
			deltaTime.tv_nsec = stopTime.tv_nsec - _startTime.tv_nsec;
		}

		totalTime += convertTimespecToSeconds(deltaTime);
	}

	return totalTime;
}

#endif // !WIN32
