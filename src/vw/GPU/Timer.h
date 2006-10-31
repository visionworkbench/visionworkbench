/****	Timer.h				****/
/****	Ian Saxton (2005)	****/

#ifndef Timer_H
#define Timer_H

#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

class Timer  {
protected:
// DATA
	timeval  startTime;
        timeval  endTime;
	long averageCount;
	double  averageValue;
	bool enabled;

public: 
// INLINE
	Timer() { Clear(); enabled = true; }
	void Clear() { averageCount = 0; averageValue = 0; }
	void Start() { if(enabled) gettimeofday(&startTime, NULL);  }
	void Stop() { if(enabled) gettimeofday(&endTime, NULL); }
	void SetEnabled(bool value) { enabled = value; }
// MEMBER
	void StopAndPrint(char* text);
	void AddToAverage();
	void Print(double scalar = 1.0);
	double ElapsedSeconds();

// VIRTUAL / OVERLOADED  
};


#endif

