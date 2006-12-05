/*
 *  Timer.cpp
 *  iFramework
 *
 *  Created by ian on 8/16/05.
 *  Copyright 2005 __MyCompanyName__. All rights reserved.
 *
 */

#include "Timer.h"



void Timer::StopAndPrint(char* text) {
  if(!enabled) return;
  Stop();
  if(text)
    printf("%s    ", text);
  Print();
}


void Timer::Print(double scalar) {
  if(!enabled) return;

	long long int deltaTimeNanos;
	double endTimeFloat;
	double startTimeFloat;
	if(averageCount==0) {
	  endTimeFloat = (endTime.tv_sec + endTime.tv_usec/1000000.0);
	  startTimeFloat = (startTime.tv_sec + startTime.tv_usec/1000000.0);
	  deltaTimeNanos = (long long int) (1000000000 * scalar * (endTimeFloat - startTimeFloat));
	}
	else {
		printf("Average Time, %i runs: ", averageCount);
		deltaTimeNanos = (long long int) (1000000000.0 * (scalar * averageValue));
	}
	int seconds = (int) (deltaTimeNanos / 1000000000);
	int ms = (int) (((long long)  deltaTimeNanos % 1000000000) / 1000000);
	int us = (int) (((long long) deltaTimeNanos % 1000000) / 1000);
	int ns = (int) ((long long)  deltaTimeNanos % 1000);
	if (scalar != 1.0) 
		printf("(Scaled by %.6f) ", scalar);
	printf("%i.%03i,%03i Seconds\n", seconds, ms, us);
}

double Timer::ElapsedSeconds() {
  double end = (endTime.tv_sec + endTime.tv_usec/1000000.0);
  double start = (startTime.tv_sec + startTime.tv_usec/1000000.0);
  return (end - start);
}


void Timer::AddToAverage() {
  //	averageValue = ((averageValue*averageCount) + (endTime - startTime)) / (double) (averageCount + 1);
  //	averageCount++;
}
