#ifndef _hr_timer_h

#define _hr_timer_h

#include <sys/time.h>







class CStopWatch {
   timeval t1, t2;
 
public:
    CStopWatch() ;
    void startTimer( ) ;
    void stopTimer( ) ;
    double getElapsedTime() ;
};




CStopWatch::CStopWatch(){

//init t1,t2
   gettimeofday(&t1, NULL);
 gettimeofday(&t2, NULL);
}

void CStopWatch::startTimer( ) {
      // start timer
    gettimeofday(&t1, NULL);
}

void CStopWatch::stopTimer( ) {
     // stop timer
    gettimeofday(&t2, NULL);
}

double CStopWatch::getElapsedTime() {
    	timeval res;
	timersub(&(timer.stop),&(timer.start),&res);
	return res.tv_sec + res.tv_usec/1000000.0; // 10^6 uSec per second


}


#endif
