#ifndef _hr_timer_h

#define _hr_timer_h


#ifdef _WIN32 || _WIN64 //windows version
#include <windows.h>//richiede librerie sdk



//High Resolution Timer Class
typedef struct
{
    LARGE_INTEGER start;
    LARGE_INTEGER stop;
} stopWatch;

class CStopWatch
{

    private:
        stopWatch timer;
        LARGE_INTEGER frequency;
        double LIToSecs( LARGE_INTEGER& L) ;
    public:
        CStopWatch() ;
        void startTimer( ) ;
        void stopTimer( ) ;
        double getElapsedTime() ;
};
double CStopWatch::LIToSecs( LARGE_INTEGER& L)
{
    return ((double)L.QuadPart / (double)frequency.QuadPart) ;
}

CStopWatch::CStopWatch()
{
    timer.start.QuadPart = 0;
    timer.stop.QuadPart = 0;
    QueryPerformanceFrequency( &frequency ) ;
}

void CStopWatch::startTimer( )
{
    QueryPerformanceCounter(&timer.start) ;
}

void CStopWatch::stopTimer( )
{
    QueryPerformanceCounter(&timer.stop) ;
}

double CStopWatch::getElapsedTime()
{
    LARGE_INTEGER time;
    time.QuadPart = 1000 * (timer.stop.QuadPart - timer.start.QuadPart); //Elapsed Time in ms
    return LIToSecs( time) ;
}


#else //linux version
#include <sys/time.h>



class CStopWatch
{
        timeval t1, t2;

    public:
        CStopWatch() ;
        void startTimer( ) ;
        void stopTimer( ) ;
        double getElapsedTime() ;
};




CStopWatch::CStopWatch()
{

    //init t1,t2
    gettimeofday(&t1, NULL);
    gettimeofday(&t2, NULL);
}

void CStopWatch::startTimer( )
{
    // start timer
    gettimeofday(&t1, NULL);
}

void CStopWatch::stopTimer( )
{
    // stop timer
    gettimeofday(&t2, NULL);
}

double CStopWatch::getElapsedTime()
{
    return  (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
}
#endif


#endif


