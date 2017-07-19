#ifndef _ExceptionHandler_h
#define _ExceptionHandler_h

#include <stdio.h>
#include <stdlib.h>

//globals
bool PAUSE_BEFORE_EXIT=true;//pause before exit?

//PROTOTYPES
void Error(const char* error,int flag);//display an error message and quit
void Warning(const char* warn);//display a warnign message
void ExitProgram(int exitflag=0);//exit program

void ExitProgram(int exitflag)//exit program
{
    //if pause is set to true it will pause before exiting
#ifdef WIN32
    if(PAUSE_BEFORE_EXIT)system("pause");

#else#linux "PAUSE"function
    if(PAUSE_BEFORE_EXIT)
    {
        printf("Press 'Enter' to exit the program");
        while (getchar() != '\n');
    }

#endif
    exit(exitflag);
}

//Utility Function
void Error(const char* error,int flag=1)//fatal error
{
    printf("FATAL ERROR:");
    printf(error);
    printf("\n");
    ExitProgram(flag);
}
void Warning(const char* warn)//warning
{
    printf("WARNING:");
    printf(warn);
    printf("\n");
}

#endif
