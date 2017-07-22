#include <stdio.h>
#include <iostream>
#include <fstream>

#include "timers/hr_time.h"

#include "triangulations/MeshGrowing.h"
#include "triangulations/MeshGrowing.cpp"
#include "util/filemanager.h"
#include "util/shell_io.h"
#include <string.h>
#include <stdint.h>


using namespace std;

//char for the header
char algoname[256] = "Ball Pivoting Algorithm";
char version[256] = "v0.1";
char author[256] = "Luigi Giaccari";
char user[256] = "WhoEver";

double* p;
int N;

int main(int argc, char* argv[])
{

    int flag;
    MESHGROWING Surface;//MEshgrowing class
    CStopWatch Timer;


    PrintHeader(algoname, version, author, user);


    cout << endl << "PROGRAM STARTED!" << endl;


    if (argc > 1) //read input when provided
    {
        cout << "Reading input!" << endl;
        ReadInputs(argc, argv);
    }
    else//otherwise asks for inputs/output
    {
        AskForInputs(1);//BPA requires mode one
    }

    cout << "Importing points...";
    FileManager.Read_Points(&p, &N, inputfile) ;


    cout << N << " Points imported" << endl;
    Timer.startTimer();
    if (N < 4)
    {
        Error("Not enough points to launch the algorithm");
    }


    cout << endl << "LAUNCHING BPA" << endl;
    Surface.ImportPoints_Pointers(p, N);
    Surface.BallPivoting(R);
    // Surface.SCBMesher();


    Timer.stopTimer();
    cout << "Total Time: " << Timer.getElapsedTime() / 1000 << " s" << endl;


    cout << endl << "OUTPUT!" << endl;
    cout << "Writing stl file...";
    FileManager.Write_stl(p, &Surface.t[0].p1, Surface.nt, outputfile);
    cout << "Done" << endl;


    //memory deallocation
    Deallocate(&p);

    cout << endl << "PROGRAM ENDED!" << endl << endl;
    cout << endl << "UNMESHABLE IS NOTHING!" << endl << endl;
    ExitProgram(0);
}
