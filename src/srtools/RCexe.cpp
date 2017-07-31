#include <stdio.h>
#include <iostream>
#include <fstream>
#include "timers/hr_time.h"
#include "triangulations/RobustCrust.h"

#include <stdint.h>
#include "util/shell_io.h"


using namespace std;

//char for the header
char algoname[256] = "Robust Crust";
char version[256] = "v0.2";
char author[256] = "Luigi Giaccari";
char user[256] = "WhoEver";

double* ptemp; //points coordinates as imported from file manager
int N;//number of points
vector<double> p;
int main(int argc, char* argv[])
{

    int flag;

    CStopWatch Timer;

    RCRUST Surface;//Power Crust class
    PrintHeader(algoname, version, author, user);


    cout << endl << "PROGRAM STARTED!" << endl;


    if (argc > 1) //read input when provided
    {
        cout << "Reading input!" << endl;
        ReadInputs(argc, argv);
    }
    else//otherwise asks for inputs/output
    {
        AskForInputs();
    }

    cout << "Importing points...";
    FileManager.Read_Points(&ptemp, &N, inputfile) ;

    p.resize(N * 3);
    Copy(&p[0], &ptemp[0], N * 3);
    Deallocate(&ptemp);


    cout << N << " Points imported" << endl;
    Timer.startTimer();
    if (N < 4)
    {
        Error("Not enough points to launch the algorithm");
    }

    cout << endl << "LAUNCHING MyRobustCrust" << endl;
    Surface.TriangulatePowerCrust(&p, N);
    /*Surface.TriangulateABPA(&p, N);*/
    int nt = Surface.t.size();

    Timer.stopTimer();
    cout << "Total Time: " << Timer.getElapsedTime() / 1000 << " s" << endl;


    cout << endl << "OUTPUT!" << endl;
    cout << "Writing stl file...";
    FileManager.Write_stl(&p[0], &Surface.t[0].p1, nt, outputfile);
    cout << "Done" << endl;


    //if(ShowSurface){
    //     cout<<endl<<"PROGRAM ENDED!"<<endl<<endl;
    //    VIEWER::StartViewer(p,N,&Surface.t[0].p1,nt);
    //    }
    //memory deallocation
    Surface.FreeMemory();


    cout << endl << "PROGRAM ENDED!" << endl << endl;
    cout << endl << "UNMESHABLE IS NOTHING!" << endl << endl;

    ExitProgram(0);
}
