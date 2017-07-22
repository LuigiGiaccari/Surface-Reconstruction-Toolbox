#include <stdio.h>
#include <iostream>
#include <fstream>
#include "timers/hr_time.h"
#include "triangulations/QuickDel2D.h"
#include "triangulations/QuickDel2D.cpp"
#include "util/filemanager.h"
#include <stdint.h>
#include "util/shell_io.h"
using namespace std;

void QuickDelWrapper(double* p, int N, int* t, int* nt) //converts input and output for QuickDel and triangulates points
{
    int i, c;
    double* pin = NULL; //input points for delaunay triangulation

    QUICKDEL2D T;//quickdel class
    Allocate(&pin, (N + 3) * 2); //N 2d points;+ 3 for service triangulation


    //copy only x and y coordinate inside pin
    c = 0;
    for (i = 0; i < N; i++)
    {
        pin[i * 2] = p[i * 3]; //3D to 2D points
        pin[i * 2 + 1] = p[i * 3 + 1];

        //jump z
    }

    /*
    double xmin,xmax,ymin,ymax,Maximum;
     xmin=Min(&pin[0], N, 2);
     ymin=Min(&pin[1], N, 2); ;
     xmax=Max(&pin[0], N, 2);
     ymax=Max(&pin[1], N, 2);



     //move to origin
     for(i=0;i<N;i++) {
         pin[i*2]-=xmin;
         pin[i*2+1]-=ymin;
     }

     //normalize to 1
     ((xmax-xmin)>(ymax-ymin))?  Maximum=xmax-xmin:Maximum=ymax-ymin;//get maximum value


     for(i=0;i<N*2;i++) {
         pin[i]/=Maximum;
     }
    */
    //launching QuickDel
    T.Options_UseBrute = true;
    T.Options_UseHCPO = true;
    T.Options_Reindex = true;
    T.Options_CheckTriangulation = false;
    T.Options_PrintStats = false;
    if (T.Triangulate(pin, N) != 0) //compute delaunay tiangulation of pin
    {
        Error("Error during Delaunay Triangulation");
    }

    //memory deallocation
    Deallocate(&pin);

    //copy output triagngulation
    c = 0;
    for (i = 0; i < T.ct; i++)
    {
        if (T.t[i].p1 < N && T.t[i].p2 < N && T.t[i].p3 < N) //don't copy the service triagnulation
        {
            t[c] = T.t[i].p1;
            c++;
            t[c] = T.t[i].p2;
            c++;
            t[c] = T.t[i].p3;
            c++;
        }
    }
    *nt = c / 3; //number of triangles


    //reindex triangulation in case UseHPCPO is ON
    if (T.Options_UseHCPO)
    {
        for (i = 0; i < *nt * 3; i++)
        {
            t[i] = T.HCPO.index[t[i]];
        }
    }

    //T.~QUICKDEL2D();//destroy triagnulation object(non c'è n'è bisogno il costruttore è chiamamto automaticamente
}


int main(int argc, char* argv[])
{
    //char for the header
    char algoname[256] = "QuickTerrain";
    char version[256] = "v0.0";
    char author[256] = "Luigi Giaccari";
    char user[256] = "WhoEver";

    int flag;
    int N, nt; //numbers of points and triangles
    double* p = NULL; //list of doubles with input points
    int* t = NULL; //list of integers with triangles indexes
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
        AskForInputs();
    }

    cout << "Importing points...";
    FileManager.Read_Points(&p, &N, inputfile) ;
    Allocate(&t, N * 3 * 2); //allocate memory for triangles (N*2)numbers of triagnles 3 number os points for each triagnle

    cout << N << " Points imported" << endl;
    if (N < 4)
    {
        Error("Not enough points to launch the algorithm");
    }



    cout << endl << "LAUNCHING QUICKDEL" << endl;
    Timer.startTimer();
    QuickDelWrapper(p, N, t, &nt);
    Timer.stopTimer();
    cout << "Number of Triangles=" << nt << endl;
    cout << "Delaunay Triangulation Time: " << Timer.getElapsedTime() / 1000 << " s" << endl;


    cout << endl << "OUTPUT!" << endl;
    cout << "Writing stl file...";
    FileManager.Write_stl(p, t, nt, outputfile);
    cout << "Done" << endl;

    cout << endl << "PROGRAM ENDED!" << endl << endl;

    cout << endl << "UNMESHABLE IS NOTHING!" << endl << endl;
    //memory deallocation
    Deallocate(&p);
    Deallocate(&t);
    ExitProgram(0);
}
