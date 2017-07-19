#ifndef _QUICKDEL2D_H_
#define _QUICKDEL2D_H_

#define STATS//enables stats printing
#define ROBUST//enables exact arithmetic

//STRUCTURES for points and triangles

struct TriangleDel
{
    int p1;//points index
    int p2;
    int p3;
    int t1;//neigh triangles
    int t2;
    int t3;
};

#ifdef STATS
struct Statistics
{
    long int incircle_call;
    long int orient_call;
    long int incircle_support_call;//indicates the number of times a support edge has been checked for flip
    long int flips;
};

extern Statistics Stats;

#endif



#ifdef STATS
Statistics Stats;//structure for statistics
#endif


#include <cstdio>
#include "util/ArraysLib.h"//some array utility
#include "triangulations/predicates.h"//robust predicates
#include "util/ExceptionHandler.h"//defines exceptions
#include "sortlib/HCPO2D.h"//to sort points in HPCO order
#include "math.h"
#include <stdlib.h>
#include <iostream>
#include <cstring>
#include <fstream>



using namespace std;




class QUICKDEL2D
{
    private:

        //variables

        double xmin, xmax, ymin, ymax; //bounding box
        bool* toflip;//states if a triangle need to be flipped
        int* flipstack;//stack for triagnles to be flipped
        int  stack_iter;//iterator for flipstack
        int first_tria;//first  triangle before flip_them_all


        //functions

        //mesh utils
        void PreProcess(double* inputp, int inputnp); //preprocessing operations before triangulation
        void ComputeBoundingBox();//preprocess points before triangulation
        void BuildServiceTriangulation();//build service triangulation
        void RemoveServiceTriangulation();//removes service triangulation
        int LocatePoint(int* cT, int i); //locate the point inside the simplex cT=(current traingle)
        //mesh kernels
        int Triangulate_Brute(double* inputp, int inputnp); //brute incremental delaunay triangulation
        int Triangulate_Global(double* inputp, int inputnp); //incremental delaunay triagnulation  points addiction
        void AddPoint(int idt, int i, bool flip = true); //add the point i insidie the simplex idt. If bool=true attemps to flip
        void AddPointOnEdge(int idt, int i, int idedge); //Adds point lying on edge
        void Flip(int T1);//Flip the triangle T1 without incircle test
        void CheckForFlip(int T1);//Check the triagnle T1 for flip (flip is performed along the first edge)
        void FlipThemAll();//Flips all new triangles and those in the flip_stack
        void CheckForFlip_Global(int T1);//just like checkforlfip but for global delaunay
        //debug
        bool IsInside(int T1, int T2); //checkes if the apex of T2 is inside th circumcircle of T1
        //memory
        void FreeMemory();//deallocate


    public:
        //variables
        Coord2D* p;//points
        TriangleDel* t;//triangles
        int np;//number of points
        int ct;//counter of triangles (final numbers of triagnle=ct+1);
        HCPO2D HCPO;//points sorter
        //Options
        bool Options_OverWritePoints;//authorize to change  points array
        bool Options_UseHCPO;//enables the HPCO sorting
        bool Options_UseBrute;//enables brute triangulation
        bool Options_Reindex;//if true preserve indexes of HPCO sort
        bool Options_CheckTriangulation;//checks triagulation at the end
        bool Options_PrintStats;//print statistics if STATS define is enabled
#ifdef STATS
        void PrintStatistics();
#endif

        QUICKDEL2D();//constructor
        ~QUICKDEL2D();//destructor

        int Triangulate(double* inputp, int inputnp); //Calls a delaunay triangulation

        void PrintOutput_mfile();//print a matlab script to visualize traigulation
        bool IsDelaunay();//checks is triagnulations is truly delaunay


};

#endif



