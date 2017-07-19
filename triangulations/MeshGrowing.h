//Mesh Growing library

//TODO

//-Seed triagnle for bpa
//-Remove duplicate points during NNgraph construction
//-Detect duplicate triagnle bug
//-Detect inverse orientation bug


#pragma once
#ifndef _MESHGROWING_h__
#define _MESHGROWING_h__




#include <iostream>
#include <fstream>
#include "structs/Coord3D.h"
#include "structs/Triangle.h"
#include "SDS/SDS3D.h"//contains the coord3D struct
#include "SDS/SDS3D.cpp"//contains the coord3D struct
#include "util/SmallPriorityQueue.h"
#include "util/ArraysLib.h"
#include "triangulations/PointEdgeMap.h"
#include "triangulations/TrianglesFunctions.h"
#include <stdint.h>//for 16 bits integers


using namespace std;
//internal constants

#define K 3 //for the seed triangle search

#ifndef M_PI
#define M_PI  3.14159265358979323846f
#endif


//Structure for edges
struct Edge
{
    //points of the edge
    int p1;
    int p2;

    //neighbour triangles
    int t1;
    int t2;
};

//internal paramenters for algorithms
struct PARAMETERS
{
    double BPA_Max_BCAngle;
    double BPA_Min_BCAngle;
    double Max_Edge_Length;
    double BPA_FlatnessToll;
    double NN_Filter_Cutdist;
    int SCB_Dotp_PriorityLevel;// 5 per il fltness test
    int SCB_SR_PriorityLevel;//2 per il search radius
    double SCB_TooSharp;
    double BPA_TooSharp_Torus;
    double BPA_TooSharp_SR;
    bool EnablePreProcessing;
    bool EnablePostProcessing;
    bool BPA_ReverseNormal;

};

struct Torus3D
{
    Coord3D Center;
    Coord3D Axis;
    double d;//plane coefficient= Dot(axis,center)
    double R;
    double r;
};



class MESHGROWING
{
private:


//GLOBALS

//Enums
    enum DeallocateMode {Mem_All,Mem_Save_t,Mem_Save_t_e};
    enum ErrorMode {Er_Ok,Er_OutOfMemory,Er_CannotAllocate};

//Variables

    int* queue;//edges queue


//entities counters
    int countt;
    int counte;


    int MAXT;//maximum number or triangles
    int MAXE;//maximum numbers of edges







    int16_t *NE;//number of connected edges
    int16_t *NT;//number of connected points
    char*nbe;//number of boundary edges for each point

    int Ndel;//number of deleted points
    int Ntdel;//number of deleted triangles

//Point edge map
    PointEdgeMap EPMap; //EPmap

//Search data structure
    SDS3D SDS;//GLtree

//Priority queue
    SMALLPRIORITYQUEUE Pqueue;
    //point id map

//          P2
//        / | \
// 	T1  P3  |  P4(newpoint) T2
//        \ | /
//		    P1

//NOrmal of the new triangle cross(v(P4-P1),v(P2-P1))
// New triagnle is store das [P1 P4 P2]
//NOrmal of the generic triangle [P1,P2,P3] cross(v(P2-P1),v(P3-P1))
//
    int P1,P2,P4,T1,ide;//temporary values
    Coord3D v21,v41;//vectors of the new triangle
    double* FlatToll;//array for flatness tolerance in SCB

    PARAMETERS P;

//Private functions Prototypes


//MEMORY FUNCTIONS
    void Memory_Deallocate(DeallocateMode Mode=Mem_All);//memory deallocation
    void Memory_Allocate();//memory allocation

//UTILITY FUNCTIONS
    void MakePublic();//Makes public important private variables
    void InitParameters();//initialize internal parameters
    void PrintParameters(int mode);//prints internal parameters
    void EdgeQueue_Fill_With_Boundary(int sr_value);//fill the queue with boundary edges

//MESHING FUNCTIONS
    void SeedTriangle(bool BPAtest);//Gets the seed triangle (BPAtest dafault=false)
    void RunFront_BPA();//check te free edge queue to expand the front
    void GetCandidates_Torus(Torus3D* Torus);
    void SelectCandidate_BPA_Torus(Torus3D* Torus);
    bool CheckEdgeConformity(int idedge,int point);
    void UpdateFront(int idedge1,int idedge2);
    void UpdateFront_Seed(int i,int idc,int id3);
    void CircumCenter(Coord3D* p1,Coord3D* p2,Coord3D* p4,Coord3D* n,double *r,Coord3D* cc);
    void GetTriangle_BPA_Torus();
    void GetTriangle_SCB();
    void GetTriangle_BPA_SR() ;
    void RunFront_SCB();
    void SelectCandidate_BPA_SR();

//GEOMETRY FUNCTIONS
    inline void TNormTriangle(int idt, Coord3D* tnorm);
    void TNormPoints(Coord3D* p1,Coord3D* p2,Coord3D* p4,Coord3D* tnorm);
    bool BallCenter(Coord3D* cc,double* cr, Coord3D* tnorm,Coord3D* BC);
    double BCAngle(Coord3D* m,Coord3D* BC,Coord3D* tnorm);
    void GetSearchPoint(Coord3D* P1,Coord3D* P2,Coord3D* tnorm,Coord3D* sp,double* sr, double kr);//gets the search point for SCBMesher
    void CircumCenter_Fast(Coord3D* p1,Coord3D* p2,Coord3D* p4,Coord3D* n,double *r,Coord3D* cc);
    bool InsideTorus(Torus3D* Torus,Coord3D* point);
    void BuildTorus(Torus3D* Torus);
    bool CylinderTest(Coord3D cc,double r2,Coord3D tnorm,bool* up0);

//PRE-PROCESSING
    void PreP_NN_MeanDist(double*meandist);
    void PreP_NN_Filter(double cutdist);




//POST-PROCESSING
    void PostP_FillQueue(int value);//Fill queue with boundary edges
    int PostP_FindNMV(char* manifold);//finds not manifold vertex
    void PostP_HoleFiller();//fills holes
    void PostP_FillTriplets();//Fill Triplets
    void PostP_HealNMV();//Deletes triagnles related with NMV
    void PostP_Restore_t();////restore triagnle id connectivity due to deleted triagnles
    void PostP_Delete_NMV_BoundTriangle_idp(int idt,int idp);//Deletes a boundary triagnle related with NMV idp
    void PostP_Delete_NMV_Triangles_All(int idp);//deletes all triangles related with the NMV idp
    void PostP_Delete_NMV_OverBound_Triangle(int idt);//Deletes the over boundary triangles idt
    int PostP_UpdateEdge_DeletedTriangle(int p1,int p2,int idt);//Updates data structure after a triagnle has been deleted
    void PostP_PostProcessing();//the gataway routine that inlcudes postprocessing operations
public:

    //Points
    int N;//number of points
    Coord3D* p;//i punti della nuvola
    Triangle* t;//Triangles
    int nt;
    Edge* e;//oversized edge structure
//Nota si assume che e abbia come valori di defautl valori negativi...
//altrimenti si rende necessario un loop di inizilizzazione
    int ne;// edges number

    double R;//The Ball Radius


//public functions prototypes
    MESHGROWING();//constructor
    ~MESHGROWING();//destructor

    void BallPivoting(double R);//contains the strategy of tessellation
    void SCBMesher();
    void ImportPoints_Pointers(double* pointer,int inputN);


};


#endif
