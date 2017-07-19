#pragma once
#ifndef _SDS3D_MAP_h__
#define _SDS3D_MAP_h__



//SDS (3D version)
#include <math.h>
#include <stdlib.h>
#include <iostream>

using namespace std;

//The structure used for 3D points
struct Coord3D
{
	    double x;
		double y;
		double z;
};



class SDS3D_MAP
{
    
	

private:
  int *First;//Access table
  int *Next;//next point into the box
  double Minx,Miny,Minz;
  int nx,ny,nz;
  Coord3D* p;
  
  
  int* idbox;
  int* first;
  int Mapsize;
  
  
  
  void Reset();
   void PrintStatistics();
  void BuildMap();
  int GetFirst(int idbox);

  public:
      
      int Np;////number of stored points
    //int np;// current nunmber of points
    int Nb;//number of boxes
    double Empty;//Fraction of empty boxes
    double density;//(number of point)/(number of boxes) default=.5
    double step;//dimensioni delle boxes
    int* idStore;
	int npts;//numbers of points found during search
	
	

//funzioni membro
    SDS3D_MAP();//constructor
	~SDS3D_MAP();//destructor
    int BuildSDS(Coord3D* inputp, int Np);
    void SearchClosest(Coord3D *pk,int* idc,double* mindist);
    void SearchClosestExclusive(Coord3D *pk,int* idc,double* mindist,int sp);
	void SearchKClosest(Coord3D *pk,int* idc,double* mindist,int k);
      void SearchRadius( Coord3D* pk, double r);
     int EmptyBallTest(Coord3D* pk,double sqrdist);
	  void RemovePoint(int idp);
	  void SearchKClosestExclusive( Coord3D *pk, int* idc, double* distances, int k,int sp);
	  void SearchCuboid( double* Cuboid);
        void SearchRadiusExclusive( Coord3D* pk, double r,int sp);
        void GetPointsInRange( double* Cuboid);      
        void Deallocate();
        
};


#endif
