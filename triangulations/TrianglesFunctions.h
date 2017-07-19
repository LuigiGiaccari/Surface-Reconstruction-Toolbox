#ifndef __TrianglesFunctions_h__
#define __TrianglesFunctions_h__


#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <iostream>
#include "structs/Coord3D.h"
#include "structs/Triangle.h"

using namespace std;

//FUNTIONS DECLARATION GOES HERE
void GenerateRandom3D(Coord3D* p,int n) ;//Generate 3 random points
void Cramer(double *A,Coord3D *b,Coord3D *x);//Cramer 3x3 solver
 void CC(Coord3D* p1,Coord3D* p2,Coord3D* p3,Coord3D* cc,double* r);//Circumcenter of a 3D triangle
 double TriangleAngle(Coord3D* P1,Coord3D* P2,Coord3D* P3,Coord3D* P4);//finds the cosine of the angle between the 2 traingles formed by the 4 points
 void SearchPointExt(Coord3D* P1,Coord3D* P2,Coord3D* P3,Coord3D* sp,double* sr);//search point extended version (UNUSED)
 void SearchPoint(Coord3D* P1,Coord3D* P2,Coord3D* P3,Coord3D* sp,double* sr, double kr);//gets the search point for SCBMesher
inline void FastCramer(double* v1, double* v2, double* v3, double* b, double* sol);//Solve a linear system. v1 v2 v3 are the system rows, b is the coefficent vector and sol is the computed solution
 void TNorm(Coord3D* p1,Coord3D* p2,Coord3D* p4,Coord3D* tnorm);

//MACROS GOES HERE


//Inverte un vettore/punto
#define Reverse(p1)	 \
	p1 x=-p1 x;\
	p1 y=-p1 y;\
	p1 z=-p1 z;


//Prodotto vettoriale
#define CrossProduct(v1,v2,vect) \
	vect x=v1 y*v2 z-v1 z*v2 y;\
    vect y=v1 z*v2 x-v1 x*v2 z;\
	vect z=v1 x*v2 y-v1 y*v2 x;

//Trova il vettore dato da due punti, in pratica fa la differenza
#define DiffPoints(p1,p2,vect) \
	vect x=p1 x-p2 x;\
    vect y=p1 y-p2 y;\
	vect z=p1 z-p2 z;

//Prodotto scalare
#define DotProduct(p1,p2,DotP) \
	DotP= p1 x*p2 x + p1 y*p2 y+p1 z*p2 z;

//Distanza al quadrato
#define SquaredDistance(p1,p2,dist) \
	dist= (p1 x-p2 x)*(p1 x-p2 x)+(p1 y-p2 y)*(p1 y-p2 y)+(p1 z-p2 z)*(p1 z-p2 z);

//Distanza fra due punti
#define Distance(p1,p2,dist) \
	dist=sqrt( (p1 x-p2 x)*(p1 x-p2 x)+(p1 y-p2 y)*(p1 y-p2 y)+(p1 z-p2 z)*(p1 z-p2 z));

//Punto medio fra 2 punti
#define MidPoint(p1,p2,pm) \
	pm x= (p1 x+ p2 x)/2;\
	pm y= (p1 y+ p2 y)/2;\
	pm z= (p1 z+ p2 z)/2;

//Baricentro di un triangolo
#define Centroid(p1,p2,p4,c) \
	c x= (p1 x+ p2 x+ p3 x)/3;\
	c y= (p1 y+ p2 y+ p3 y)/3;\
	c z= (p1 z+ p2 z+ p3 z)/3;


//normalizza un coord3D data la sua lunghezza
#define Normalize(p1,leng) \
	 leng=sqrt((p1 x * p1 x)+(p1 y * p1 y)+(p1 z * p1 z)) ; \
	p1 x= p1 x/leng;\
	p1 y= p1 y/leng;\
	p1 z= p1 z/leng;

//trova il terzo punto dato un triangolo ed un edge
#define Setdiff(t,e,idp)\
     if ((t p1!=e p1) && (t p1!=e p2) )  { idp=t p1;}\
	 else if ((t p2!=e p1) && (t p2!=e p2) )  { idp=t p2;}\
	 else if ((t p3!=e p1) && (t p3!=e p2) )  { idp=t p3;}

//length of a vector
#define Length(v,leng)\
         leng=sqrt((v x * v x)+(v y * v y)+(v z * v z)) ;







 void Cramer(double *A,Coord3D *b,Coord3D *cc)
{//Solve a system of the Ax=b form, the solution is stored in Coord3D cc
//inputs are A[9] and Coord3D b;

double Det=A[0]*(A[4]*A[8]-A[7]*A[5])
          +A[3]*(-A[1]*A[8]+A[7]*A[2])
		  +A[6]*(A[1]*A[5]-A[2]*A[4]);

double Detx=b->x*(A[4]*A[8]-A[7]*A[5])
          +A[3]*(-b->y*A[8]+A[7]*b->z)
		  +A[6]*(b->y*A[5]-b->z*A[4]);

double Dety=A[0]*(b->y*A[8]-A[7]*b->z)
          +b->x*(-A[1]*A[8]+A[7]*A[2])
		  +A[6]*(A[1]*b->z-A[2]*b->y);

double Detz=A[0]*(A[4]*b->z-b->y*A[5])
          +A[3]*(-A[1]*b->z+b->y*A[2])
		  +b->x*(A[1]*A[5]-A[2]*A[4]);

cc->x=Detx/Det;
cc->y=Dety/Det;
cc->z=Detz/Det;

}

void GenerateRandom3D(Coord3D* p,int n)
{//function to generate random points in 0-1 range



// loop to generate random points not normalized
    int i;//counter;
	double tempmax=0;//temporary integer maximum random point;
	srand( (unsigned)time( NULL ) );



    for( i = 0;   i <n;i++ )
	{
		//X
		p[i].x=rand();
		if (p[i].x>tempmax)
		{tempmax=p[i].x;}

		//Y
		p[i].y=rand();
		if (p[i].y>tempmax)
		{tempmax=p[i].y;}

		//Z
		p[i].z=rand();
		if (p[i].z>tempmax)
		{tempmax=p[i].z;}



	}


    //loop to normalize
    for( i = 0;   i <n;i++ )
	{
		p[i].x=p[i].x/tempmax;

		p[i].y=p[i].y/tempmax;

		p[i].z=p[i].z/tempmax;
	}
}

void CC(Coord3D* p1,Coord3D* p2,Coord3D* p4,Coord3D* n,double *r,Coord3D* cc)

{ //Ritorna il quadrato del raggio del sphera D2.5D

	//Versione non ottimizzata per il calcolo del circocentro
	//possibili ottimizzazioni:
	//evitare il calcolo di v21 e v31 andando direttamente a scrivere il sistema a
	// scrittura simbolica del sistema, qualcosa si può semplificare sicuro

	//Questa routine viene usata solo nello start front, nel run front uso una versione "inline"
	//che evita di calcolare le normali 2 volte

 double A[9];
 Coord3D b,v21,v41,m;
 double rtemp;

 //get tnorm
DiffPoints(p2->,p1->,v21.);
DiffPoints(p4->,p1->,v41.);
CrossProduct(v41.,v21.,n->);
Normalize(n->,rtemp);

 //Assembling the system
 A[0]=n->x;
 A[1]=v21.x;
 A[2]=v41.x;

  A[3]=n->y;
 A[4]=v21.y;
 A[5]=v41.y;

  A[6]=n->z;
 A[7]=v21.z;
 A[8]=v41.z;

 //vettore dei termini noti
 DotProduct(n->,p1->,b.x);//primo elemento vettore termini noti

 MidPoint(p1->,p2->,m.);//salva in m il puntomedio

 DotProduct(v21.,m.,b.y);//secondo elemento vettore termini noti

 MidPoint(p1->,p4->,m.);//salva in m il puntomedio

 DotProduct(v41.,m.,b.z);//terzo elemento vettore termini noti


//Solve the system
Cramer(&A[0],&b,cc);

//CircumRadius
SquaredDistance(cc->,p1->,*r);
//take the minum distance
//avoids the triangle to self destruct during Dtest
SquaredDistance(cc->,p2->,rtemp);
if (rtemp<*r)
{*r=rtemp;}
SquaredDistance(cc->,p4->,rtemp);
if (rtemp<*r)
{*r=rtemp;}



}


double TriangleAngle(Coord3D* P1,Coord3D* P2,Coord3D* P3,Coord3D* P4)
{
//Points orientation
//    p1
//   / | \
// p3  |  p4
//   \ | /
//    p2
  //triangoli adiacenti all'edge

  Coord3D v12,v31,v41,tnorm1,tnorm2;
  double leng;//leng to normalize
  double alpha;






  //trova i vettori
  DiffPoints(P1->,P2->,v12.);//find the vector p1-p2
  DiffPoints(P3->,P1->,v31.);//find the vector p3-p1
  DiffPoints(P4->,P1->,v41.);//find the vector p3-p1

  //calcola le normali
  CrossProduct(v12.,v31.,tnorm1.);
   CrossProduct(v12.,v41.,tnorm2.);

   //normalizza

  Normalize(tnorm1.,leng);
  Normalize(tnorm2.,leng);


  //trova l'angolo con il DotProduct
 DotProduct(tnorm1.,tnorm2.,alpha);

 return alpha;
}



void SearchPoint(Coord3D* P1,Coord3D* P2,Coord3D* tnorm,Coord3D* sp,double* sr,double kr)
{



	//memory allocation

	Coord3D v21,cosdir;
    double leng;


  //trova i vettori
  DiffPoints(P2->,P1->,v21.);//find the vector p2-p1


   //calcola la direzione del search point
 CrossProduct(v21.,tnorm->,cosdir.);

   //normalizza

  Normalize(cosdir.,leng);

   //trova il punto medio dell'edge (ricicla v21)
  MidPoint(P1->,P2->,v21.);

  //sqrt(3)/2=0.866025403784439

  //Calcola il search radius come lunghezza dell'edge
  Distance(P1->,P2->,*sr);

  *sr=(*sr+*sr*kr)*.5;
  //calcola il search point in grandendo un pò sqrt(3)/2=0.866025403784439
  sp->x=v21.x+cosdir.x**sr*0.866025403785;
  sp->y=v21.y+cosdir.y**sr*0.866025403785;
  sp->z=v21.z+cosdir.z**sr*0.866025403785;


  Distance(P1->,sp->,*sr);//ricalcola il search radius per il caso maxr>1
  *sr=*sr*.9999999;
}


void SearchPointExt(Coord3D* P1,Coord3D* P2,Coord3D* tnorm,Coord3D* sp,double* sr)
{
//Points orientation
//    p1
//   / | \
// p3  |  sp
//   \ | /
//    p2


	//memory allocation

	Coord3D v21,cosdir;
    double leng;


  //trova i vettori
  DiffPoints(P2->,P1->,v21.);//find the vector p2-p1


   //calcola la direzione del search point
 CrossProduct(v21.,tnorm->,cosdir.);

   //normalizza

  Normalize(cosdir.,leng);

   //trova il punto medio dell'edge (ricicla v21)
  MidPoint(P1->,P2->,v21.);

  //sqrt(3)/2=0.866025403784439

  //Calcola il search radius come lunghezza dell'edge
  Distance(P1->,P2->,*sr);

  //*sr=*sr*kr;
  //calcola il search point in grandendo un pò sqrt(3)/2=0.866025403784439
  sp->x=v21.x+cosdir.x**sr*0.866025403785;
  sp->y=v21.y+cosdir.y**sr*0.866025403785;
  sp->z=v21.z+cosdir.z**sr*0.866025403785;


  Distance(P1->,sp->,*sr);//ricalcola il search radius per il caso maxr>1
  *sr=*sr*.9999999;
}


void TNorm(Coord3D* p1,Coord3D* p2,Coord3D* p4,Coord3D* tnorm)
    {

     Coord3D v21,v41;
     double leng;
 //computing the normal of the triagnle
                    //Get normal of new triangle
                      DiffPoints(p2->, p1->, v21.);
                      DiffPoints(p4->, p1->, v41.);
                      CrossProduct(v21., v41., tnorm->);
                      Normalize(tnorm->, leng);//la normale potrebbe essere normalizzata  afine dtest pensaci!!!
    }


inline void FastCramer(double* v1, double* v2, double* v3, double* b, double* sol) {//Solve a linear system. v1 v2 v3 are the system rows, b is the coefficent vector and sol is the computed solution


    double detx, dety, detz, det;
    double det01, det12, det20;

    det01=v2[0]*v3[1]-v2[1]*v3[0];
    det12=v2[1]*v3[2]-v2[2]*v3[1];
    det20=v2[2]*v3[0]-v2[0]*v3[2];

    det=v1[0]*det12+
            v1[1]*det20+
            v1[2]*det01;

#ifdef _DEBUG
    if (det<1e-14 && det>-1e-14)cout<<"Small det Found: "<< det<<endl;
#endif

    detx=b[0]*det12+
            -b[1]*(v1[1]*v3[2]-v1[2]*v3[1])+
            b[2]*(v1[1]*v2[2]-v1[2]*v2[1]);
    dety=b[0]*det20
            -b[1]*(v1[2]*v3[0]-v1[0]*v3[2])+
            b[2]*(v1[2]*v2[0]-v1[0]*v2[2]);
    detz=b[0]*det01
            -b[1]*(v1[0]*v3[1]-v1[1]*v3[0])+
            b[2]*(v1[0]*v2[1]-v1[1]*v2[0]);

    sol[0]=detx/det;
    sol[1]=dety/det;
    sol[2]=detz/det;

}
#endif
