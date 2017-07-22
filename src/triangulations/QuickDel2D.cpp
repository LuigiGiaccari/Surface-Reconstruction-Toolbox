#include "triangulations/QuickDel2D.h"




//TODO
//-p4 non può creare conflitti se appartiene al service
//-La macro swapTriangle può essere ottimizzata. per T1 sappiamo qual'è il triangolo da swappare
//-trasformare ct+1 e ct+2 in ct1 e ct2
//-controllare quanta memoria viene allocate per i triangoli
//-testare il point location come quello di shechuck (con la linea divisioria in caso gli edge sono entrambi positivi all'orient test)
//-testare la remove service come shcechuck (girando sul perimetro della triangolazione di servizio)
//-testare la poin location ONEDGE e un addpoint specifico
//-Ottimizzare i parametri HPCO e identificare eventualmente metodi più efficienti

//SOME MACROS


//idt had T1 now it has T2
#define SwapTriangles(idt, T1, T2)\
    if(idt>=0)\
    { if(t[idt].t1==T1){t[idt].t1=T2;}\
        else if(t[idt].t2==T1){t[idt].t2=T2;}\
        else {t[idt].t3=T2;}\
    }





//constructor
QUICKDEL2D::QUICKDEL2D()  //nitialize all data
{
    ct = 0;
    np = 0;
    p = NULL;
    t = NULL;
    toflip = NULL;
    flipstack = NULL;
    stack_iter = 0;
    first_tria = 0;
    Options_OverWritePoints = true; //authorize to change  points array
    Options_UseHCPO = true; //enables the HPCO sorting
    Options_UseBrute = true; //enables brute triangulation
    Options_Reindex = true; //if true preserve indexes of HPCO sort
    Options_CheckTriangulation = false; //checks triagulation at the end
#ifdef STATS
    Stats.incircle_call = 0;
    Stats.orient_call = 0;
    Stats.incircle_support_call = 0;
    Stats.flips = 0;
#endif

}

//destructor
QUICKDEL2D::~QUICKDEL2D()
{
    FreeMemory();
}



void QUICKDEL2D::FreeMemory()
{
    //Deallocate(&p);

    Deallocate(&toflip);
    Deallocate(&flipstack);
    Deallocate(&t);
}

void QUICKDEL2D::ComputeBoundingBox()
{

    MinMax(&p[0].x, np, &xmax, &xmin, 2); //minmax x
    MinMax(&p[0].y, np, &ymax, &ymin, 2); //minmax y
}


void QUICKDEL2D::PreProcess(double* inputp, int inputnp)
{
    //preliminary operations like:
    //  -converting data structure
    //  -Allocate memory
    //	-finding the bounding box
    //  -initialize routines for exact arithmetics
    //  -build service triangulation


    //copy to class data structure
    np = inputnp;
    p = (Coord2D*)inputp; //pass double sìto Coord2D data structures
    const int maxt = np * 2 + 1; //maximum number of triangles (the first one plus 2 for every point added)

    //allocate memory
    Allocate(&t, maxt);
    Alls(&t[0].p1, maxt * 6, -1); //initiate triangualtion to -1

    ComputeBoundingBox();// get bounding box

    BuildServiceTriangulation();





#ifdef ROBUST
    exactinit();//initilaize routines for exact aritmetics
#endif
}

void QUICKDEL2D::BuildServiceTriangulation()
{
    // adds the first triangle that inclose all other points
    //add the service triangulation (notice this requires an extra space inside
    double W = xmax - xmin;
    double H = ymax - ymin;

    if (H > W)
    {
        W = H;    //finds the largest dimension
    }

    //we must ensure the dataset do not fall outside of the service

    //first 3 points in counterclockwise order
    p[np].x = .5 * (xmax + xmin);
    p[np].y = ymax + 3 * W;
    p[np + 1].x = xmin - 3 * W;
    p[np + 1].y = ymin - 1 * W;
    p[np + 2].x = xmax + 3 * W;
    p[np + 2].y = ymin - 1 * W;

    t[0].p1 = np;
    t[0].p2 = np + 1;
    t[0].p3 = np + 2;
    ct = 1; //we now have one triangle


}

//locate the enclosing triangle for each point.
// this function detecs duplicated points and points lying on edges
// the return flag is:
//
// -0 points isnide simplex
// -1 duplicate point
//  1 point on edge
//
int QUICKDEL2D::LocatePoint(int* cT, int i)
{

    // cT current triangle
    // i point to locate




    //int counter=0;//counter to avoid overloop
    double o1, o2, o3;//orientation

    int pe2 = -1; //point of the edge where we come from.
    //it is used to skip the orient call on the same edge where we come from
    //notice that since edges are sorted clockwise, we just need a point to recognize an edge
    //initialize to -1 for first iteration

    while (1)
    {

        int p1 = t[*cT].p1;
        int p2 = t[*cT].p2;
        int p3 = t[*cT].p3;

        o1 = o2 = o3 = 1; //initialize for duplicate points check(BUG corretto: prima era fuori dal while e veniva inizializzato solo all'inizio)
        //check edge 1 (p2,p3)

        if (p1 != pe2)
        {

#ifdef ROBUST
            o1 = orient2d(&p[p1].x, &p[p2].x, &p[i].x );
#else
            o1 = orient2dfast(&p[p1].x, &p[p2].x, &p[i].x );
#endif

            if (o1 < 0)
            {
                pe2 = p2;
                *cT = t[*cT].t1;
                if (*cT >= 0)
                {
                    continue;   //%go to next triangle
                }
                else
                {
                    Error("Error in point location");
                }
            }
        }
        //check edge 1 (p2,p3)
        if (p2 != pe2)
        {
#ifdef ROBUST
            o2 = orient2d(&p[p2].x, &p[p3].x, &p[i].x);
#else
            o2 = orient2dfast(&p[p2].x, &p[p3].x, &p[i].x);
#endif


            if (o2 < 0)
            {
                pe2 = p3;
                *cT = t[*cT].t2;
                if (*cT >= 0)
                {
                    continue;   //%go to next triangle
                }
                else
                {
                    Error("Error in point location");
                }
            }
        }


        //check edge 1 (p3,p1)

        if (p3 != pe2)
        {
#ifdef ROBUST
            o3 = orient2d(&p[p3].x, &p[p1].x, &p[i].x);
#else
            o3 = orient2dfast(&p[p3].x, &p[p1].x, &p[i].x);
#endif

            if (o3 < 0)
            {
                pe2 = p1;
                *cT = t[*cT].t3;
                if (*cT >= 0)
                {
                    continue;   //%go to next triangle
                }
                else
                {
                    Error("Error in point location");
                }
            }
        }

        //check for duplicate points can be optimized
        if ( o1 == 0 || o2 == 0 || o3 == 0) //check for some degenerancies
        {
            if (o1 == 0 && o2 == 0 || o1 == 0 && o3 == 0 || o2 == 0 && o3 == 0) //duplicate point !
            {
                Warning("duplicate point found\n");
                /*printf("i: %4.0d %4.4f %4.4f \n",i,p[i].x,p[i].y);
                printf("p1: %4.0d %4.4f %4.4f \n",p1,p[p1].x,p[p1].y);
                printf("p2: %4.0d %4.4f %4.4f \n",p2,p[p2].x,p[p2].y);
                printf("p3: %4.0d %4.4f %4.4f \n",p3,p[p3].x,p[p3].y);
                 printf("o1:%4.15f \n",o1);
                  printf("o2:%4.15f \n",o2);
                   printf("o3:%4.15f \n",o3);*/

                if (o1 == 0 && o2 == 0 && o3 == 0)
                {
                    Error("Degenerate Facet during Delaunay Triangulation");    //il triangolo è composto da 3 punti colineari
                }

                return -1;//return duplicate points flag
            }
            else if (o1 == 0) //points lying on edge
            {
                return 1;   //return lying on edge 1 flag
            }
            else if (o2 == 0) //points lying on edge
            {
                return 2;   //return lying on edge 1 flag
            }
            else if (o3 == 0) //points lying on edge
            {
                return 3;   //return lying on edge 1 flag
            }
        }
        return 0; //we found the enclosing simplex

    }

}

//Adds the point i inside the triangle idt and checks for flips if bool flip=true
void QUICKDEL2D::AddPoint(int idt, int i, bool flip)
{
    //Triangle idt before isnertion of the point i
    //	  p1
    //T1  /\  T3
    //   /i \
    //p2/____\p3
    //    T2
    //T1 is the triagnle belonging t edge p1-p2; T2->p2-p3; T3->p3-p1
    //        p1
    //       /|\
    //      / | \
    //     /  |  \
    //    / / i\  \
    //   //      \ \
    //p2/___________\p3


    //The data structure was
    //   T[idt]=[ p1,p2,p3,T1,T2,T3]
    //and  becomes
    //   T[idt]=[p1,p2,i,ct,idt,ct1]
    //   T[ct]=[p2,p3,i,T2,ct1,idt]
    //   T[ct1]=[p3,p1,i,T3,idt,ct]

    /* Unoptimized swaps
     * //the third new triangle
     * int p1,p2,p3,T1,T2,T3;
     * p1=t[idt].p1;;p2=t[idt].p2;p3=t[idt].p3;T1=t[idt].t1;T2=t[idt].t2;T3=t[idt].t3;
     * //first we build 3 new triangles
     * //the first one replace idt;
     * t[idt].p1=p1;t[idt].p2=p2;t[idt].p3=i;
     * t[idt].t1=T1;
     * t[idt].t2=ct;t[idt].t3=ct1;
     * //the second one
     * t[ct].p1=p2;t[ct].p2=p3;t[ct].p3=i;
     * t[ct].t1=T2;t[ct].t2=ct1;t[ct].t3=idt;
     * //the third one
     * t[ct1].p1=p3;t[ct1].p2=p1;t[ct1].p3=i;
     * t[ct1].t1=T3;t[ct1].t2=idt;t[ct1].t3=ct;
     *
     *
     * //T2 had idt now it has ct
     * SwapTriangles(T2,idt,ct)
     *
     * //T3 had idt now it has ct1
     * SwapTriangles(T3,idt,ct1)
     *
     */

    //Optimized version


    int ct1 = ct + 1; //local variable for ct+1
    //the third one
    t[ct1].p1 = t[idt].p3;
    t[ct1].p2 = t[idt].p1;
    t[ct1].p3 = i; //points of the 3th triangle
    t[ct1].t1 = t[idt].t3;
    t[ct1].t2 = idt;
    t[ct1].t3 = ct; //neigh triagnles of the 3th triangle

    //the second one
    t[ct].p1 = t[idt].p2;
    t[ct].p2 = t[idt].p3;
    t[ct].p3 = i; //points of the 2th triangle
    t[ct].t1 = t[idt].t2;
    t[ct].t2 = ct1;
    t[ct].t3 = idt; //points of the 2th triangle



    //T2 had idt now it has ct
    SwapTriangles(t[idt].t2, idt, ct)

    //T3 had idt now it has ct1
    SwapTriangles(t[idt].t3, idt, ct1)


    //the first one replace idt; (//change the first at last)
    //t[idt].p1=p1;t[idt].p2=p2;
    t[idt].p3 = i;
    //t[idt].t1=T1;
    t[idt].t2 = ct;
    t[idt].t3 = ct1;



    //end optimized version



    // flip the edges if needed
    if (flip)
    {
        CheckForFlip(idt);
        CheckForFlip(ct);
        CheckForFlip(ct1);
    }
    else//update stack
    {
        flipstack[stack_iter] = idt;
        stack_iter++;//add triagnle to stack
        flipstack[stack_iter] = ct;
        stack_iter++;//add triagnle to stack
        flipstack[stack_iter] = ct1;
        stack_iter++;//add triagnle to stack
        toflip[ct] = toflip[ct1] = toflip[idt] = true; //this triangles needs to checked for flip
    }


    ct = ct + 2; //two more triangles

}

void QUICKDEL2D::Flip(int T1)  //Flip the triangle T1 without incircle test
{


    //Assigning local variables (it can be avoided but it is much more readeable and just alittle slower)
    int p1, p2, p3, p4, T2, T3, T4, T5, T6;
    p1 = t[T1].p1;
    p2 = t[T1].p2;
    p3 = t[T1].p3;
    T2 = t[T1].t1; //T3 T4 moved after the incircle




    //find the position of p2 to get the position of p4 T5 T6
    //now we only need p4, this comparison can be repetead after the incircle (slightly optimized)
    if (t[T2].p1 == p2)
    {
        p4 = t[T2].p3;
        T5 = t[T2].t2;
        T6 = t[T2].t3;
    }
    else if (t[T2].p2 == p2)
    {
        p4 = t[T2].p1;
        T5 = t[T2].t3;
        T6 = t[T2].t1;
    }
    else
    {
        p4 = t[T2].p2;
        T5 = t[T2].t1;
        T6 = t[T2].t2;
    }








    T3 = t[T1].t2;
    T4 = t[T1].t3;

    //flip!
    //t[T1].p1=p1;
    t[T1].p2 = p4;
    //t[T1].p3=p3;
    t[T1].t1 = T5;
    t[T1].t2 = T2;
    //t[T1].t3=T4;
    t[T2].p1 = p4;
    t[T2].p2 = p2;
    t[T2].p3 = p3;
    t[T2].t1 = T6;
    t[T2].t2 = T3;
    t[T2].t3 = T1;

    //we also need to update something into T4,T5

    //T3 had T1 now it has T2
    SwapTriangles(T3, T1, T2)

    // T5 had T2 now it has T1
    SwapTriangles(T5, T2, T1)

    //we also need to check p1-p4 & p4-p2
    CheckForFlip(T1);//recursive call
    CheckForFlip(T2);//recursive call


}

//check for flips the edge between triangle T1 & T2
// the code is optized in order to have the flip always on the first edge of T1 (p1-p2)
void QUICKDEL2D::CheckForFlip(int T1)
{
    //before
    //	  p3
    //    /\
    // T4/T1\T3
    //p1/____\p2
    //  \    /
    // T5\T2/T6
    //    \/
    //	  p4

    //after


    //we may end here
    if (t[T1].t1 < 0) //no triangle
    {
        return;
    }

    //Assigning local variables (it can be avoided but it is much more readeable and just alittle slower)
    int p1, p2, p3, p4, T2, T3, T4, T5, T6;
    p1 = t[T1].p1;
    p2 = t[T1].p2;
    p3 = t[T1].p3;
    T2 = t[T1].t1; //T3 T4 moved after the incircle




    //find the position of p2 to get the position of p4 T5 T6
    //now we only need p4, this comparison can be repetead after the incircle (slightly optimized)
    if (t[T2].p1 == p2)
    {
        p4 = t[T2].p3;
        T5 = t[T2].t2;
        T6 = t[T2].t3;
    }
    else if (t[T2].p2 == p2)
    {
        p4 = t[T2].p1;
        T5 = t[T2].t3;
        T6 = t[T2].t1;
    }
    else
    {
        p4 = t[T2].p2;
        T5 = t[T2].t1;
        T6 = t[T2].t2;
    }


    if (p4 > np)
    {
        return;    //service points wont make the incircle fail
    }


#ifdef STATS
    if (p1 >= np || p2 >= np || p3 >= np) //we are flipping a support edge
    {
        Stats.incircle_support_call++;
    }
#endif

    // incircle test
# ifdef ROBUST
    double r = incircle( &p[p1].x, &p[p2].x, &p[p3].x, &p[p4].x);
#else
    double r = incirclefast( &p[p1].x, &p[p2].x, &p[p3].x, &p[p4].x);
#endif

    //DEBUG LINES
    /*if(orient2d( &p[p1].x, &p[p2].x, &p[p4].x)==0.0)
    {printf("point on edge returned r= %4.15f \n",r);}
    if(orient2d( &p[p1].x, &p[p2].x, &p[p3].x)==0.0)
    {printf("Degenerate Facet r= %4.15f \n",r);
     printf("p1  %4.15f %4.15f\n",p[p1].x,p[p1].y);
     printf("p2  %4.15f %4.15f\n",p[p2].x,p[p2].y);
     printf("p3  %4.15f %4.15f\n",p[p3].x,p[p3].y);
      printf("p4  %4.15f %4.15f\n",p[p4].x,p[p4].y);
     Error("Degenerate facet");
    }*/
    //END DEBUG LINES

    //points is inside a flip is requested
    if (r > 0)
    {
#ifdef STATS
        Stats.flips++;
#endif

        T3 = t[T1].t2;
        T4 = t[T1].t3;

        //flip!
        //t[T1].p1=p1;
        t[T1].p2 = p4;
        //t[T1].p3=p3;
        t[T1].t1 = T5;
        t[T1].t2 = T2;
        //t[T1].t3=T4;
        t[T2].p1 = p4;
        t[T2].p2 = p2;
        t[T2].p3 = p3;
        t[T2].t1 = T6;
        t[T2].t2 = T3;
        t[T2].t3 = T1;

        //we also need to update something into T4,T5

        //T3 had T1 now it has T2
        SwapTriangles(T3, T1, T2)

        // T5 had T2 now it has T1
        SwapTriangles(T5, T2, T1)

        //we also need to check p1-p4 & p4-p2
        CheckForFlip(T1);//recursive call
        CheckForFlip(T2);//recursive call
    }
}

void QUICKDEL2D::AddPointOnEdge(int idt, int i, int idedge) //Adds point lying on edge
{
    int ct1 = ct + 1; //local variable for ct+1
    //the third one
    t[ct1].p1 = t[idt].p3;
    t[ct1].p2 = t[idt].p1;
    t[ct1].p3 = i; //points of the 3th triangle
    t[ct1].t1 = t[idt].t3;
    t[ct1].t2 = idt;
    t[ct1].t3 = ct; //neigh triagnles of the 3th triangle

    //the second one
    t[ct].p1 = t[idt].p2;
    t[ct].p2 = t[idt].p3;
    t[ct].p3 = i; //points of the 2th triangle
    t[ct].t1 = t[idt].t2;
    t[ct].t2 = ct1;
    t[ct].t3 = idt; //points of the 2th triangle



    //T2 had idt now it has ct
    SwapTriangles(t[idt].t2, idt, ct)

    //T3 had idt now it has ct1
    SwapTriangles(t[idt].t3, idt, ct1)


    //the first one replace idt; (//change the first at last)
    //t[idt].p1=p1;t[idt].p2=p2;
    t[idt].p3 = i;
    //t[idt].t1=T1;
    t[idt].t2 = ct;
    t[idt].t3 = ct1;



    //end optimized version



    // flip the edges if needed

    if (idedge != 1)
    {
        CheckForFlip(idt);
    }
    else
    {
        Flip(idt);
    }

    if (idedge != 2)
    {
        CheckForFlip(ct);
    }
    else
    {
        Flip(ct);
    }

    if (idedge != 3)
    {
        CheckForFlip(ct1);
    }
    else
    {
        Flip(ct1);
    }


    ct = ct + 2; //two more triangles

}

void QUICKDEL2D::FlipThemAll()//flips all those triagnles that need to be flipped
{
    int idt, i;

    //now flip from first_tria to last_tria

    //for(idt=first_tria;idt<ct;idt++){CheckForFlip_Global(idt);}
    //first_tria=ct;//first_tria for next FlipThemAll call

    //flip those in the stack
    for (i = 0; i < stack_iter; i++)
    {
        CheckForFlip_Global(flipstack[i]);
    }
    stack_iter = 0; //stack is now empty


    //DEBUGLINES
    //if (IsDelaunay()){cout<<"After flip Triangulation is delaunay"<<endl;}
    //else {cout<<"After flip Triangulation is NOT delaunay"<<endl;}

    /*DEBUG LINES (see if some triangle is still to flip
    for(i=0;i<ct;i++)
    {
     if(toflip[i])
     {
    	 Warning("Triangle has not been checked for flip");
         cout<<"Id triangle="<<i<<endl;
    	 break;
     }
    }
    */
}


//check for flips the edge between triangle T1 & T2
// the code is optized in order to have the flip always on the first edge of T1 (p1-p2)
void QUICKDEL2D::CheckForFlip_Global(int T1)
{
    //before
    //	  p3
    //    /\
    // T4/T1\T3
    //p1/____\p2
    //  \    /
    // T5\T2/T6
    //    \/
    //	  p4

    //after
    toflip[T1] = false; //no more to flip

    //we may end here
    if (t[T1].t1 < 0) //no triangle
    {
        return;
    }

    //we may end here
    if (toflip[t[T1].t1])//this triangle will be checked for flip by his neighbour
    {
        return;
    }

    //Assigning local variables (it can be avoided but it is much more readeable and just alittle slower)
    int p1, p2, p3, p4, T2, T3, T4, T5, T6;
    p1 = t[T1].p1;
    p2 = t[T1].p2;
    p3 = t[T1].p3;
    T2 = t[T1].t1; //T3 T4 moved after the incircle




    //find the position of p2 to get the position of p4 T5 T6
    //now we only need p4, this comparison can be repetead after the incircle (slightly optimized)
    if (t[T2].p1 == p2)
    {
        p4 = t[T2].p3;
        T5 = t[T2].t2;
        T6 = t[T2].t3;
    }
    else if (t[T2].p2 == p2)
    {
        p4 = t[T2].p1;
        T5 = t[T2].t3;
        T6 = t[T2].t1;
    }
    else
    {
        p4 = t[T2].p2;
        T5 = t[T2].t1;
        T6 = t[T2].t2;
    }


    if (p4 > np)
    {
        return;    //service points wont make the incircle fail
    }


#ifdef STATS
    if (p1 >= np || p2 >= np || p3 >= np) //we are flipping a support edge
    {
        Stats.incircle_support_call++;
    }
#endif

    // incircle test
# ifdef ROBUST
    double r = incircle( &p[p1].x, &p[p2].x, &p[p3].x, &p[p4].x);
#else
    double r = incirclefast( &p[p1].x, &p[p2].x, &p[p3].x, &p[p4].x);
#endif


    //points is inside a flip is requested
    if (r > 0)
    {
#ifdef STATS
        Stats.flips++;
#endif

        T3 = t[T1].t2;
        T4 = t[T1].t3;

        //flip!
        //t[T1].p1=p1;
        t[T1].p2 = p4;
        //t[T1].p3=p3;
        t[T1].t1 = T5;
        t[T1].t2 = T2;
        //t[T1].t3=T4;
        t[T2].p1 = p4;
        t[T2].p2 = p2;
        t[T2].p3 = p3;
        t[T2].t1 = T6;
        t[T2].t2 = T3;
        t[T2].t3 = T1;

        //we also need to update something into T4,T5

        //T3 had T1 now it has T2
        SwapTriangles(T3, T1, T2)

        // T5 had T2 now it has T1
        SwapTriangles(T5, T2, T1)

        //we also need to check p1-p4 & p4-p2
        CheckForFlip_Global(T1);//recursive call
        CheckForFlip_Global(T2);//recursive call
    }
}











#ifdef STATS
void QUICKDEL2D::PrintStatistics()
{
    //compute the ratio depending on the number of points
    float incircle_call_ratio = float(Stats.incircle_call) / float(np);
    float orient_call_ratio = float(Stats.orient_call) / float(np);
    float incircle_support_call_ratio = float(Stats.incircle_support_call) / float(Stats.incircle_call);
    float flips_ratio = float(Stats.flips) / float(np);

    printf(" Incircle calls ratio: %4.4f\n", incircle_call_ratio);
    printf(" Incircle suppport calls ratio: %4.4f\n", incircle_support_call_ratio);
    printf(" Flips ratio: %4.4f\n", flips_ratio);
    printf(" Orient calls ratio: %4.4f\n", orient_call_ratio);

}
#endif


void QUICKDEL2D::PrintOutput_mfile()
{
    int i;

    ofstream myfile;
    myfile.open ("Output.m");
    //POINTS

    myfile << "p=[\n";

    for (i = 0; i < np + 3; i++)
    {
        myfile << p[i].x << " " << p[i].y << endl;
    }
    myfile << "];\n";

    //TRIANGLES
    myfile << "t=[\n";

    for (i = 0; i < ct; i++)
    {
        myfile << t[i].p1 + 1 << " " << t[i].p2 + 1 << " " << t[i].p3 + 1 << endl; //plus one for matlab notation
    }
    myfile << "];\n" << endl;
    myfile << "np=" << np << ";" << endl;
    myfile << "ind=any(t(:,1:3)>np,2);" << endl;
    myfile << " t(ind,:)=[];" << endl;

    myfile << "triplot(t,p(:,1),p(:,2),'color','g');" << endl;

    myfile.close();



}
bool QUICKDEL2D::IsDelaunay()//cheks wheter the triangulation is delaunay
{
    int i;
    //loop trough all triangles
    for (i = 0; i < ct; i++)
    {
        if (IsInside(i, t[i].t1))
        {
            return false;
        }
        if (IsInside(i, t[i].t2))
        {
            return false;
        }
        if (IsInside(i, t[i].t3))
        {
            return false;
        }
    }



    return true;
}


bool QUICKDEL2D::IsInside(int T1, int T2) //checkes if the apex of T2 is inside the circumcircle of T1
{
    //we may end here
    if (T2 < 0) //no triangle
    {
        return false;
    }

    //Assigning local variables (it can be avoided but it is much more readeable and just alittle slower)
    int p1, p2, p3, p4;
    p1 = t[T1].p1;
    p2 = t[T1].p2;
    p3 = t[T1].p3;




    //Get P4 (the apex)
    if (t[T2].p1 == p2)
    {
        p4 = t[T2].p3;
    }
    else if (t[T2].p2 == p2)
    {
        p4 = t[T2].p1;
    }
    else
    {
        p4 = t[T2].p2;
    }


    if (p4 > np)
    {
        return false;    //service points wont make the incircle fail
    }

    // incircle test
# ifdef ROBUST
    double r = incircle( &p[p1].x, &p[p2].x, &p[p3].x, &p[p4].x);
#else
    double r = incirclefast( &p[p1].x, &p[p2].x, &p[p3].x, &p[p4].x);
#endif
    if (r > 0)
    {
        Warning("Triangulations is not Delaunay");
        cout << "Triangle " << T1 << " conflicts with point " << p4 << endl;
        return true;
    }
    else
    {
        return false;
    }
}

void QUICKDEL2D::RemoveServiceTriangulation()//removes service triangulation
{}


int QUICKDEL2D::Triangulate(double* inputp, int inputnp) //incremental delaunay triangulation wih HPCO
{
    int flag = 1; //default to error
    //Preprocessing of points
    if (Options_UseHCPO)
    {
        HCPO.sortHCPO2D((Coord2D*)inputp, inputnp, true); //the value of 2 is empirical and not optimized
    }



    if (!Options_Reindex)
    {
        HCPO.~HCPO2D();    //destroy object
    }


    //choose triangulation method
    if (Options_UseBrute)
    {
        flag = Triangulate_Brute(inputp, inputnp);
    }
    else
    {
        flag = Triangulate_Global(inputp, inputnp);
    }

    if (flag != 0)
    {
        return flag;
    }

    if (Options_CheckTriangulation)
    {
        if (IsDelaunay())
        {
            cout << "Triangulation is Delaunay" << endl;
        }
        else
        {
            cout << "Triangulation is not Delaunay" << endl;
        }
    }

#ifdef STATS
    if (Options_PrintStats)
    {
        PrintStatistics();
    }
#endif


    return flag;
}


int QUICKDEL2D::Triangulate_Brute(double* inputp, int inputnp)  //traignulate using a brute incremental algorithm
{
    //inputp must have (inputnp+3)*2 values  wheras +3 stands for the service triangulation
    int i, flag;
    int printstats = 100000;

    PreProcess(inputp, inputnp);//preprocessing operations

    int EnclosingT = 0; //current enclosing triagnle
    for (i = 0; i < np; i++)
    {
        //loop trough all points

        //finds the enlcosing simplex for the ith point
        flag = LocatePoint(&EnclosingT, i);
        if (flag < 0)
        {
            continue;   //duplicate point skip!
        }

        if (flag == 0)
        {
            AddPoint(EnclosingT, i);   //adds the poit i inside the triangle idt
        }
        else
        {
            AddPointOnEdge(EnclosingT, i, flag);
        }

        //AddPoint(EnclosingT, i);
        /*
        //Debug Lines
        if (i==printstats){
        	cout<<"i="<<i<<endl;
        	PrintStatistics();
        	printstats=printstats+1e5;
        }*/
    }

    RemoveServiceTriangulation();


    return 0;
}


int QUICKDEL2D::Triangulate_Global(double* inputp, int inputnp) //incremental delaunay triangulation wih HPCO
{
    //inputp must have (inputnp+3)*2 values  wheras +3 stands for the service triangulation
    int i, flag;
    int printstats = 1e5;

    PreProcess(inputp, inputnp);//preprocessing operations

    //Memory allocationpns only valid for Global method
    AllocateAndInit(&toflip, (np + 3) * 2 + 1, false);
    Allocate(&flipstack, (np + 3) * 2 + 1);

    int EnclosingT = 0; //current enclosing triagnle
    for (i = 0; i < np; i++)

    {
        //loop trough all points

        flag = LocatePoint(&EnclosingT, i); //finds the enlcosing simplex
        if (flag < 0)
        {
            continue;   //duplicate point skip!
        }

        if (toflip[EnclosingT])
        {
            FlipThemAll();//call flip them all
            flag = LocatePoint(&EnclosingT, i); //Relocate point
        }

        if (flag == 0)
        {
            AddPoint(EnclosingT, i);   //adds the poit i inside the triangle idt
        }
        else
        {
            AddPointOnEdge(EnclosingT, i, flag);
        }

        //AddPoint(EnclosingT, i);
        /*
        //Debug Lines
        if (i==printstats){
        cout<<"i="<<i<<endl;
        PrintStatistics();
        printstats=printstats+1e5;
        }*/
    }

    FlipThemAll();//last call to flip them all


    RemoveServiceTriangulation();


    return 0;
}
