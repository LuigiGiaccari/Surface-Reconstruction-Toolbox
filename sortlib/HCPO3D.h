//class to sort an array of doubles in HPCO (Hierarchical Column Prime Order)
//the insertions order helps point location in Delaunay triangulation
//
// it uses the sedgesort_index library
//
// index of sorting are stored in the index array until the distructor is called
//
//

//TODO
// -decide wheter to swap x y in case of "thin" bounding box
// -in case the overflow causes a last x sort of more or less nx/2 join or not xsorts
// -Optimize the reindex-phase

#include  <math.h>
#include "structs/Coord3D.h"
#include "structs/Coord2D.h"
#include "util/ArraysLib.h"
#include "sortlib/sedgesort_index.h"


using namespace std;

class HCPO3D
{
private:
    //variables

    bool xreverse,yreverse,zreverse;//defines wheter to reverse order of sorting
    double* temparray;//this array is a temporary store for all sorts we are going to perform
    Coord3D*p;//array of points to order
    int np;//number of points

    //functions
    void sortCPOyx(int i1,int i2);//sort points between index i1 and i2 in CPO order along x & y
    void sortCPOzyx(int i1,int i2);//sort points between index i1 and i2 in CPO order along z & x & y
    bool IsSorted(double* arr,int n);//debug routine to check wheter the array is sorted
    void sortY(int i1,int i2);//sort points in i1 i2 range in Y direction (versus is defined by yreverse)
    void sortX(int i1,int i2);//sort points in i1 i2 range in X direction (versus is defined by xreverse)
    void sortZ(int i1,int i2);//sort points in i1 i2 range in Z direction (versus is defined by zreverse)
    void Shuffle();//shuffle points they are not random
    void Reindex();//sorts the p array according to indexes
public:
    void sortHCPO2D(Coord3D*p,int np,bool shuffle);//sort 2D points in HPCO order
    void sortHCPO3D(Coord3D*p,int np,bool shuffle);//sort points in HPCO order
    void FreeMemory();//Free memory to reset the object status
    float nx_factor, nyx_factor;//relative factor on numbers of points related to dimensions proportion
    float nt_factor;
    int* index;//indexes of HPCO public cause useful for reindexing
    HCPO3D();//constructor
    ~HCPO3D();//desstructor
};


HCPO3D::HCPO3D()//initdata
{
    index=NULL;
    temparray=NULL;
    np=0;
    xreverse=yreverse=false;
    nx_factor=1;//x versus y subdivision factor
    nyx_factor=1;//x versus y subdivision factor
    nt_factor=1;//number of points to add versus umber of triagnles
}
HCPO3D::~HCPO3D()//desstructor
{
    FreeMemory();
}
void HCPO3D::FreeMemory()
{
    np=0;
    xreverse=yreverse=zreverse=false;
    nx_factor=1;//x versus y subdivision factor
    nt_factor=2;//number of points to add versus umber of triagnles
    Deallocate(&index);
    Deallocate(&temparray);
}

void HCPO3D::sortCPOyx(int i1,int i2)//sort points between index i1 and i2 in CPO order
{
    int nextstep;
    int start,end;
    double xmin,xmax,ymin,ymax;//bounding box

    if (i2==i1)return;//no need to sort one point

    int step=i2-i1+1;
    //Get bounding box
    MinMax(&p[i1].x,step,&xmax,&xmin,2);//minmax x
    MinMax(&p[i1].y,step,&ymax,&ymin,2);//minmax y
    int nx=nx_factor*sqrt(step*(ymax-ymin)/(xmax-xmin));//number of elements to sort in x x
    if (nx==0)nx=1;

    //sort points in y direction
    sortY(i1,i2);
    //sort buckets in x direction

    start=i1;
    end=i1+nx-1;
    while(1)
    {
        sortX(start,end);

        if(i2==end)break;//we are done
        start=end+1;
        end=start+nx-1;

        nextstep=end+nx;//predict overflow in the next step
        if (nextstep>i2)end=i2;
    }

}


void HCPO3D::sortCPOzyx(int i1,int i2)//sort points between index i1 and i2 in CPO order
{
    int nextstep;
    int start,end;
    double xmin,xmax,ymin,ymax,zmin,zmax;//bounding box
    int step=i2-i1+1;
    //Get bounding box
    MinMax(&p[i1].x,step,&xmax,&xmin,3);//minmax x
    MinMax(&p[i1].y,step,&ymax,&ymin,3);//minmax y
    MinMax(&p[i1].z,step,&zmax,&zmin,3);//minmax x

    double V=(zmax-zmin)*(xmax-xmin)*(ymax-ymin);//volume
    double squared_edge=pow(V/double(step),.33333333333333333);//dimensione del quadrato medio riferita al numero attuale di punti=step
    int nx=(xmax-xmin)/squared_edge+.5;
    if(nx==0)nx=1;
    int ny=(ymax-ymin)/squared_edge+.5;
    if(ny==0)ny=1;
    int nz=(zmax-zmin)/squared_edge+.5;
    if(nz==0)nz=1;






    int nyx=nyx_factor*nx*ny;//number of elements to sort in yx direction
    if (nyx==0)nyx=1;

    //sort points in z direction
    sortZ(i1,i2);
    //sort buckets in x direction

    start=i1;
    end=i1+nyx-1;
    while(1)
    {
        sortCPOyx(start,end);

        if(i2==end)break;//we are done
        start=end+1;
        end=start+nyx-1;

        nextstep=end+nyx;//predict overflow in the next step
        if (nextstep>i2)end=i2;
    }

}


void HCPO3D::sortX(int i1,int i2) //sort points in i1 i2 range in X direction (versus is defined by xreverse)
{
    int i,step;
    //First copies x value
    if(xreverse)//copies with negative sign
    {
        for(i=i1; i<=i2; i++)
        {
            temparray[i]=-p[index[i]].x;
        }
    }
    else
    {
        for(i=i1; i<=i2; i++)
        {
            temparray[i]=p[index[i]].x;
        }
    }

    step=i2-i1+1;//number of elements we have to sort
    temparray[i2+1]=HUGE_VAL;//insert sentinel for sedge sort
    sedgesort_index (&temparray[i1],&index[i1],step);//call sorting routine
    xreverse=!xreverse;//invert xreverse

    //DEBUG LINES
    //if (!IsSorted(&temparray[i1],step))Error("Unsorted array");
    //cout<<"sortx:  ";
    //for(i=i1;i<=i2;i++)
    // {cout<<p[index[i]].y<<" ";}
    // cout<<endl;
}





void HCPO3D::sortY(int i1,int i2)//sort points in i1 i2 range in X direction (versus is defined by xreverse)
{
    int i,step;
    //First copies x value
    if(yreverse)//copies with negative sign
    {
        for(i=i1; i<=i2; i++)
        {
            temparray[i]=-p[index[i]].y;
        }
    }
    else
    {
        for(i=i1; i<=i2; i++)
        {
            temparray[i]=p[index[i]].y;
        }
    }

    step=i2-i1+1;//number of elements we have to sort
    temparray[i2+1]=HUGE_VAL;//insert sentinel for sedge sort
    sedgesort_index (&temparray[i1],&index[i1],step);//call sorting routine
    yreverse=!yreverse;//invert yreverse

    //DEBUG LINES
    //if (!IsSorted(&temparray[i1],step))Error("Unsorted array");
    //cout<<"sorty:  ";
    //for(i=i1;i<=i2;i++)
    // {cout<<p[index[i]].y<<" ";}
    // cout<<endl;
}

void HCPO3D::sortZ(int i1,int i2)//sort points in i1 i2 range in Z direction (versus is defined by Zreverse)
{
    int i,step;
    //First copies x value
    if(zreverse)//copies with negative sign
    {
        for(i=i1; i<=i2; i++)
        {
            temparray[i]=-p[index[i]].z;
        }
    }
    else
    {
        for(i=i1; i<=i2; i++)
        {
            temparray[i]=p[index[i]].z;
        }
    }

    step=i2-i1+1;//number of elements we have to sort
    temparray[i2+1]=HUGE_VAL;//insert sentinel for sedge sort
    sedgesort_index (&temparray[i1],&index[i1],step);//call sorting routine
    zreverse=!zreverse;//invert yreverse

    //DEBUG LINES
    //if (!IsSorted(&temparray[i1],step))Error("Unsorted array");
    //cout<<"sorty:  ";
    //for(i=i1;i<=i2;i++)
    // {cout<<p[index[i]].y<<" ";}
    // cout<<endl;
}


void HCPO3D::sortHCPO2D(Coord3D*inputp,int inputnp,bool shuffle)//sort points in HCPO order
{
    int start,end,nextstep,i,current_nt;

    //Assign variables
    p=inputp;
    np=inputnp;

    //Allocate memory
    Allocate(&temparray,np+1);//+1 for sentinel

    if(index=NULL)
    {
        Allocate(&index,np);//if index is not provided allocate memory
        //initiate index
        for(i=0; i<np; i++)index[i]=i;
    }


    if(shuffle)Shuffle();//in case they are not random

    start=1;//start from 1 first point is already considered in HPCO order
    end=3;//first time we sort 3 points

    while(1)
    {
        sortCPOyx(start,end);//sort bucket in cpo order
        if(end==np-1)break;//bail out
        start=end+1;
        current_nt=(end*2+1);
        end=start+nt_factor*current_nt;//k*current number of traingles with the given set of points
        if(start==end)end=start+1;//avoid this condition
        nextstep=end+1+nt_factor*(end*2+1);//compute nextstep to predict overflow
        if (nextstep>np-1)end=np-1;
    }


    Reindex();//finally sort the parray

    //Memory deallocation
    Deallocate(&temparray);//we deallocate only temparray.Index may be usefull for reindex outside the class
}



void HCPO3D::sortHCPO3D(Coord3D*inputp,int inputnp,bool shuffle)//sort points in HCPO order
{
    int start,end,nextstep,i,current_nt;

    //Assign variables
    p=inputp;
    np=inputnp;

    //Allocate memory
    Allocate(&temparray,np+1);//+1 for sentinel

    if(index==NULL)
    {
        Allocate(&index,np);//if index is not provided allocate memory
        //initiate index
        for(i=0; i<np; i++)index[i]=i;
    }


    if(shuffle)Shuffle();//in case they are not random

    start=1;//start from 1 first point is already considered in HPCO order
    end=4;//first time we sort 3 points

    while(1)
    {
        sortCPOzyx(start,end);//sort bucket in cpo order
        if(end==np-1)break;//bail out
        start=end+1;
        current_nt=(end*3+1);//predicted numbers of tetraedrons
        end=start+nt_factor*current_nt;//k*current number of traingles with the given set of points
        if(start==end)end=start+1;//avoid this condition
        nextstep=end+1+nt_factor*(end*3+1);//compute nextstep to predict overflow
        if (nextstep>np-1)end=np-1;
    }


    Reindex();//finally sort the parray

    //Memory deallocation
    Deallocate(&temparray);//we deallocate only temparray.Index may be usefull for reindex outside the class
}



void HCPO3D::Reindex()//sorts the p array according to indexes
{
    int i;
    Coord3D* newp=NULL;
    Allocate(&newp,np);//buils a new aarray of points

    Copy(newp,p,np);//copy contents

    //reindex
    for(i=0; i<np; i++)p[i]=newp[index[i]];

    Deallocate(&newp);//free memory

}
void HCPO3D::Shuffle()//Shuffle points order
{
    for (int i=0; i<np; i++)
    {
        int r = i + (rand() % (np-i)); // Random remaining position.
        int temp = index[i];
        index[i] = index[r];
        index[r] = temp;
    }

}
bool HCPO3D::IsSorted(double* v,int N)//check wheter the array is dorted or not
{
    int i;
    for(i=0; i<N-1; i++)
    {
        if (v[i]>v[i+1])
        {
            return false;
        }
    }

    return true;

}
