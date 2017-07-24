//Library that uses a PowerCrust based algorithm for surface recon

#ifndef _ROBUSTCRUST_H_
#define _ROBUSTCRUST_H_

#ifndef M_PI
#define M_PI 3.14159265358979323846 //pigreco
#endif

/*speed up debug mode*/
#define _ITERATOR_DEBUG_LEVEL 0

//INCLUDES
#include <iostream>
#include "structs/Triangle.h"
#include "structs/Coord3D.h"
#include "sortlib/HCPO3D.h"
#include "triangulations/TrianglesFunctions.h"
#include "util/ArraysLib.h"
#include "util/filemanager.h"
#include "tetgen1.5.0/predicates.cxx"
#include "tetgen1.5.0/tetgen.cxx"
#include <queue>
#include <stack>


static const int TET_EDGE_NODE[6][2] = {
   { 0, 1 },{ 1, 2 },{ 2, 0 },{ 0, 3 },{ 1, 3 },{ 2, 3 }
};

/* also indicates edge rotation */
static const int TET_EDGE_FACE[6][2] = {
    { 2, 3 },{ 0, 3 },{ 1, 3 },{ 1, 2 },{ 2, 0 },{ 0, 1 }
};

static const int TET_FACE[4][3] = {
    { 1, 2, 3 },{ 2, 0, 3 },{ 0 ,1, 3 },{ 0, 2, 1 }
};


//STRUCTURES
struct Tetrahedra
{
    int p1;
    int p2;
    int p3;
    int p4;
};

struct Neigh
{
    int T1;
    int T2;
    int T3;
    int T4;
};
struct Model3D
{

    double mx;
    double Mx;
    double my;
    double My;
    double mz;
    double Mz;
    Coord3D c;//center
    double r;// half the distance of the maximu enclosing cuboid edge
};

struct PQueueNode
{
    float Ifact;//Ifactor
    float plevel;//priority level
    int T1;
    int T2;
    bool operator <(const PQueueNode& rhs)
    {
        return plevel < rhs.plevel;
    }
};

//operaor to compare 2 priority queue nodes
class CompareNode
{
    public:
        bool operator()(PQueueNode& t1, PQueueNode& t2)
        {
            return t1.plevel < t2.plevel;
        }
};

enum PC_TET_STATUS
{
    PC_UNKNOWN,
    PC_INSIDE,
    PC_OUTSIDE
};

class RCRUST
{
        //GLOBALS

    public:

        vector<Triangle> t;
        FILE_MANAGER FileManager;

        //public functions
        int Triangulate(vector<double>* inputp, int innp);
        void FreeMemory();
        void CopyTriangles(int* triangles);
        RCRUST();
        ~RCRUST();

    private:
        tetgenio in, out;//for delaunnay triangulation
        Coord3D* p;
        int N;//number of cloud points
        int Nshield;//number of shiedl points
        Tetrahedra* T;
        Neigh* Tneigh;
        int NT;
        int ct;//triangles iterator
        Model3D Model;//infos about the model
        priority_queue<PQueueNode, vector<PQueueNode>, CompareNode > queue;
        PC_TET_STATUS* TetStatus;

        //private functions
        inline double Ifact(Coord3D* cc1, double* r1, Coord3D* cc2, double* r2);
        void TetraCC(int idT, Coord3D* cc, double* r);
        void AnalyzeModel();
        void Marking();
        void BuildSurface();
        void InitQueue();
        double BallCenter(Coord3D* cc, double* cr2, Coord3D* tnorm, Coord3D* BC);
        void AnalyzeIntersection(int T1, int T2);
        void OutsideTetraCC(Tetrahedra* T1, int T2, Coord3D* cc2, double* r2);
        float ComputePlevel(float Ifact);
        void MarkTetra(int T1, int T2, float Ifact);
        int DelaunayTriangulation();
        void HandleError(int flag);
        void FindFacet(int T1, int T2, Triangle* facet);
        void GenerateTriangles();
        void LookForTriangle( int T1, int T2);
        void Clustering();
        void ManifoldExtraction();
        double TetraVol(int idT);
        void CopyPoints(vector<double>* inputp, int innp);
        void AddShield(Coord3D* p, int N);
};

RCRUST::RCRUST()//constructor
{
    ct = 0;
    TetStatus = NULL;
    T = NULL;
    Tneigh = NULL;
    p = NULL;
}

RCRUST::~RCRUST()//desstructor
{
    //deallocate memory
    FreeMemory();
}

void RCRUST::FreeMemory()//desstructor
{
    //deallocate memory

    Deallocate(&TetStatus);

    //free the memory from delaunay triangulation
    out.deinitialize();//delete delaunay triangulation
    // Deallocate(&T);
    // Deallocate(&p); //bisogna rmuovere gli shield points
    t.clear();

}


int RCRUST::Triangulate(vector<double>* inputp, int innp)
{
    int exitcode;

    CopyPoints(inputp, innp); //copy pointers

    AnalyzeModel();

    AddShield(&p[N], Nshield);

    exitcode = DelaunayTriangulation();

    BuildSurface();
    //out.deinitialize();//delete delaunay triangulation

    return 0;
}

void RCRUST::CopyPoints(vector<double>* inputp, int innp)
{
    int i;
    N = innp;
    Nshield = N * .01;

    //if points are just a few we still have to use a few shield points
    if (N < 10000)
    {
        Nshield = 100;    //Use at least 100 shield points
    }
    if (N < 1000)
    {
        Nshield = 50;    //Use at least 50 shield points
    }

    //Adding space for shield points
    inputp->reserve((N + Nshield) * 3);

    //passing the pointer
    p = (Coord3D*) &inputp->front();

}


void RCRUST::AddShield(Coord3D* ps, int n)
{
    /*    1.  Choose z uniformly distributed in [−1,1].
    2. Choose t uniformly distributed on [0, 2 π).
    3. Let r = √(1−z2).
    4. Let x = r * cos(t).
    5. Let y = r * sin(t).
    */
    int i;//,flag;
    double* u = NULL;
    Allocate(&u, n);
    double* t = NULL;
    Allocate(&t, n);
    double* r = NULL;
    Allocate(&r, n);

    //Seed the random generatom to have always the same random numbers

    Random(u, n, (double) 0, (double) 1);
    Random(t, n, (double) 0, (double)2 * M_PI );

    //z equals u;
    for (i = 0; i < n; i++)
    {
        ps[i].z = u[i];
    }
    //r = √(1−z2);
    for (i = 0; i < n; i++)
    {
        r[i] = sqrt(1 - u[i] * u[i]);
    }
    // x and y
    for (i = 0; i < n; i++)
    {
        ps[i].x = r[i] * cos(t[i]);
        ps[i].y = r[i] * sin(t[i]);
    }

    //setting center and radius
    for (i = 0; i < n; i++)
    {
        ps[i].x = ps[i].x * Model.r + Model.c.x;
        ps[i].y = ps[i].y * Model.r + Model.c.y;
        ps[i].z = ps[i].z * Model.r + Model.c.z;

    }


    Deallocate(&u);
    Deallocate(&r);
    Deallocate(&t);
}

int RCRUST::DelaunayTriangulation()
{
    //compute the delaunay triagnulariotn of the dataset

    // HCPO3D hpco;

    //HPCO sort

    //hpco.sortHCPO3D(p,N,true);//pre shuflle and sort points
    //Launch tessellation
    in.firstnumber = 0; //first number in labeling
    in.pointlist = &p[0].x; //assign points list
    in.numberofpoints = N + Nshield;

    tetgenbehavior behav;//behaviour class for tetgen
    //behav.insertaddpoints=true;
    //behav.plc=true;
    //behav.conformdel = true;
    behav.nonodewritten = true;
    behav.noelewritten = false;
    behav.nofacewritten = true;
    behav.edgesout = false;
    behav.facesout = false;
    behav.neighout = true;
    behav.voroout = false;

    tetrahedralize(&behav, &in, &out);

    //copying pointers to local data structure

    T = (Tetrahedra*)out.tetrahedronlist;
    Tneigh = (Neigh*)out.neighborlist;
    NT = out.numberoftetrahedra;

#if 0
    std::cout << "Nt= " << NT << std::endl;
    std::cout << "Maximum value= " << Max(&T[0].p1, NT * 4) << std::endl;
    std::cout << "Maximum value= " << Max(&Tneigh[0].T1, NT * 4) << std::endl;
#endif
    return 0;
}


void RCRUST::TetraCC(int idT, Coord3D* cc, double* r)
{
    //compute the circumcenter of a tetraedron


    //tetra points
    Coord3D p1 = p[T[idT].p1];
    Coord3D p2 = p[T[idT].p2];
    Coord3D p3 = p[T[idT].p3];
    Coord3D p4 = p[T[idT].p4];

    //tetra vectors
    Coord3D v12, v13, v14;
    DiffPoints(p1., p2., v12.)
    DiffPoints(p1., p3., v13.)
    DiffPoints(p1., p4., v14.)

    //midpoints
    Coord3D m12, m13, m14;
    MidPoint(p1., p2., m12.)
    MidPoint(p1., p3., m13.)
    MidPoint(p1., p4., m14.)

    //coefficent vector
    double d[3];
    DotProduct(v12., m12., d[0])
    DotProduct(v13., m13., d[1])
    DotProduct(v14., m14., d[2])

    //Solve the system
    FastCramer(&v12.x, &v13.x, &v14.x, d, &cc->x); //circumcenter

    Distance(cc->, p1., *r) //circumradius



}


inline double RCRUST::Ifact(Coord3D* cc1, double* r1, Coord3D* cc2, double* r2)
{
    //compute intersection factor between 2 tetrahedrons

    double distcc;
    double Ifact;

    SquaredDistance(cc1->, cc2->, distcc);

    Ifact = (*r1** r1 + *r2** r2 - distcc) / (2 * *r1** r2); //Carnot teorem
#if 0
    if (Ifact >= 1 || Ifact <= -1)
    {
        Ifact = Ifact; //wrong Ifactor
    }
#endif

    return Ifact;

}

void RCRUST::AnalyzeModel()
{
    //Mark Tetraedrons as inside or outside

    //model bounding box
    Model.mx = Min(&p[0].x, N, 3);
    Model.Mx = Max(&p[0].x, N, 3);
    Model.my = Min(&p[0].y, N, 3);
    Model.My = Max(&p[0].y, N, 3);
    Model.mz = Min(&p[0].z, N, 3);
    Model.Mz = Max(&p[0].z, N, 3);

    //model center
    Model.c.x = (Model.Mx + Model.mx) / 2;
    Model.c.y = (Model.My + Model.my) / 2;
    Model.c.z = (Model.Mz + Model.mz) / 2;


    //model radius (look for the greatest dx dy d
    double dx, dy, dz;
    dx = Model.c.x - Model.Mx;
    dy = Model.c.y - Model.My;
    dz = Model.c.z - Model.Mz;

    Model.r = 2 * sqrt(dx * dx + dy * dy + dz * dz); //radius =radius*2

}


void RCRUST::InitQueue()
{
    //scan all convex hull facets and initiate his tetraedrons as inside or outside

    int i;
    // Coord3D cc;//double r;
    // PQueueNode Node;


    //FLag as checked and outside all shield tetraedrons
    for (i = 0; i < NT; i++)
    {
        if (T[i].p1 >= N || T[i].p2 >= N || T[i].p3 >= N || T[i].p4 >= N)
        {
            TetStatus[i] = PC_OUTSIDE;
        }
        else
        {
            TetStatus[i] = PC_UNKNOWN;
        }
    }
    //add to queue all neighbours of shield tetraedrons

    for (i = 0; i < NT; i++)
    {
        if (TetStatus[i] == PC_OUTSIDE)
        {
            if (Tneigh[i].T1 >= 0)
            {
                AnalyzeIntersection(i, Tneigh[i].T1);
            }
            if (Tneigh[i].T2 >= 0)
            {
                AnalyzeIntersection(i, Tneigh[i].T2);
            }
            if (Tneigh[i].T3 >= 0)
            {
                AnalyzeIntersection(i, Tneigh[i].T3);
            }
            if (Tneigh[i].T4 >= 0)
            {
                AnalyzeIntersection(i, Tneigh[i].T4);
            }
        }

    }

}


void RCRUST::AnalyzeIntersection(int T1, int T2)
{
    // Analyze the intersection between two tetrahedra. Data are stored in the priority queue
    // T1,T2 is always >0; The marking flow goes from T1 to T2
    Coord3D cc1, cc2;
    double r1, r2;

    //stop for debug
#if 0
    if (T1 == -7581 && T2 == 12)
    {
        T1 = T1;
    }
#endif

    PQueueNode Node;

    if (TetStatus[T2] != PC_UNKNOWN)
    {
        return;    //no need to check go forward
    }

    TetraCC(T1, &cc1, &r1);
    TetraCC(T2, &cc2, &r2);

    Node.Ifact = Ifact(&cc1, &r1, &cc2, &r2);
    Node.T1 = T1;
    Node.T2 = T2;
    Node.plevel = ComputePlevel(Node.Ifact);

    queue.push(Node);//add to queue

}

float RCRUST::ComputePlevel(float Ifact)
{
    float plevel;

    if (Ifact < 0)
    {
        plevel = -Ifact;
    }
    else
    {
        plevel = Ifact;
    }

    return plevel;
}

void RCRUST::FindFacet(int T1, int T2, Triangle* facet)
{
    //retunrs the facet shared by T1 and T2, T1 is a not deleted tetraedron
    //getting the position of T2
    int p1, p2, p3;
    int face_id;
    int *T1ptr = &T[T1].p1;

    //    Coord3D tnorm;
    if (Tneigh[T1].T1 == T2)
    {
        face_id = 0;
    }
    else if (Tneigh[T1].T2 == T2)
    {
        face_id = 1;
    }
    else if (Tneigh[T1].T3 == T2)
    {
        face_id = 2;
    }
    else if (Tneigh[T1].T4 == T2)
    {
        face_id = 3;
    }
    else
    {
        Error("Facet error");
    }

     p1 = TET_FACE[face_id][0];
     p2 = TET_FACE[face_id][1];
     p3 = TET_FACE[face_id][2];

     facet->p1 = T1ptr[p1];
     facet->p2 = T1ptr[p2];
     facet->p3 = T1ptr[p3];

}


void RCRUST::Marking()
{
    //mark all tetraedrons as inside or outside
    PQueueNode Node;
    int c = 0;
    while (!queue.empty())
    {
        Node = queue.top();
        queue.pop();//get and delete the last element form the queue

        if (Node.T2 < 0 || TetStatus[Node.T2] != PC_UNKNOWN)
        {
            continue;    //no need to check outside tetra our already checked
        }

        MarkTetra(Node.T1, Node.T2, Node.Ifact);

    }

}
void RCRUST::MarkTetra(int T1, int T2, float Ifact)
{
    //Mark the tetraedrons T2 Accoording to the Ifact and the status of T1
    //T2 is always >=0, T1<0 is assumed as deleted
    //after marking analyze the connection with his neighbours
    PC_TET_STATUS status;

    if (TetStatus[T2] != PC_UNKNOWN)
    {
        return;    //no need to recheck
    }

    if (T1 >= 0)
    {
        status = TetStatus[T1];
    }
    else
    {
        status = PC_OUTSIDE; //outside tetraedrons are deleted
    }



    if (Ifact > 0) //intersection marl as equal
    {
        TetStatus[T2] = status;
    }
    else
    {
        if (TetStatus[T1] == PC_INSIDE)
        {
            TetStatus[T2] = PC_OUTSIDE;
        }
        else
        {
            TetStatus[T2] = PC_INSIDE;
        }
    }

    //analyize connection with his neighbours
    if (Tneigh[T2].T1 >= 0)
    {
        AnalyzeIntersection(T2, Tneigh[T2].T1);
    }
    if (Tneigh[T2].T2 >= 0)
    {
        AnalyzeIntersection(T2, Tneigh[T2].T2);
    }
    if (Tneigh[T2].T3 >= 0)
    {
        AnalyzeIntersection(T2, Tneigh[T2].T3);
    }
    if (Tneigh[T2].T4 >= 0)
    {
        AnalyzeIntersection(T2, Tneigh[T2].T4);
    }

    return;
}


void RCRUST::GenerateTriangles()
{
    //finds the intersection between deleted and not tetraedra

    for(int i = 0; i< NT; i++)
    {
        LookForTriangle( i, Tneigh[i].T1);
        LookForTriangle( i, Tneigh[i].T2);
        LookForTriangle( i, Tneigh[i].T3);
        LookForTriangle( i, Tneigh[i].T4);
    }

}

/* Generates a triangle between T1 INSIDE and T2 OUTSIDE*/
void RCRUST::LookForTriangle(int T1, int T2)
{
    Triangle facet;
    bool deleted2;

    if ( T1 >= 0 &&
         T2 >= 0 &&
         TetStatus[T1] == PC_INSIDE &&
         TetStatus[T2] == PC_OUTSIDE)
    {
        FindFacet(T1, T2, &facet);

       if (facet.p1 < N && facet.p2 < N && facet.p3 < N)
       {
           t.push_back(facet);
       }//new triangle is generated}
       else
       {
           cout << "Erroneus facet " << T1 << " " << T2 << endl;
       }
    }

}

void RCRUST::BuildSurface()
{
    //All operation to extract the surface
    //int flag;
    AllocateAndInit(&TetStatus, NT, PC_UNKNOWN); //true for checked tetraedroms

    InitQueue();

    Marking();

    Clustering();

    GenerateTriangles();

    cout << "Genrated: " << t.size() << " triangles From " << NT << " Tetraedrons" << endl;

    //remove from memory
    Deallocate(&TetStatus);
}


void RCRUST::CopyTriangles(int* triangles)
{
    int i = 0;
    int c = 0;
    for (i = 0; i < t.size(); i++)
    {
        triangles[c] = t[i].p1;
        triangles[c + 1] = t[i].p2;
        triangles[c + 2] = t[i].p3;
        c = c + 3;
    }

}

void RCRUST::Clustering()
{
    //set  to not deleted the biggest cluster of tetraedrons

    int i;
    vector<int> Stack;
    int* NumTetra = NULL; //number of tetraedrons inside a cluster
    int* ClustId = NULL;
    AllocateAndInit(&ClustId, NT, -1);
    int idClust;// id of the current cluster
    int BiggerClust;
    int T1, T2;

    idClust = 0;
    for(i = 0; i < NT; i++)
    {
        T1 = i;

        if (ClustId[T1] >= 0)
            continue;

        //Start from tetra 0
        Stack.push_back(T1);
        while (!Stack.empty())
        {
            T1 = Stack.back();
            Stack.pop_back();//get and delete element
            ClustId[T1] = idClust;
            T2 = Tneigh[T1].T1;
            if (T2 >= 0 && 
                TetStatus[T2] == PC_INSIDE &&
                ClustId[T2] == -1)
            {
                Stack.push_back(T2);
            }
            T2 = Tneigh[T1].T2;
            if (T2 >= 0 &&
                TetStatus[T2] == PC_INSIDE &&
                ClustId[T2] == -1)
            {
                Stack.push_back(T2);
            }
            T2 = Tneigh[T1].T3;
            if (T2 >= 0 &&
                TetStatus[T2] == PC_INSIDE &&
                ClustId[T2] == -1)
            {
                Stack.push_back(T2);
            }
            T2 = Tneigh[T1].T4;
            if (T2 >= 0 &&
                TetStatus[T2] == PC_INSIDE &&
                ClustId[T2] == -1)
            {
                Stack.push_back(T2);
            }
        }

        idClust++;//go to next cluster
    }

    //counting the number of tetraedrons inside each cluster
    AllocateAndInit(&NumTetra, idClust, 0);
    for (i = 0; i < NT; i++)
    {
        if (ClustId[i] >= 0)
        {
            NumTetra[ClustId[i]]++;
        }
    }
    //finding the biggerCLuster
    BiggerClust = 0;
    for (i = 1; i < idClust; i++)
    {
        if (NumTetra[i] > NumTetra[BiggerClust])
        {
            BiggerClust = i;
        }
    }

    cout << "Found " << idClust << "Clusters. Bigger cluster had: " << NumTetra[BiggerClust] << " Tetraedrons" << endl;
    //flagging as deleted all tetraedrons non belonging to the bigger cluster
    for (i = 0; i < NT; i++)
    {
        if (ClustId[i] == BiggerClust)
        {
            TetStatus[i] = PC_INSIDE;
        }
        else
        {
            TetStatus[i] = PC_OUTSIDE;
        }
    }

}

void RCRUST::ManifoldExtraction()
{
    //Extract the Manifold performing  a walking operation
}

double RCRUST::TetraVol(int idT)
{
    //Gets the volume of the tetraedron idt

    //tetra points
    Coord3D p1 = p[T[idT].p1];
    Coord3D p2 = p[T[idT].p2];
    Coord3D p3 = p[T[idT].p3];
    Coord3D p4 = p[T[idT].p4];

    //tetra vectors
    Coord3D v12, v13, v14;
    DiffPoints(p1., p2., v12.)
    DiffPoints(p1., p3., v13.)
    DiffPoints(p1., p4., v14.)

    double* v1, *v2, *v3;
    v1 = &v12.x;
    v2 = &v13.x;
    v3 = &v14.x;

    double det = v1[0] * v2[1] * v3[2] - v2[2] * v3[1] +
                 v1[1] * v2[2] * v3[0] - v2[0] * v3[2] +
                 v1[2] * v2[0] * v3[1] - v2[1] * v3[0];

    return det / 6;
}

#endif
