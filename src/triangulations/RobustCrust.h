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
#include <set>
#include "csparse/cs.h"

static const int TET_EDGE_NODE[6][2] = {
   { 0, 1 },{ 1, 2 },{ 2, 0 },{ 0, 3 },{ 1, 3 },{ 2, 3 }
};

/* also indicates edge rotation with tets*/
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
    int operator [](int idx) {
        int *ptr = &this->p1;
        return ptr[idx];
    }
};

struct Neigh
{
    int T1;
    int T2;
    int T3;
    int T4;
    int operator [](int idx) {
        int *ptr = &this->T1;
        return ptr[idx];
    }
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

struct ABPAQueueNode
{
    int T_in;/*inside*/
    int T_out;
    double plevel;
    int p1;
    int p2;

    bool operator<(const ABPAQueueNode &rhs) const 
    {
        if (T_in != rhs.T_in)
        {
            return T_in < rhs.T_in;
        }
        else
        {
            return T_out < rhs.T_out;
        }
    }

    bool operator>(const ABPAQueueNode &rhs) const
    {
        if (T_in != rhs.T_in)
        {
            return T_in > rhs.T_in;
        }
        else
        {
            return T_out > rhs.T_out;
        }
    }
};

//operator to compare 2 priority queue nodes
class CompareABPANode
{
public:
    bool operator()(ABPAQueueNode& t1, ABPAQueueNode& t2)
    {
        /*reverting < as we want negative ifact high priority*/
        return t1.plevel > t2.plevel;
    }
};

typedef priority_queue<PQueueNode, vector<PQueueNode>, CompareNode > PC_Queue;
typedef priority_queue<ABPAQueueNode, vector<ABPAQueueNode>, CompareABPANode > ABPA_Queue;

class RCRUST
{
        //GLOBALS

    public:

        vector<Triangle> t;
        FILE_MANAGER FileManager;

        //public functions
        int TriangulatePowerCrust(vector<double>* inputp, int innp);
        int TriangulateABPA(vector<double>* inputp, int innp);
        int RCRUST::TriangulateSpectral(vector<double>* inputp, int innp);
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
        PC_Queue queue;
        PC_TET_STATUS* TetStatus;

        //private functions
        inline double IfactFromCC(Coord3D* cc1, double* r1, Coord3D* cc2, double* r2);
        inline double IfactFromTetra(int T1, int T2);
        bool ComputeTetraCC(int idT, Coord3D* cc, double* r);
        void AnalyzeModel();
        void Marking();
        void BuildSurface();
        void InitQueue();
        double BallCenter(Coord3D* cc, double* cr2, Coord3D* tnorm, Coord3D* BC);
        void AnalyzeIntersection(int T1, int T2);
        float ComputePlevel(float Ifact);
        void MarkTetra(int T1, int T2, float Ifact);
        int DelaunayTriangulation();
        void HandleError(int flag);
        void FindFacet(int T1, int T2, Triangle* facet);
        void GenerateTriangles();
        void LookForTriangle( int T1, int T2);
        void Clustering();
        double TetraVol(int idT);
        void CopyPoints(vector<double>* inputp, int innp);
        void AddShield(Coord3D* p, int N);

        /*ABPA*/
        void ABPA();
        ABPAQueueNode RCRUST::ABPA_get_seed();
        bool RCRUST::TetraIsInfinite(int *tetra);
        int RCRUST::TetEdgeRotationForward(int Tid, int p0, int p1);
        int RCRUST::TetEdgeRotationBackward(int Tid, int p0, int p1);
        int RCRUST::TetEdgeMapOnTetra(int* tet, int p0, int p1);
        void RCRUST::ABPA_propagate(ABPAQueueNode *node, ABPA_Queue *queue);
        ABPAQueueNode RCRUST::RunPivotOnEdge(int Tin, int Tout, int p1, int p2, bool mark);
        int RCRUST::CollectPivotingTetras(int Tin, int Tout, int p1, int p2, vector<int> *forward);
        ABPAQueueNode RCRUST::FindMinIfactPivot(vector<int> &tetras, int Tin, bool mark);
        void RCRUST::SetTetStatusAndUpdateQueue(int Tid, PC_TET_STATUS status, ABPA_Queue *queue);
        inline double RCRUST::IfactFromTetraFaceCC(int T1, int T2);
        void RCRUST::FaceCircumCenter(Coord3D* p1, Coord3D* p2, Coord3D* p4, Coord3D* n, double* r, Coord3D* cc);

        /*SPECTRAL*/
        void SpectralPropagation();

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


bool RCRUST::ComputeTetraCC(int idT, Coord3D* cc, double* r)
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
    if (FastCramer(&v12.x, &v13.x, &v14.x, d, &cc->x))
    {
        Distance(cc->, p1., *r) //circumradius
        return true;
    }

    return false;

}

bool RCRUST::TetraIsInfinite(int *tetra)
{
    for (int i = 0; i < 4; i++)
    {
        if (tetra[i] >= N)
            return true;
    }
    return false;
}

inline double RCRUST::IfactFromCC(Coord3D* cc1, double* r1, Coord3D* cc2, double* r2)
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

// T1,T2 is always >0; The marking flow goes from T1 to T2
inline double RCRUST::IfactFromTetra(int T1, int T2)
{
    Coord3D cc1, cc2;
    double r1, r2;
    double ifact;

    if (!ComputeTetraCC(T1, &cc1, &r1))
        return 0.0;

    if (!ComputeTetraCC(T2, &cc2, &r2))
        return 0.0;

    ifact = IfactFromCC(&cc1, &r1, &cc2, &r2);

    return ifact;
}

inline double RCRUST::IfactFromTetraFaceCC(int T1, int T2)
{
    Coord3D cc, norm;
    double cr;
    Triangle facet;

    FindFacet(T1, T2, &facet);

    FaceCircumCenter( &p[facet.p1], &p[facet.p2], &p[facet.p3], &norm, &cr, &cc);

    return cr;
}

// Analyze the intersection between two tetrahedra. Data are stored in the priority queue
void RCRUST::AnalyzeIntersection(int T1, int T2)
{
    // T1,T2 is always >0; The marking flow goes from T1 to T2

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

    Node.Ifact = IfactFromTetra(T1, T2);
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
    int face_id = -1;

    for (int i = 0; i < 4; i++)
    {
        if (Tneigh[T1][i] == T2)
        {
            face_id = i;
            break;
        }
    }

    if(face_id == -1)
    {
        Error("Facet error");
    }

     p1 = TET_FACE[face_id][0];
     p2 = TET_FACE[face_id][1];
     p3 = TET_FACE[face_id][2];

     facet->p1 = T[T1][p1];
     facet->p2 = T[T1][p2];
     facet->p3 = T[T1][p3];

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

/*
 Triangualtion routine for the volumetric robust Crust
*/
int RCRUST::TriangulatePowerCrust(vector<double>* inputp, int innp)
{
    int exitcode;

    CopyPoints(inputp, innp); //copy pointers

    AnalyzeModel();

    AddShield(&p[N], Nshield);

    exitcode = DelaunayTriangulation();

    BuildSurface();

    return 0;
}


ABPAQueueNode RCRUST::ABPA_get_seed()
{
    double ifact, min_ifact = DBL_MAX;
    ABPAQueueNode node;

    for (int i = 0; i < N; i ++)
    {
        int *Tn = &Tneigh[i].T1;

        for (int j = 0; j < 4; j++)
        {
            if (i > Tn[j] &&
                    Tn[j] >= 0)
            {
                ifact = IfactFromTetra(i, Tn[j]);

                if (ifact < min_ifact)
                {
                    min_ifact = ifact;
                    node.T_in = i;
                    node.T_out = Tn[j];
                }
            }
        }
    }

    return node;
}


int RCRUST::TetEdgeMapOnTetra(int* tet, int p0, int p1)
{
    for(int i = 0; i< 6; i++)
    {
        int p0_tet_id = TET_EDGE_NODE[i][0];
        int p1_tet_id = TET_EDGE_NODE[i][1];

        if( tet[p0_tet_id] == p0 &&
            tet[p1_tet_id] == p1)
        {
            return i;
        }
        else if( tet[p0_tet_id] == p1 &&
                 tet[p1_tet_id] == p0)
        {
            return i;
        }
    }

    Error(" Cannot map edge on Tetra\n");
    return 0;
}

int RCRUST::TetEdgeRotationForward(int Tid, int Tprev, int tet_edge_id)
{
    int id = TET_EDGE_FACE[tet_edge_id][0];
    int Tn_id = Tneigh[Tid][id];

    if( Tneigh[Tid][id] != Tprev)
    {
        return Tn_id;
    }
    else
    {
        id = TET_EDGE_FACE[tet_edge_id][1];
        Tn_id = Tneigh[Tid][id];
        return Tn_id;
    }
}

int RCRUST::TetEdgeRotationBackward(int Tid, int p0, int p1)
{
    return 0;
}

void RCRUST::FaceCircumCenter(Coord3D* p1, Coord3D* p2, Coord3D* p4, Coord3D* n, double* r, Coord3D* cc)
{
    Coord3D b, m;
    double leng;
    double rtemp;
    Coord3D v21, v41;

    //Getting the normal
    DiffPoints(p2->, p1->, v21.);
    DiffPoints(p4->, p1->, v41.);
    CrossProduct(v41., v21., n->);
    Normalize(n->, leng);

    //Assembling the system

    DotProduct(n->, p1->, b.x);

    MidPoint(p1->, p2->, m.);

    DotProduct(v21., m., b.y);

    MidPoint(p1->, p4->, m.);

    DotProduct(v41., m., b.z);

    //Solve the system
    FastCramer(&n->x, &v21.x, &v41.x, &b.x, &cc->x);

    //CircumRadius
    SquaredDistance(cc->, p1->, *r);
    //take the minum distance
    //avoids the triangle to self destruct during Dtest
    SquaredDistance(cc->, p2->, rtemp);
    if (rtemp < *r)
    {
        *r = rtemp;
    }
    SquaredDistance(cc->, p4->, rtemp);
    if (rtemp < *r)
    {
        *r = rtemp;
    }

#ifdef _DEBUG
    double d1, d2, d3;
    Distance(cc->, p1->, d1);
    Distance(cc->, p2->, d2);
    Distance(cc->, p4->, d3);
#endif


}

ABPAQueueNode RCRUST::FindMinIfactPivot(vector<int> &tetras, int Tin, bool mark)
{
    double min_ifact = DBL_MAX, ifact;
    int pivot = -1;
    ABPAQueueNode node;

    node.T_in = -1;
    node.T_out = -1;

    for (int i = 0; i < tetras.size(); i++)
    {
        int Tid1 = tetras[i];
        int Tid2;
        if (i < tetras.size() - 1)
        {
            Tid2 = tetras[i + 1];
        }
        else
        {
            Tid2 = Tin;
        }

#if 1
        ifact = IfactFromTetra(Tid1, Tid2);
#else
        ifact = IfactFromTetraFaceCC(Tid1, Tid2);
#endif

        if (ifact < min_ifact)
        {
            pivot = i;
            min_ifact = ifact;
            node.plevel = ifact;
            node.T_out = Tid1;
            node.T_in = Tid2;
        }
    }

    if (mark)
    {
        for (int i = 0; i < tetras.size(); i++)
        {
            int Tid = tetras[i];
            if (i <= pivot)
            {
                SetTetStatusAndUpdateQueue(Tid, PC_OUTSIDE, NULL);
            }
            else
            {
                SetTetStatusAndUpdateQueue(Tid, PC_INSIDE, NULL);
            }
        }
    }

    return node;
}

int RCRUST::CollectPivotingTetras(int Tin, int Tout, int p1, int p2, vector<int> *forward)
{
    int edge_id_on_current;
    int nextT, currentT,prevT;

    assert(TetStatus[Tin] == PC_INSIDE);
    assert(TetStatus[Tout] == PC_OUTSIDE);

    forward->push_back(Tout);
    prevT = Tin;
    currentT = Tout;
    while(1)
    {
        edge_id_on_current = TetEdgeMapOnTetra((int*)&T[currentT], p1, p2);
        nextT = TetEdgeRotationForward(currentT, prevT, edge_id_on_current);

        if (nextT < 0)
            break;

        /*if meet one inside bail out*/
        /*if meet one outside bail out*/
        if (TetStatus[nextT] == PC_INSIDE)
            break;

        forward->push_back(nextT);
        prevT = currentT;
        currentT = nextT;
    }

    return nextT;
}


void RCRUST::SetTetStatusAndUpdateQueue(int Tid, PC_TET_STATUS status, ABPA_Queue *queue)
{
    ABPAQueueNode node;

    if (TetStatus[Tid] == status)
        return; /* nothing changed*/

    if(TetStatus[Tid] != PC_INSIDE)
    {
        TetStatus[Tid] = status;
        for (int i = 0; i < 4; i++)
        {
            int neighT = Tneigh[Tid][i];
            if (neighT >= 0 &&
                TetStatus[neighT] != PC_UNKNOWN &&
                TetStatus[neighT] != status)
            {
                if (status == PC_INSIDE)
                {
                    node.T_in = Tid;
                    node.T_out = neighT;
                }
                else
                {
                    node.T_out = Tid;
                    node.T_in = neighT;
                }
                if( queue != NULL)
                {
                    node.plevel = IfactFromTetra(Tid, neighT);
                    queue->push(node);
                    assert(TetStatus[node.T_in] == PC_INSIDE);
                    assert(TetStatus[node.T_out] == PC_OUTSIDE);
                }
            }
        }
    }
    else
    {
        Error("Trying to  change inside tet\n");
    }
}

ABPAQueueNode RCRUST::RunPivotOnEdge(int Tin, int Tout, int p1, int p2, bool mark)
{
    /*find edge and pivot orientation*/
    int tet_out_edge_id = TetEdgeMapOnTetra((int*)&T[Tout], p1, p2);
    int lastT;
    vector<int> tetras;
    ABPAQueueNode node;

    assert(TetStatus[Tin] == PC_INSIDE);
    assert(TetStatus[Tout] == PC_OUTSIDE);

    /*to invalidate the node*/
    node.T_out = Tout;
    node.T_in = Tin;

    /*collect pivot tetras and pivot status*/
    lastT = CollectPivotingTetras(Tin, Tout, p1, p2, &tetras);

    /*mark tetras depending on the pivot status*/
    if (lastT < 0 || TetStatus[lastT] == PC_OUTSIDE)
    {
#if 0
        Error("out of hull");
        //met convex hull
        for (int i = 0; i < tetras.size(); i++)
        {
            int Tid = tetras[i];

            SetTetStatusAndUpdateQueue(Tid, PC_OUTSIDE, NULL);
        }
#endif
    }
    else if (TetStatus[lastT] == PC_INSIDE)
    {
        node = FindMinIfactPivot(tetras, lastT, mark);
    }
    else
    {
        Error(" Unsupported Pivot");
    }

    return node;
}


void RCRUST::ABPA_propagate(ABPAQueueNode *node, ABPA_Queue *abpa_queue)
{
    Triangle facet;
    ABPAQueueNode new_node;


    assert(TetStatus[node->T_in] == PC_INSIDE);
    assert(TetStatus[node->T_out] == PC_OUTSIDE);

    FindFacet(node->T_in, node->T_out, &facet);

    /*pivot edges (orientation of edge?)*/
    node->p1 = facet.p1;
    node->p2 = facet.p2;
    new_node = RunPivotOnEdge(node->T_in, node->T_out, node->p1, node->p2, false);
    node->plevel = new_node.plevel;
    /* does the pivot change the surface?*/
    if (TetStatus[new_node.T_in] != PC_INSIDE ||
        TetStatus[new_node.T_out] != PC_OUTSIDE)
    {
        abpa_queue->push(*node);
    }


    node->p1 = facet.p2;
    node->p2 = facet.p3;
    new_node = RunPivotOnEdge(node->T_in, node->T_out, node->p1, node->p2, false);
    node->plevel = new_node.plevel;
    /* does the pivot change the surface?*/
    if (TetStatus[new_node.T_in] != PC_INSIDE ||
        TetStatus[new_node.T_out] != PC_OUTSIDE)
    {
        abpa_queue->push(*node);
    }

    node->p1 = facet.p3;
    node->p2 = facet.p1;
    new_node = RunPivotOnEdge(node->T_in, node->T_out, node->p1, node->p2, false);
    node->plevel = new_node.plevel;
    /* does the pivot change the surface?*/
    if (TetStatus[new_node.T_in] != PC_INSIDE ||
        TetStatus[new_node.T_out] != PC_OUTSIDE)
    {
        abpa_queue->push(*node);
    }

}

void RCRUST::ABPA()
{
    ABPAQueueNode node;
    ABPA_Queue abpa_queue;
    int count = 0;


    AllocateAndInit(&TetStatus, NT, PC_UNKNOWN); //true for checked tetraedroms

    node = ABPA_get_seed();

    SetTetStatusAndUpdateQueue(node.T_in, PC_INSIDE, NULL);
    SetTetStatusAndUpdateQueue(node.T_out, PC_OUTSIDE, NULL);

    ABPA_propagate(&node, &abpa_queue);

    while (!abpa_queue.empty())
    {
        //propagate//
        node = abpa_queue.top();
        abpa_queue.pop();//get and delete the last element form the queue

        /*cout << "priority "<< node.plevel  << endl;*/
#if 1
        if (node.plevel > -0.5)
            continue;
#endif
        /* run pivot and do marking if node is still valid*/
        if (TetStatus[node.T_in] == PC_INSIDE &&
            TetStatus[node.T_out] == PC_OUTSIDE)
        {
            ABPAQueueNode new_node = RunPivotOnEdge(node.T_in, node.T_out, node.p1, node.p2, true);

            /* run pivot on the newly marked and add to queue*/
            if (TetStatus[new_node.T_in] == PC_INSIDE &&
                TetStatus[new_node.T_out] == PC_OUTSIDE)
            {
                ABPA_propagate(&new_node, &abpa_queue);
            }
        }

        count++;
    }

    /*extract surface between in out Tetras*/
    GenerateTriangles();

    /*free*/
    Deallocate(&TetStatus); //true for checked tetrahedrons
}

/*
Triangulation routine for the Adaptive Ball Pivoting Algorithm
*/
int RCRUST::TriangulateABPA(vector<double>* inputp, int innp)
{
    int exitcode;

    CopyPoints(inputp, innp); //copy pointers

    AnalyzeModel();

    AddShield(&p[N], Nshield);

    exitcode = DelaunayTriangulation();

    ABPA();

    return 0;

}

/*
Triangulation routine with spectral
*/
int RCRUST::TriangulateSpectral(vector<double>* inputp, int innp)
{
    int exitcode;

    CopyPoints(inputp, innp); //copy pointers

    AnalyzeModel();

    AddShield(&p[N], Nshield);

    exitcode = DelaunayTriangulation();

    SpectralPropagation();

    return 0;

}

void RCRUST::SpectralPropagation()
{
    double *w = NULL;
    ABPAQueueNode node;
    AllocateAndInit(&w, NT, 0.0);
    AllocateAndInit(&TetStatus, NT, PC_UNKNOWN); //true for checked tetraedroms

    node = ABPA_get_seed();

#if 0

    cs *T = cs_spalloc(NT, NT, NT*4, 1.0, 1);
    for (int i = 0; i < NT; i++)
    {
        double tot_ifact = 0.0;
        double new_w = 0.0;
        double ifact[4] = { 0.0 };

        cs_entry(T, i, i, 1.0);

        if (i == node.T_in)
        {
            w[i] = -1.0;
            continue;
        }

        if (i == node.T_out)
        {
            w[i] = 1.0;
            continue;
        }

        for (int j = 0; j < 4; j++)
        {
            int Tn_id = Tneigh[i][j];

            if (Tn_id >= 0)
            {
                ifact[j] = IfactFromTetra(i, Tn_id);
                tot_ifact += abs(ifact[j]);
            }
        }
        for (int j = 0; j < 4; j++)
        {
            int Tn_id = Tneigh[i][j];

            if (Tn_id >= 0)
            {
                if (Tn_id == node.T_in)
                {
                    w[i] += -ifact[j] / tot_ifact;
                }
                else if (Tn_id == node.T_out)
                {
                    w[i] += ifact[j] / tot_ifact;
                }
                else
                {
                    double val = ifact[j] / tot_ifact;
                    cs_entry(T, i, Tn_id, val);
                }
                
            }
        }
    }
    cs *A = cs_compress(T);
    int ok = cs_lusol(0, A, w, 1e-1);
    T = cs_spfree(T);
    T = cs_spfree(A);
#else

    double *ifact_cache = NULL;
    Allocate(&ifact_cache, NT*4); //true for checked tetraedroms

    for (int i = 0; i < NT; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            int Tn_id = Tneigh[i][j];

            if (Tn_id >= 0)
            {
                ifact_cache[i*4+j] = IfactFromTetra(i, Tn_id);
            }
        }
    }

    for (int iter = 0; iter < 1000; iter++)
    {
        double tot_w = 0.0;
        w[node.T_in] = -1.0;
        w[node.T_out] = 1.0;
        for (int i = 0; i < NT; i++)
        {
            double tot_ifact = 0.0;
            double new_w = 0.0;
            double ifact[4] = { 0.0 };

            if (i == node.T_in ||
                i == node.T_out)
                continue;

            for (int j = 0; j < 4; j++)
            {
                int Tn_id = Tneigh[i][j];
                if (Tn_id >= 0)
                {
                    tot_ifact += abs(ifact_cache[i * 4 + j]);
                }
            }

            if (tot_ifact > epsilon)
            {
                for (int j = 0; j < 4; j++)
                {
                    int Tn_id = Tneigh[i][j];
                    if (Tn_id >= 0)
                    {
                        new_w += ifact_cache[i * 4 + j] / tot_ifact * w[Tn_id];
                    }
                }
                w[i] = new_w;
                tot_w += w[i];
            }

        }
        cout << "iter: " << iter <<" tot_w : "<<tot_w<< endl;
    }
#endif
    for (int i = 0; i < NT; i++)
    {
        if (w[i] > 0.0)
        {
            TetStatus[i] = PC_OUTSIDE;
        }
        else
        {
            TetStatus[i] = PC_INSIDE;
        }
    }

    Deallocate(&w);

    /*extract surface between in out Tetras*/
    GenerateTriangles();

    Deallocate(&TetStatus);
}


#endif
