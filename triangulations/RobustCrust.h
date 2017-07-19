//Library that uses a PowerCrust based algorithm for surface recon

#ifndef _ROBUSTCRUST_H_
#define _ROBUSTCRUST_H_

#ifndef M_PI
#define M_PI 3.14159265358979323846 //pigreco
#endif


//INCLUDES
#include <iostream>
#include "structs/Triangle.h"
#include "triangulations/Coord3D.h"
#include "sortlib/HCPO3D.h"
#include "triangulations/TrianglesFunctions.h"
#include "util/ArraysLib.h"
#include "util/filemanager.h"
#include "tetgen1.4.3/predicates.cxx"
#include "tetgen1.4.3/tetgen.cxx"
#include <queue>
#include <stack>


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
        RCRUST();//constructor
        ~RCRUST();//desstructor

    private:
        tetgenio in, out;//for delaunnay triangulation
        Coord3D* p;
        int N;//numbe rof cloud points
        int Nshield;//number of shiedl points
        Tetrahedra* T;
        Neigh* Tneigh;
        int NT;
        int ct;//triangles iterator
        Model3D Model;//infos about the model
        priority_queue<PQueueNode, vector<PQueueNode>, CompareNode > queue;
        bool* deleted;
        bool* checked;


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
        void LookForTriangle(vector<int>* Stack, int T1, int T2);
        int FindUnClustered();
        void Clustering();
        void ManifoldExtraction();
        void FillWedges();
        void DeleteSlivers();
        double TetraVol(int idT);
        void CopyPoints(vector<double>* inputp, int innp);
        void AddShield(Coord3D* p, int N);
};

RCRUST::RCRUST()//constructor
{
    ct = 0;
    checked = NULL;
    deleted = NULL;
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



    Deallocate(&checked);
    Deallocate(&deleted);

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
    behav.conformdel = true;
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

    std::cout << "Nt= " << NT << std::endl;
    std::cout << "Maximum value= " << Max(&T[0].p1, NT * 4) << std::endl;
    std::cout << "Maximum value= " << Max(&Tneigh[0].T1, NT * 4) << std::endl;



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
    //compute intersection factor between 2 tetredrons

    double distcc;
    double Ifact;

    SquaredDistance(cc1->, cc2->, distcc);

    Ifact = (*r1** r1 + *r2** r2 - distcc) / (2 * *r1** r2); //Carnot teorem
#ifdef _DEBUG
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
            checked[i] = deleted[i] = true;
        }

    }
    //add to queue all neighbours of shield tetraedrons

    for (i = 0; i < NT; i++)
    {
        if (deleted[i])
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

    /* //Old version


    Model.r*=10;
    for (i=0;i<NT;i++)
        {
           if(Tneigh[i].T1<0)AnalyzeIntersection(Tneigh[i].T1,i);
           if(Tneigh[i].T2<0)AnalyzeIntersection(Tneigh[i].T2,i);
           if(Tneigh[i].T3<0)AnalyzeIntersection(Tneigh[i].T3,i);
           if(Tneigh[i].T4<0)AnalyzeIntersection(Tneigh[i].T4,i);
        }

    //New version
    double Vol;
     for (i=0;i<NT;i++)
        {
           if(Tneigh[i].T1<0 || Tneigh[i].T2<0 || Tneigh[i].T3<0|| Tneigh[i].T4<0)//Tetredron is boundary
               {Vol=TetraVol(i);
                if (Vol<1e-14 &&Vol >-1e-14)continue;//jump slivers
               TetraCC(i,&cc,&r);
           if(r>Model.r) {Node.Ifact=1;//Tetraedrons with circumcenter bigger  than the model radius can not belogn to the surface
                    Node.T1=-1;
                    Node.T2=i;
                    Node.plevel=ComputePlevel(Node.Ifact);
                    queue.push(Node);}//add to queue


               }
        }
        */

}


void RCRUST::AnalyzeIntersection(int T1, int T2)
{
    // Analyze the intersection between two tetrahedra. Data are stored in the priority queue
    // T1,T2 is always >0; The marking flow goes from T1 to T2
    Coord3D cc1, cc2;
    double r1, r2;

    //stop for debug
#ifdef _DEBUG
    if (T1 == -7581 && T2 == 12)
    {
        T1 = T1;
    }
#endif

    PQueueNode Node;

    if (checked[T2])
    {
        return;    //no need to check go forward
    }

    TetraCC(T1, &cc1, &r1);
    TetraCC(T2, &cc2, &r2);

    //NOTA: sperimentare questa procedura!
    //if(r2>Model.r) Node.Ifact=1;//Tetraedrons with circumcenter bigger  than the model radius can not belogn to the surface

    Node.Ifact = Ifact(&cc1, &r1, &cc2, &r2);
    Node.T1 = T1;
    Node.T2 = T2;
    Node.plevel = ComputePlevel(Node.Ifact);

    queue.push(Node);//add to queue

}

float RCRUST::ComputePlevel(float Ifact)
{
    //definisce la strategia di priorità in base all'Ifact
    float plevel;

    //Il livello di priorità deve favorire le intersezioni molto forti Ifact>0
    // rispetto alle intersezioni molto deboli, Ifact<0
    // Ovviamente un intersezione molto debole è da preferire ad una intersezione poco forte

    //plevel varia fra 1 (altissima priorità) e 0 (scarsissima priorità)
    if (Ifact < 0) //separation
    {
        ;//intersezione molto debole rendila poco debole
        plevel = -Ifact; //rendilo positivo fra 0 e 1,=0 per Ifact=-1

    }
    else //intersection
    {
        plevel = Ifact; //rendilo positivo fra 0 e 1,=0 per Ifact=1
    }

    return plevel;
}


void RCRUST::FindFacet(int T1, int T2, Triangle* facet)
{
    //retunrs the facet shared by T1 and T2, T1 is a not deleted tetraedron
    //getting the position of T2
    int p1, p2, p3, p4;
    //    Coord3D tnorm;
    if (Tneigh[T1].T1 == T2)
    {
        p1 = T[T1].p2;
        p2 = T[T1].p3;
        p3 = T[T1].p4;
        p4 = T[T1].p1;
    }
    else if (Tneigh[T1].T2 == T2)
    {
        p1 = T[T1].p1;
        p2 = T[T1].p3;
        p3 = T[T1].p4;
        p4 = T[T1].p2;
    }
    else if (Tneigh[T1].T3 == T2)
    {
        p1 = T[T1].p1;
        p2 = T[T1].p2;
        p3 = T[T1].p4;
        p4 = T[T1].p3;
    }
    else if (Tneigh[T1].T4 == T2)
    {
        p1 = T[T1].p1;
        p2 = T[T1].p2;
        p3 = T[T1].p3;
        p4 = T[T1].p4;
    }
    else
    {
        Error("Facet error");
    }

    //Check the orientation of p4 and i case swap triagnle
    double o = orient3d(&p[p1].x, &p[p2].x, &p[p3].x, &p[p4].x);

    if (o > 0)
    {
        facet->p1 = p1;
        facet->p2 = p2;
        facet->p3 = p3;
    }
    else
    {
        facet->p1 = p2;
        facet->p2 = p1;
        facet->p3 = p3;
    }



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

        if (Node.T2 < 0 || checked[Node.T2])
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
    bool status;

    if (checked[T2])
    {
        return;    //no need to recheck
    }

    if (T1 >= 0)
    {
        status = deleted[T1];
    }
    else
    {
        status = true; //outside tetraedrons are deleted
    }



    if (Ifact > 0) //intersection marl as equal
    {
        deleted[T2] = status;
    }
    else
    {
        deleted[T2] = !status; //not interection mark as different
    }

    //T2 is now checked
    checked[T2] = true;
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

    vector<int> Stack;
    int T1;

    Alls(checked, NT, false); //Reset

    //Start from tetra 0
    Stack.push_back(0);

    while (!Stack.empty())
    {
        T1 = Stack.back();
        Stack.pop_back();//get and delete element


        LookForTriangle(&Stack, T1, Tneigh[T1].T1);
        LookForTriangle(&Stack, T1, Tneigh[T1].T2);
        LookForTriangle(&Stack, T1, Tneigh[T1].T3);
        LookForTriangle(&Stack, T1, Tneigh[T1].T4);
        checked[T1] = true; //T1 is now checked
    }

}

void RCRUST::LookForTriangle(vector<int>* Stack, int T1, int T2)
{
    //if THey have different deleted a traignle is generated
    //T1 Always>0

    Triangle facet;
    bool deleted2;

    if (T2 >= 0)
    {
        if (checked[T2])
        {
            return;   //already checked skip
        }
        deleted2 = deleted[T2];
        Stack->push_back(T2);
    }//add to stack}
    else
    {
        if (checked[T1])
        {
            return;
        }
        deleted2 = true;
    }//outside tetraedrons are deleted



    if (deleted[T1] != deleted2) //need to generate a triagnle
    {
        //Finds the factes using using the inside tetraedrons to set outwar normal
        if (deleted[T1])
        {
            FindFacet(T2, T1, &facet);
        }
        else
        {
            FindFacet(T1, T2, &facet);
        }
        if (facet.p1 < N && facet.p3 < N && facet.p3 < N )
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
    AllocateAndInit(&checked, NT, false); //true for checked tetraedroms
    AllocateAndInit(&deleted, NT, false); //true for deleted tetraedroms


    InitQueue();

    Marking();

    //DeleteSlivers();

    Clustering();

    //FillWedges();


    GenerateTriangles();

    cout << "Genrated: " << t.size() << " triangles From " << NT << " Tetraedrons" << endl;



    //remove from memory
    Deallocate(&checked);
    Deallocate(&deleted);

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

    Alls(checked, NT, false); //Reset



    idClust = 0;
    while (1)
    {
        T1 = FindUnClustered();
        if (T1 < 0)
        {
            break;
        }
        //Start from tetra 0

        Stack.push_back(T1);
        while (!Stack.empty())
        {
            T1 = Stack.back();
            Stack.pop_back();//get and delete element
            checked[T1] = true; //T1 is now checked
            ClustId[T1] = idClust;
            T2 = Tneigh[T1].T1;
            if (T2 >= 0 && !deleted[T2] && !checked[T2])
            {
                Stack.push_back(T2);
            }
            T2 = Tneigh[T1].T2;
            if (T2 >= 0 && !deleted[T2] && !checked[T2])
            {
                Stack.push_back(T2);
            }
            T2 = Tneigh[T1].T3;
            if (T2 >= 0 && !deleted[T2] && !checked[T2])
            {
                Stack.push_back(T2);
            }
            T2 = Tneigh[T1].T4;
            if (T2 >= 0 && !deleted[T2] && !checked[T2])
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
            deleted[i] = false;
        }
        else
        {
            deleted[i] = true;
        }
    }

}


int RCRUST::FindUnClustered()
{
    //finds the first unclustered tetraedron id
    //returns -1 if all points are clustered
    int i;

    for (i = 0; i < NT; i++)
    {
        if (!checked[i] && !deleted[i])
        {
            return i;
        }
    }
    return -1;

}

void RCRUST::ManifoldExtraction()
{
    //Extract the Manifold performing  a walking operation


}
void RCRUST::FillWedges()
{
    //Fills wedge tetraedrons which have 3 triangles on the surface
    int i;
    char c;//int8
    for (i = 0; i < NT; i++)
    {
        if (!deleted[i])
        {
            continue;    //already filled
        }
        c = 0;
        if (Tneigh[i].T1 < 0)
        {
            continue;
        }
        else if (!deleted[Tneigh[i].T1])
        {
            c++;
        }
        if (Tneigh[i].T2 < 0)
        {
            continue;
        }
        else if (!deleted[Tneigh[i].T2])
        {
            c++;
        }
        if (Tneigh[i].T3 < 0)
        {
            continue;
        }
        else if (!deleted[Tneigh[i].T3])
        {
            c++;
        }
        if (Tneigh[i].T4 < 0)
        {
            continue;
        }
        else if (!deleted[Tneigh[i].T4])
        {
            c++;
        }
        if (c == 3)
        {
            deleted[i] = false;    //fill wedge
        }
    }
}


void RCRUST::DeleteSlivers()
{
    //delete all tetraedrons with almost null value
    int i;
    double det;
    for (i = 0; i < NT; i++)
    {
        det = TetraVol(i);
        if (det < 1e-14 && det > -1e-14)
        {
            deleted[i] = true;
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

#endif
