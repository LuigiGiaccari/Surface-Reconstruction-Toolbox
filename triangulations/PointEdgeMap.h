//Nuova versione della PointEdgeMap che evita l'uso dei vector

#ifndef __PointEdgeMap_h__
#define __PointEdgeMap_h__


#define InsertNode(p1,p2)  \
    idFirst=First[p1];\
    First[p1]=counter;\
    Map[counter].Next=idFirst;\
    Map[counter].Point2Edge=idedge;\
    Map[counter].Point2Point=p2;

#define SearchEdge(p1,p2)  \
    id=First[p1];\
    while(id>=0)\
    {\
        if (Map[id].Point2Point==p2)\
        {return Map[id].Point2Edge;}\
        id=Map[id].Next;\
    }

struct PENode
{
    int Point2Point;
    int Point2Edge;
    int Next;
};

class PointEdgeMap
{
    private:

        int* First;//Get the first node to check
        PENode* Map;//Store for PEnodes
        int counter;//numero di nodi inseriti nella mappa

    public:

        PointEdgeMap();
        ~PointEdgeMap();

        int BuildPointEdgeMap(int Np, int MAXE);
        void AddEdge(int p1, int p2, int idedge);
        int GetEdge(int p1, int p2);
        void GetEdges(int p1);
        void Deallocate();
};


PointEdgeMap::PointEdgeMap()//constructor
{
    //do nothing just allocate memory
}

//Costruisce la pointedgeMap
int PointEdgeMap::BuildPointEdgeMap(int Np, int MAXE)
{
    //return:
    //     0: ok
    //     1: out of memory
    int i;
    counter = 0;

    First = new int[Np];
    if (First == NULL)
    {
        return 1;
    }
    for (i = 0; i < Np; i++)
    {
        First[i] = -1;
    }

    Map = new PENode[MAXE];
    if (Map == NULL)
    {
        delete [] First;    //out of memory (delete first and quit)
        return 1;
    }

    for (i = 0; i < MAXE; i++)
    {
        Map[i].Next = -1;
    }


    return 0;
}


PointEdgeMap::~PointEdgeMap()//destructor
{
    Deallocate();
}

void PointEdgeMap::Deallocate()//destructor
{
    if (First != NULL)
    {
        delete [] First;
    }
    First = NULL;
    if (Map != NULL)
    {
        delete [] Map;
    }
    Map = NULL;

}

//returns the edge id or -1
int PointEdgeMap::GetEdge(int p1, int p2)
{
    int id;

    if (p1 > p2) //sort p1 and p2
    {
        SearchEdge(p1, p2)
    }
    else
    {
        SearchEdge(p2, p1)
    }
    return -1;

}


// add new edge to the map
void PointEdgeMap::AddEdge(int p1, int p2, int idedge)
{
    int idFirst;

    if (p1 > p2) //sort p1 and p2
    {

        InsertNode(p1, p2)
    }
    else
    {
        InsertNode(p2, p1)
    }
    counter++;

}




#endif
