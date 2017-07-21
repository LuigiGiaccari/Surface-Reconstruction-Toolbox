#ifndef __TripletsMap_h__
#define __TripletsMap_h__
#endif

//includes
#ifndef __PointEdgeMap_h__
#include "PointEdgeMap.h"
#endif

struct TripletsMapNode
{
    int s;//second
    int t;//third
};
typedef TripletsMapNode Node;

class TripletsMap
{
    private:

        Node* Map;//node storage
        int Np;//number of inserted nodes in the map
        PointEdgeMap* PEMap;
        Edge* e;

    public:

        TripletsMap();//constructor
        ~TripletsMap();//destructor

        void TripletsMap::BuildTripletsMap(int N, PointEdgeMap* PEMap, Edge* ein);
        int TripletsMap::AnalyzeEdge(int p1, int p2);

};


TripletsMap::TripletsMap()//constructor
{
    //do nothing just allocate memory
}

//Costruisce la pointedgeMap
void TripletsMap::BuildTripletsMap(int N, PointEdgeMap* PEMapIn, Edge* ein)
{
    int i;
    Np = N;

    PEMap = PEMapIn;

    e = ein;

    Map = new Node[Np];

    //initialize to -1
    for (i = 0; i < Np; i++)
    {
        Map[i].s = Map[i].t = -1;
    }

}


TripletsMap::~TripletsMap()//destructor
{
    delete [] Map;//deallocate
}


int TripletsMap::AnalyzeEdge(int p1, int p2);
{
    int p4;
    int idedge;

    if (Map[p1].s > 0) //point was alrady in the map
    {
        //edge p1-p4 was in the map
        p4 = Map[p1].s;
    }
    else if (Map[p2].s > 0)
    {
        //edge p1-p4 was in the map
        p4 = Map[p2].s;
    }
    else//add points in the map
    {
        Map[p1].s = p2;
        Map[p2].s = p1;
        return -1;
    }//no triangle found}

    idedge = PEMap-> GetEdge( p2, p4);

    //the edge exists and it is a boundary one, we just found a candidate triangle
    if (idedege > 0 && e[idedge].t2 < 0)
    {
        return p4;
    }

}