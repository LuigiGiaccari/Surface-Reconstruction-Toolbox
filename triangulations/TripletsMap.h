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



    Node* Map;//Store per i nodi
    int Np;//numero di nodi inseriti nella mappa
    PointEdgeMap* PEMap;
    Edge* e;

public:


    TripletsMap();//constructor
    ~TripletsMap();//destructor



//funzioni membro
    void TripletsMap::BuildTripletsMap(int N, PointEdgeMap* PEMap, Edge* ein);
    int TripletsMap::AnalyzeEdge(int p1,int p2);


};


TripletsMap::TripletsMap()//constructor
{
    //do nothing just allocate memory
}

//Costruisce la pointedgeMap
void TripletsMap::BuildTripletsMap(int N, PointEdgeMap* PEMapIn,Edge* ein)
{
    int i;
    Np=N;//numero di punti

    PEMap=PEMapIn;//salva il puntatore alla point 2 edge map

    e=ein;
//Dimensiono le mappe in base al numero di punti (inizializzo con -1 )

    Map=new Node[Np];

    //initializza to -1
    for(i=0; i<Np; i++)
    {
        Map[i].s=Map[i].t=-1;
    }

}


TripletsMap::~TripletsMap()//destructor
{
    delete [] Map;//deallocate
}


int TripletsMap::AnalyzeEdge(int p1,int p2);
{
    int p4;
    int idedge;

    if (Map[p1].s>0)//il punto era già allocato nella mappa
    {
        //edge p1-p4 was in the map
        p4=Map[p1].s;//trova il terzo punto chiamato P4 secondo le precedenti convenzioni
    }
    else if(Map[p2].s>0)
    {
        //edge p1-p4 was in the map
        p4=Map[p2].s;//trova il terzo punto chiamato P4 secondo le precedenti convenzioni
    }
    else//aggiungi i punti alla mappa
    {
        Map[p1].s=p2;
        Map[p2].s=p1;
        return -1;
    }//no triangle found}


    idedge=PEMap-> GetEdge( p2, p4);

    if (idedege>0 && e[idedge].t2<0)//se l'edge esiste ed è boundary abbiamo trovato un triangolo
    {
        return p4;
    }

}