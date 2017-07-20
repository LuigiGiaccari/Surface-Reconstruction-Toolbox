//this class implements a dynamic hash tecnique for fast enclosing triangle location
//
// the map is defined for a maximum of N elements. The hashing function can be changed according the current number of elements
//It is a good choice to change the hash function everytime the current number reaches a n^2 value
//this tecnique should garantee a more accurate location of target when querying the map.
#ifndef _TMAP_2D_H
#define _TMAP_2D_H
#include "quickdel_datas.h"



class TMAP_2D
{


    private:
        int* Map; //main hashinf array
        int max_edge_size;
        int max_number_elements;
        int current_number_elements;
        int current_edge_size;
        int current_edge_break;
        void Change_Hash_Func();//changes the hash function

    public:
        TMAP();
        ~TMAP();
        void BuildTMap(int number_of_points, int number_of_elements ); //build the map according to the maximum number of points and maximum number of map position
        void AddTriangle(int idx, int idy, int idt); //adds the triagnle idt into the coordinate idx idy
        void GetTriangle(int idx, int idy); //get traignle from the map

};

#endif