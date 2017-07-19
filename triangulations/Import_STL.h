//A library to manage STL files

//Returned values flags:
// 0: ok
// 1: out of memory
// 2: file can not be openened
// 3: wrong file format
// above 100: special errors flags
// 101: (v2pt()) no vertex data

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
using namespace std;









class STL_FILE_MANAGER
{

private:


//Structure the make the computation easier
    struct Triangle
    {
        int p1;
        int p2;
        int p3;
    };
    struct Coord3D
    {
        float x;
        float y;
        float z;

        bool operator<(const Coord3D& p2) const
        {
            if (x<p2.x)
            {
                return true;   //x1<x2
            }

            if (x!=p2.x)
            {
                return false;   //x1>x2
            }

            if (y<p2.y)
            {
                return true;   //y1<y2
            }

            if (y!=p2.y)
            {
                return false;   //y1>y2
            }

            if (z<p2.z)
            {
                return true;
            }
            return false;//equals or z2<z1
        }

    };


    struct Vertex3D
    {
        float x;
        float y;
        float z;
        int id;

        bool operator<(const Vertex3D& p2) const
        {
            if (x<p2.x)
            {
                return true;   //x1<x2
            }

            if (x!=p2.x)
            {
                return false;   //x1>x2
            }

            if (y<p2.y)
            {
                return true;   //y1<y2
            }

            if (y!=p2.y)
            {
                return false;   //y1>y2
            }

            if (z<p2.z)
            {
                return true;
            }
            return false;//equals or z2<z1
        }

    };


//FUNCTIONS PROTOTYPES
    int CheckFileType(char* file_name);//cheks wheter a file si binary or ascii
    int CountLinesNumber(char* file_name,int* linesnum);//counts the number of lnes into an stl file
    int ImportAscii(char* file_name);//imports ascii stl file
    int ImportBinary(char* nome_file);//import binary stl files


public:

    char* Header;
    char Type;// n=nofile;b=binary file;a=ascii file;
    int NP;//number of unique points
    int NT;//number of triangles

    Coord3D* v;//list of vertex like stl format
    Triangle *t;//list of triangles indexes to points
    Coord3D *p;//list of unique points
    Coord3D *tnorm;//normals of triangles

    STL_FILE_MANAGER();//constructor
    ~STL_FILE_MANAGER();//destructor
    int STL_Import(char* nome_file);//import an stl file (binary or ascii)
    int v2pt();//gets unique points and triangles form vertex data
    void FreeMemory();//Free the memory


};