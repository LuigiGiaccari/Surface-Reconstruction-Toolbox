#include "Import_STL.h"
/////////////////////////////////////////////
//     costruttore   //
/////////////////////////////////////////////

STL_FILE_MANAGER::STL_FILE_MANAGER()
{
    NT = NULL;
    NP = NULL;
    v = NULL;
    t = NULL;
    p = NULL;
    Header = NULL;
    Type = 'n';
}
/////////////////////////////////////////////
//     distruttore   //
/////////////////////////////////////////////
STL_FILE_MANAGER::~STL_FILE_MANAGER()
{
    FreeMemory();
}

void STL_FILE_MANAGER::FreeMemory()
{
    //deallocate memory
    if (v != NULL)
    {
        delete [] v;
        v = NULL;
    }
    if (p != NULL)
    {
        delete [] p;
        p = NULL;
    }
    if (t != NULL)
    {
        delete [] t;
        t = NULL;
    }
    if (tnorm != NULL)
    {
        delete [] tnorm;
        tnorm = NULL;
    }
    if (Header != NULL)
    {
        delete [] Header;
        Header = NULL;
    }
    NP = 0;
    NT = 0;
    Type = 'n';
}

//The core function to import stl files
int STL_FILE_MANAGER::STL_Import(char* file_name)

{
    int flag;

    //Check whether the file is binary or ascii
    flag = CheckFileType(file_name);

    Type = 'a';
    //import file
    if (Type == 'b')
    {
        flag = ImportBinary(file_name);
    }
    else if (Type == 'a')
    {
        flag = ImportAscii(file_name);
    }


    return 0;
}

//Read stl data from a binary file
///   return 0; ok
///   return 1; //warning('stlread:nodata','No data in STL file.');
///   return 2; option different from zero this will cause a crash
// return 3 //can not open file
int STL_FILE_MANAGER::ImportBinary(char* file_name)
{
    size_t nread;//number of read bytes, to check fread runs correctly
    //float data[12];//temporary store for data
    unsigned __int16 option;//for the stl option
    FILE* pFile;//pointer to the stl file
    int i;//internal counter


    Header = new char[80]; //allocating memory for header
    if (Header == NULL)
    {
        FreeMemory();
        return 1;
    }//out of memory

    //Opening the stl binary file
    pFile = fopen(file_name, "rb");
    if (pFile == NULL)
    {
        return 3;   //can not open file
    }
    nread = fread(Header, sizeof(char), 80, pFile); //read the header

    //Bytes 81-84 are an unsigned 32-bit integer specifying the number of faces
    // that follow.
    nread = fread(&NT, sizeof(unsigned int), 1, pFile); //number of faces

    if (NT <= 0)
    {
        return 1;   //warning('stlread:nodata','No data in STL file.')
    }


    //allocate memory

    v = new Coord3D[NT * 3];
    tnorm = new Coord3D[NT];
    if ( v == NULL || tnorm == NULL)
    {
        FreeMemory();
        return 1;
    }//out of memory

    //loop trough all facets and import data

    for (i = 0; i < NT; i++)
    {

        nread = fread(&tnorm[i].x, sizeof(float), 3, pFile); //reading normals

        nread = fread(&v[i * 3].x, sizeof(float), 3, pFile); //reading first vertex
        nread = fread(&v[i * 3 + 1].x, sizeof(float), 3, pFile); //reading second vertex
        nread = fread(&v[i * 3 + 2].x, sizeof(float), 3, pFile); //reading third vertex


        //Skip the option WARNING IF OPTION !=0 THIS WILL CAUSE A CRASH
        nread = fread(&option, sizeof(unsigned __int16), 1, pFile);
        //if (option>=0) {return 2;} //can not open files with options}

    }

    return 0;
}



// Not unique verte to unique points /triangle data structure
int STL_FILE_MANAGER::v2pt()
{


    // Il sort3 può essere ottimizzato


    //macro to swap values
#define swap(x, y,temp) \
    temp = x; x = y; y = temp;

    //macro to sort 3 coordnates
#define Sort3(a,b,c,temp)\
    if (b < a) {swap(b, a,temp);}\
    if (c < b) {swap(c, b,temp);}\
    if (b < a) {swap(b, a,temp);}



    int i;//counters
    int c;//unique points counter
    int m, n; //conscutive vertex after sorting


    if (NT <= 0)
    {
        return 101;   //no vertes data
    }

    int* reindex = new int[NT * 3]; //map vertex index-->point index   pointid=reindex[i]
    Vertex3D* vcopy = new Vertex3D[NT * 3]; //copy for sorting
    bool* equal = new bool[NT * 3]; //to check if points are equal
    float temp;//temporary float for swap operations


    if ( vcopy == NULL || reindex == NULL || equal == NULL)
    {
        if (reindex != NULL)
        {
            delete [] reindex;
        }
        if (equal != NULL)
        {
            delete [] equal;
        }
        if (vcopy != NULL)
        {
            delete [] vcopy;
        }
        FreeMemory();
        return 1;
    }//out of memory


    //copying teh v into vcopy
    for (i = 0; i < NT * 3; i++)
    {
        vcopy[i].x = v[i].x;
        vcopy[i].y = v[i].y;
        vcopy[i].z = v[i].z;
        vcopy[i].id = i;
    }



    //Sort the vertex by rows
    for (i = 0; i < NT * 3; i++)
    {
        Sort3(vcopy[i].x, vcopy[i].y, vcopy[i].z, temp)
    }


    //soting the vertexes
    sort(&vcopy[0], &vcopy[NT * 3 - 1]);

    //print the first sorted values
    //for (i=0;i<10;i++)
    //  {c=index[i];
    //cout<<vcopy[c].x<<" "<<vcopy[c].y<<" "<<vcopy[c].z<<endl;}

    //count the number of unique points and build the reindex array

    //initilaize the first data and uniquepoints counter
    c = 0;
    reindex[vcopy[0].id] = 0;
    equal[0] = false;
    for (i = 1; i < NT * 3; i++)
    {
        m = i;
        n = i - 1;
        //compare two consecutive points
        equal[i] = ( (vcopy[m].x == vcopy[n].x) && (vcopy[m].y == vcopy[n].y) && (vcopy[m].z == vcopy[n].z));


        if (!equal[i])//vertex are not unique we can increase counter
        {
            c++;
        }

        // not unique vertex i equals unique points c;
        m = vcopy[i].id; //not unique vertex id
        reindex[m] = c; //that vertex equals the unique point c;
    }
    NP = c + 1; //setting the unique points number

    //build the unique points array

    p = new Coord3D[NP]; //we now know the number of unique points
    t = new Triangle[NT];

    if ( p == NULL || t == NULL)
    {
        delete [] reindex;
        delete [] equal;
        delete [] vcopy;
        FreeMemory();
        return 1;
    }//out of memory


    //saving the unique points value

    //first point
    for (i = 0; i < NT * 3; i++)
    {
        if (!equal[i])//to copy just once the unique points
        {
            m = vcopy[i].id; //id of unique vertex taken from sorted index
            c = reindex[m]; //id unique point
            p[c] = v[m]; //unique point m equals vertex i according to the map reindex
        }
    }

    //Writing the triangle data pointing to unique points
    c = 0; //now used as general counter
    for (i = 0; i < NT; i++)
    {
        t[i].p1 = reindex[c];
        t[i].p2 = reindex[c + 1];
        t[i].p3 = reindex[c + 2];
        c = c + 3;
    }


    //delallocate memory
    delete [] reindex;
    delete [] equal;
    delete [] vcopy;

    return 0;
}







//Read stl data from an ascii file
///   return 0; ok
///   return 1; //warning('stlread:nodata','No data in STL file.');
///   return 2; option different from zero this will cause a crash
// return 3 //can not open file
int STL_FILE_MANAGER::ImportAscii(char* file_name)
{
    FILE* pFile;//pointer to the stl file

    const int bufsize = 100;
    char line[bufsize];//buffer for reading

    int nlines;//number of lines
    float x, y, z; //temp values
    int cv, cn; //counters for vertex and normals
    int flag, i, j;
    char FirstChar;
    //Counting the number of points
    flag = CountLinesNumber(file_name, &nlines);

    //retrieving the number of triangles from the number of lines

    NT = (nlines - 3) / 7; //removing header and footer each facets has 7 lines (-3 for the header footer and end of file)

    pFile = fopen(file_name, "r"); //open the file
    if (!pFile)
    {
        return 2;   //file can not be opened
    }

    //allocating memory
    Header = new char[80];
    v = new Coord3D[NT * 3];
    tnorm = new Coord3D[NT];

    if (Header == NULL || v == NULL || tnorm == NULL)
    {
        FreeMemory();
        return 1;
    }//out of memory

    fgets(Header, 80, pFile); //get the file header (solid someting)

    //importing points the number of points
    cv = 0;
    cn = 0;
    do
    {
        fgets(line, bufsize, pFile); //get the line characters

        //getting the first non nll char
        for (j = 0; j < bufsize; j++)
        {
            if (line[j] != ' ')
            {
                FirstChar = line[j];
                break;
            }
        }
        //il primo carattere non nullo è j;

        if (FirstChar == 'v') //vertex data
        {
            i = sscanf(&line[j + 6], "%f%f%f", &x, &y, &z); //j+6 to jump "vertex"

            if (i != 3) //check if the data is readenn correclty

            {
                fclose(pFile);//close the file
                FreeMemory();
                return 3;
            }//wrong file format



            v[cv].x = x;
            v[cv].y = y;
            v[cv].z = z;
            cv++;
        }
        else if (FirstChar == 'f') //facet data
        {
            i = sscanf(&line[j + 12], "%f %f %f", &x, &y, &z); //j+12 to jump   facet normal
            if (i != 3) //check if the data is readenn correclty

            {
                fclose(pFile);//close the file
                FreeMemory();
                return 3;
            }//wrong file format

            tnorm[cn].x = x;
            tnorm[cn].y = y;
            tnorm[cn].z = z;
            cn++;
        }

    }
    while (!feof(pFile));

    fclose(pFile);//close the file


    if (cn != NT || cv != NT * 3) //problem in reading the file
    {
        FreeMemory();
        return 3;
    }//wrong file format




    return 0;
}


//retunrs the number of lines into a text file
//retunrs 0 if the file cna not be opened
int STL_FILE_MANAGER::CountLinesNumber(char* file_name, int* size)
{
    ifstream inData;

    inData.open( file_name );
    //if ( !inData )
    //return 0;


    *size = 0;
    string line;
    while ( !inData.eof() )
    {
        getline(inData, line);
        *size = *size + 1;
    }
    return 0;
}



int STL_FILE_MANAGER::CheckFileType(char* file_name)
{

    FILE* pFile;
    const int N = 30; //numbers of letters to be read
    char LastLine[N];//line to be read
    char* CompString = "endsolid";
    int i, j;
    bool found;


    pFile = fopen(file_name, "rb");
    if (!pFile)
    {
        return 2;//file ca not be openend
    }



    fseek(pFile, -N * sizeof(char), SEEK_END); //set the read to start N char before the end of file
    fread(&LastLine[0], sizeof(char), N, pFile);


    found = false;
    j = 0;
    for (i = 0; i < N; i++)
    {
        if (LastLine[i] == CompString[j])
        {
            j++;
            if (j >= 8) //length of "endsolid"
            {
                found = true;
                break;
            }
        }
        else
        {
            j = 0;
        }
    }


    if (found)
        //string found file is ascci
    {
        Type = 'a';
    }
    else
    {
        //string not found file is binary
        Type = 'b';
    }


    fclose(pFile);
    return 0;
}