#pragma once
#ifndef __SMALLPRIORITYQUEUE_h__
#define __SMALLPRIORITYQUEUE_h__
//class that implesnt an afficent priority queue for a few priorities
//Level of priority are considered from 0 to N. Highest priority is 0;

//NOTE: this implementation is not robust: it does not consider errors check
// including memory errors or array overflow.
using namespace std;

class SMALLPRIORITYQUEUE
{
    private:
        int* queue;
        int Nlevel;//numbers of level of priority
        int* First;
        int iq;//queue iterator
    public:
        SMALLPRIORITYQUEUE();
        ~SMALLPRIORITYQUEUE();
        bool IsEmpty();
        bool Pop(int* id);
        void Push(int id, int plevel);
        void PushHigh(int id);
        void Set(int* queuepointer, int Nlevelinput);
        void Reset();
        void Empty();
};

//COnstructor that points to the already allocated array queue and set the input level numbers;
SMALLPRIORITYQUEUE::SMALLPRIORITYQUEUE()
{
    First = NULL;
    iq = -1;
    Nlevel = 0;
    queue = NULL;
}

SMALLPRIORITYQUEUE::~SMALLPRIORITYQUEUE()
{
    Reset();
}

void SMALLPRIORITYQUEUE::Set(int* queuepointer, int Nlevelinput) //Set the loaction and priority levels for the queue
{
    int i;
    queue = queuepointer;
    Nlevel = Nlevelinput;
    First = new int[Nlevel];
    for (i = 0; i < Nlevel; i++)
    {
        First[i] = 0;    //initialize
    }
    iq = -1;
}

void SMALLPRIORITYQUEUE::Reset()
{
    Nlevel = 0;
    iq = -1;
    if (First != NULL)
    {
        delete [] First;
    }
    First = NULL;
}


//query for empty queue
bool SMALLPRIORITYQUEUE::IsEmpty()
{
    if (iq < 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}


//add an highest priority element
void SMALLPRIORITYQUEUE::PushHigh(int id)
{
    iq++;
    queue[iq] = id;
}



//add a generic plevel priority element
//Warning: it deos not work for an highest priority element
void SMALLPRIORITYQUEUE::Push(int id, int plevel)
{
    int i;
    if (plevel >= Nlevel)
    {
        return;    //priority too low do not insert into the queue
    }
    iq++;
    queue[iq] = queue[First[0]]; //pushing the highest priority elemnt before swappping values
    for (i = 0; i < plevel; i++) //swap values form lowest to highest priority
    {
        queue[First[i]] = queue[First[i + 1]];
        First[i]++;
    }
    queue[First[plevel]] = id; //finally insert the value into the queue

}
//pop the highest priority element
bool SMALLPRIORITYQUEUE::Pop(int* id)
{
    int i;
    if (iq < 0)
    {
        return false;
    }
    *id = queue[iq];

#ifdef _DEBUG

    if (*id < 0)
    {
        cout << "BUG" << endl;
    }
#endif
    iq--;//decrrease iterator after popping
    if (iq < 0)
    {
        return true;   //per quando la coda arriva all'ultimo elemnto e poi si ricarica
    }



    for (i = 0; i < Nlevel; i++)
    {
        if (First[i] > iq)
        {
            First[i] = iq;
        }
        else
        {
            break;
        }
    }
    return true;
}

//query for empty queue
void SMALLPRIORITYQUEUE::Empty()//rest the iterator queue to first position
//Priority levels are untouched
{
    iq = -1; //set iterator to empty position
}
#endif