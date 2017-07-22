//Function to manage:
//- memory allocation for arrays;
//- basic operations on arrays (max,min,mean...);
//- Errors and Warnings


#pragma once
#ifndef __ArraysLib_h__
#define __ArraysLib_h__

#include <iostream>
#include "util/ExceptionHandler.h"
using namespace std;



//Functions Prototypes

template <class arr>void Allocate (arr**  a, int N); //allocate memory of an array

template <class arr>void AllocateAndInit (arr** a, int N, arr value);

template <class arr>void Deallocate (arr** a);//Deallocate Memory

template <class arr>void Alls(arr* a, int N, arr value) ; //Set all the array to a given value

template <class arr>void Random(arr* a, int N, arr min, arr max) ; //Generates random values between min and max

template <class arr> arr Min(arr* a, int N, int step = 1); //gets  minum

template <class arr> arr Max(arr* a, int N, int step = 1); //gets maximum

template <class arr> void MinMax(arr* a, int N, arr* max, arr* min, int step = 1); //gets maximum and minum

template <class arr>void PrintArray(arr* a, int N, int col, int step = 1); //Prints all the values of an array using the rows col tabulation

template <class arr>arr Mean(arr* a, int N, int step = 1); //gets the meanvalue

template <class arr1, class arr2> void Copy(arr1* destination, arr2* source, int N); //copies 2 arrays




//- memory allocation for arrays;
//CODE STARTS HERE

template <class arr>
void Allocate(arr** a, int N)
{

    if (*a != NULL)
    {
        Error("Can not Allocate");;
    }
    *a = new arr[N];
    if (*a != NULL)
    {
        return;
    }
    else
    {
        Error("Out of Memory");
    }
}

template <class arr>
void AllocateAndInit (arr** a, int N, arr value)
{



    if (*a != NULL)
    {
        Error("Can Not Allocate");;
    }
    *a = new arr[N];
    if (*a != NULL)
    {
        Alls(*a, N, value);
        return;
    }
    else
    {
        Error("Out of Memory");
    }
}


template <class arr>
void Deallocate (arr** a)
{
    if (*a == NULL)
    {
        return;
    }
    delete [] *a;
    *a = NULL;
    return;
}


template <class arr>
void Alls(arr* a, int N, arr value)
{
    int i;
    for (i = 0; i < N; i++)
    {
        a[i] = value;
    }
    return;

}


template <class arr>
void Random(arr* a, int N, arr min, arr max) //Generates random values between min and max
{
    int i;
    for (i = 0; i < N; i++)
    {
        a[i] = min + (arr)rand() / RAND_MAX * max;
    }
}

template <class arr>
arr Max(arr* a, int N, int step) //Gets The maximum value of the array
{
    //the variable step is use to iterate at step distance
    arr max = a[0];
    int i;
    int c = step; //jumping the first element
    for (i = 1; i < N; i++)
    {
        if (a[c] > max)
        {
            max = a[c];
        }
        c = c + step;
    }
    return max;
}

template <class arr>
arr Min(arr* a, int N, int step) //Gets The minum value of the array
{
    //the variable step is use to iterate at step distance
    arr min = a[0];
    int i;
    int c = step;
    for (i = 1; i < N; i++)
    {
        if (a[c] < min)
        {
            min = a[c];
        }
        c = c + step;
    }
    return min;
}


template <class arr>
void MinMax(arr* a, int N, arr* max, arr* min, int step) //Gets The maximum and minium value of the array
{
    //the variable step is use to iterate at step distance
    *min = a[0];
    *max = a[0] ;
    int i;
    int c = step;
    for (i = 1; i < N; i++)
    {
        if (a[c] < *min)
        {
            *min = a[c];
        }
        if (a[c] > *max)
        {
            *max = a[c];
        }
        c = c + step;
    }
}


template <class arr>
void PrintArray(arr* a, int N, int col, int step) //Prints all the values of an array using the rows col tabulation
{
    //the variable step is use to iterate at step distance

    int i, r, c;
    i = 0;
    r = 0;
    c = 0;

    while (1)
    {
        for (c = 0; c < col; c++) //loop trough rows

        {
            cout << " " << a[i] << " ";
            i += step;
            if (i >= N)
            {
                cout << endl;
                return;
            }
        }
        cout << endl;
    }

    return;
}

template <class arr>
arr Mean(arr* a, int N, int step) //gets the meanvalue
{
    //the variable step is use to iterate at step distance

    int i;
    int c = 0;
    arr mean = 0;
    for (i = 0; i < N; i++)
    {
        mean += a[c];
        c = c + step;
    }
    return mean / N;
}
template <class arr1, class arr2>
void Copy(arr1* destination, arr2* source, int N)
{
    int i;
    for (i = 0; i < N; i++)
    {
        destination[i] = source[i];
    }



}





#endif
