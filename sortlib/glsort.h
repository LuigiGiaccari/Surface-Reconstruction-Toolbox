#include <iostream>
#include "sortlib\sortn.h"

#define INSERTSORT 25//insetion sort value
#define SORTN 7//


typedef double Item;

#define STAT
#ifdef STAT
int QS = 0;
int IS = 0;
int S2 = 0;
int S3 = 0;
int S4 = 0;
int S5 = 0;
int S6 = 0;
int S7 = 0;
int S8 = 0;
#endif


void sortN(double* a, int bucket)
{
    //Binary search to locate the correct function
    double temp;

    if (bucket > 4) //5-6-7
    {
        if (bucket == 7)
        {
            sort7(a);
#ifdef STAT
            S7++;
#endif
        }
        //5-6
        else if (bucket == 6)
        {
            sort6(a);
#ifdef STAT
            S6++;
#endif
        }
        else
        {
            sort5(a);
#ifdef STAT
            S5++;
#endif
        }

    }
    else//3-4-2
    {
        if (bucket == 4)
        {
            sort4(a);
#ifdef STAT
            S4++;
#endif
        }
        else if (bucket == 3)
        {
            sort3(a);
#ifdef STAT
            S3++;
#endif
        }
        else if (a[1] < a[0])
        {
            temp = a[0];
            a[0] = a[1];
            a[1] = temp;
        }
#ifdef STAT
        S2++;
#endif

    }

}


void insertionSort(double* x, int length)
{
#ifdef STAT
    IS++;
#endif


    double key;
    int i;

    for (int j = 1; j < length; j++)
    {
        key = x[j];
        i = j - 1;
        while (x[i] > key && i >= 0)
        {
            x[i + 1] = x[i];
            i--;
        }
        x[i + 1] = key;
    }
}









void quickersort(double*  a, int lower, int upper)
{


#ifdef STAT
    QS++;
#endif

    int	i, j, bucket = upper - lower + 1;
    double	temp, pivot;





    SWAP(a[lower], a[(upper + lower) / 2]);
    i = lower;
    j = upper + 1;
    pivot = a[lower];

    while (1)
    {
        /*
         * ------------------------- NOTE --------------------------
         * ignoring BIG NOTE above may lead to an infinite loop here
         * ---------------------------------------------------------
         */
        do
        {
            i++;
        }
        while (LT(a[i], pivot));
        do
        {
            j--;
        }
        while (GT(a[j], pivot));
        if (j < i)
        {
            break;
        }
        SWAP(a[i], a[j]);
    }
    SWAP(a[lower], a[j]);//pu teh median at the center of the list






    /*  recursion */
    bucket = j - lower + 1;
    if (bucket > INSERTSORT)
    {
        quickersort(a, lower, j);   //quick sort (c'era un più uno a j che ho levato)
    }
    else if (bucket > SORTN)
    {
        insertionSort(&a[lower], bucket);  //insert sort
    }
    else if (bucket > 1)
    {
        sortN(&a[lower], bucket);
    }



    bucket = upper - i + 1;
    if (bucket > INSERTSORT)
    {
        quickersort(a, i, upper);   //quick sort
    }
    else if (bucket > SORTN)
    {
        insertionSort(&a[i], bucket);  //insert sort
    }
    else if (bucket > 1)
    {
        sortN(&a[i], bucket);
    }


}


void GLSort(double*  a, int lower, int upper)
{

    int bucket = upper - lower + 1;

    if (bucket > INSERTSORT)
    {
        quickersort(a, lower, upper);   //quick sort (c'era un più uno a j che ho levato)
    }
    else if (bucket > SORTN)
    {
        insertionSort(&a[lower], bucket);  //insert sort
    }
    else if (bucket > 1)
    {
        sortN(a, bucket);
    }
}

#ifdef STAT

void PrintGLSortStat()
{
    cout << "FUNCTIONS CALLS" << endl;
    cout << "QS= " << QS << endl;
    cout << "IS= " << IS << endl;
    cout << "S3= " << S3 << endl;
    cout << "S4= " << S4 << endl;
    cout << "S5= " << S5 << endl;
    cout << "S6= " << S6 << endl;
    cout << "S7= " << S7 << endl;
    cout << "S8= " << S8 << endl;


}
#endif

