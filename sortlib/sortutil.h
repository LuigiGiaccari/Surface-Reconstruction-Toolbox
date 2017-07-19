


void PrintVect(double*p,int N)
{

    int i;
    printf("\n\n\n");
    for(i=0; i<N; i++)
    {
        printf("%4.2f \n",p[i]);
    }
}


void GenerateRandom(double* p,int n)
{
    //function to generate random points in 0-1 range


    // loop to generate random points not normalized
    int i;//counter;
    int tempmax=0;//temporary integer maximum random point;
    //srand( (unsigned)time( NULL ) );

    for( i = 0;   i <n; i++ )
    {

        p[i]=rand();
        if (p[i]>tempmax)
        {
            tempmax=p[i];
        }

    }


    //loop to normalize
    for( i = 0;   i <n; i++ )
    {
        p[i]=p[i]/tempmax;

    }
}



bool IsSorted(Item* v,int N)
{
    int i;
    bool sorted;
    for(i=0; i<N-1; i++)
    {
        if (v[i]>v[i+1])
        {
            return sorted=false;
        }
    }

    return sorted=true;

}



void Copy(double* p,double* p0,int N)
{
    int i;
    for(i=0; i<N; i++)
    {
        p[i]=p0[i];
    }
}