//JUST LIKE THE SEDGE SORT LIBRARY BUT ALSO SORT INDEXES
// ho fatto le seguenti modifiche:
//
// -cambiato l'insertion sort
// -cambiato SWAP MACRO


#ifndef CUTOFF
#  define CUTOFF 15//insertion sort tolerance
#endif

//Swaps both integer and doubles values between indexes i1 and i2;
#define SWAP(i1,i2) tempd = a[i1];a[i1] = a[i2]; a[i2] = tempd;\
    tempi = index[i1]; index[i1] = index[i2]; index[i2] = tempi;



#define GT(x, y) ((x) > (y))
#define LT(x, y) ((x) < (y))
#define GE(x, y) ((x) >= (y))
#define LE(x, y) ((x) <= (y))
#define EQ(x, y) ((x) == (y))
#define NE(x, y) ((x) != (y))

/*
 * This is the SWAP macro to swap between two keys.
 * Replace these macros by YOUR swap macro.
 * e.g. if you are sorting an array of pointers to strings
 * You can define it as:
 *
 *	#define SWAP(x, y) temp = (x); (x) = (y); (y) = temp
 *
 * Bug: 'insort()' doesn't use the SWAP macro.
 */

/*-------------------------------------------------------------------
 *		  This file shouldn't be touched.
 *	     For customizable parameters, see 'sort.h'
 *-----------------------------------------------------------------*/

/*
 |  void  insort (array, len)
 |  KEY_T  array[];
 |  int    len;
 |
 |  Abstract:	Sort array[0..len-1] into increasing order.
 |
 |  Method:	Optimized insertion-sort (ala Jon Bentley)
 */

void  insort_index (double* a, int* index, int len)

{
    int	i, j;
    double	tempd;//temp double
    int tempi;//temp integer


    for (i = 1; i < len; i++)
    {
        /* invariant:  array[0..i-1] is sorted */
        j = i;
        /* customization bug: SWAP is not used here */
        tempd = a[j];
        tempi = index[j];
        while (j > 0 && GT(a[j - 1], tempd))
        {
            a[j] = a[j - 1];
            index[j] = index[j - 1];
            j--;
        }
        a[j] = tempd;
        index[j] = tempi;
    }
}






/*
 |  void  partial_quickersort (array, lower, upper)
 |  KEY_T  array[];
 |  int    lower, upper;
 |
 |  Abstract:
 |	Sort array[lower..upper] into a partial order
 |     	leaving segments which are CUTOFF elements long
 |     	unsorted internally.
 |
 |  Efficiency:
 |	Could be made faster for _worst_ cases by selecting
 |	a pivot using median-of-3. I don't do it because
 |	in practical cases my pivot selection is arbitrary and
 |	thus pretty random, your mileage may vary.
 |
 |  Method:
 |	Partial Quicker-sort using a sentinel (ala Robert Sedgewick)
 |
 |  BIG NOTE:
 |	Precondition: array[upper+1] holds the maximum possible key.
 |	with a cutoff value of CUTOFF.
 */

void  partial_quickersort_index (double* a, int* index, int lower, int upper)

{
    int	i, j;
    double	tempd, pivot;
    int tempi;


    if (upper - lower > CUTOFF)
    {

        SWAP(lower, (upper + lower) / 2);
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
            SWAP(i, j);
        }
        SWAP(lower, j);
        partial_quickersort_index (a, index, lower, j - 1);
        partial_quickersort_index (a, index, i, upper);
    }
}


/*
 |  void  sedgesort (array, len)
 |  KEY_T  array[];
 |  int    len;
 |
 |  Abstract:
 |	Sort array[0..len-1] into increasing order.
 |
 |  Method:
 |	Use partial_quickersort() with a sentinel (ala Sedgewick)
 |	to reach a partial order, leave the unsorted segments of
 |	length == CUTOFF to a simpler low-overhead, insertion sort.
 |
 |	This method is the fastest in this collection according
 |	to my experiments.
 |
 |  BIG NOTE:
 |	precondition: array[len] must hold a sentinel (largest
 |	possible value) in order for this to work correctly.
 |	An easy way to do this is to declare an array that has
 | 	len+1 elements [0..len], and assign MAXINT or some such
 |	to the last location before starting the sort (see sorttest.c)
 */
void  sedgesort_index (double* a, int* index, int len)

{
    /*
     * ------------------------- NOTE --------------------------
     * ignoring BIG NOTE above may lead to an infinite loop here
     * ---------------------------------------------------------
     */
    partial_quickersort_index (a, index, 0, len - 1);
    insort_index (a, index, len); //insertion sort on all the array
}





