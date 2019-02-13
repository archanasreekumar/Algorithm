//Instructions to run the c program:
//gcc filename.c -o filename
//./filename
//output file 'time of sorts.txt' will be generated in the same folder
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

//insertion sort function
void insertionSort(int arr[], int size)
{
   int i, key, j;
   for (i = 1; i < size; i++)
   {
       key = arr[i];
       j = i-1;
 
       while (j >= 0 && arr[j] > key)
       {
           arr[j+1] = arr[j];
           j = j-1;
       }
       arr[j+1] = key;
   }
}
 
//selection sort function


void swap(int *xp, int *yp)
{
    int temp = *xp;
    *xp = *yp;
    *yp = temp;
}


void selectionSort(int arr[], int size)
{
    int i, j, min;

    
    for (i = 0; i < size-1; i++){

        min = i;
        for (j = i+1; j < size; j++)//finding min element
          if (arr[j] < arr[min])
            min = j;

 
        swap(&arr[min], &arr[i]);
    }
}


//bubble sort function

void swap2(int *xp2, int *yp2)
{
    int temp2 = *xp2;
    *xp2 = *yp2;
    *yp2 = temp2;
}


void bubbleSort(int arr[], int size)
{
   int i, j;
   for (i = 0; i < size-1; i++)
       for (j = 0; j < size-i-1; j++)
           if (arr[j] > arr[j+1])
              swap2(&arr[j], &arr[j+1]);
}

//merge sort functions

void merge(int arr[], int l, int m, int r) //function to merge the sorted arrays
{ 
    int i, j, k; 
    int n1 = m - l + 1; 
    int n2 =  r - m; 
  
    int L[n1], R[n2]; //initializing 2 temporary arrays
  
    for (i = 0; i < n1; i++) //copying data to L from main array
        L[i] = arr[l + i]; 
    for (j = 0; j < n2; j++) //copying data to R from main array
        R[j] = arr[m + 1+ j]; 
  
    //merging  L and R to get final sorted array
    i = 0; // Initial index of first subarray 
    j = 0; // Initial index of second subarray 
    k = l; // Initial index of merged subarray 
    while (i < n1 && j < n2) //comparing elements of L&R and storing it in main array
    { 
        if (L[i] <= R[j]) 
        { 
            arr[k] = L[i]; 
            i++; 
        } 
        else
        { 
            arr[k] = R[j]; 
            j++; 
        } 
        k++; 
    } 
  
    while (i < n1) //copying remaining elements of L if any
    { 
        arr[k] = L[i]; 
        i++; 
        k++; 
    } 
  
    while (j < n2) //copying remaining elements of R if any
    { 
        arr[k] = R[j]; 
        j++; 
        k++; 
    } 
} 


void mergeSort(int arr[], int l, int r) //function for merging
{ 
    if (l < r) 
    { 
       
        int m = l+(r-l)/2; //to find the median;avoids overflow for large l &h
  
        
        mergeSort(arr, l, m); //sorting first half of array
        mergeSort(arr, m+1, r); //sorting second half of array
  
        merge(arr, l, m, r); //to merge the sorted sub-arrays
    } 
} 
   

//quick sort function

void Swap1(int *xp2, int *yp2)
{
    int temp2 = *xp2;
    *xp2 = *yp2;
    *yp2 = temp2;
}

  int Median3( int A[ ], int Left, int Right )
        {
            int Center = ( Left + Right ) / 2;

            if( A[ Left ] > A[ Center ] )
                Swap1( &A[ Left ], &A[ Center ] );
            if( A[ Left ] > A[ Right ] )
                Swap1( &A[ Left ], &A[ Right ] );
            if( A[ Center ] > A[ Right ] )
                Swap1( &A[ Center ], &A[ Right ] );

            // A[ Left ] <= A[ Center ] <= A[ Right ] 

            Swap1( &A[ Center ], &A[ Right - 1 ] );  // Hide pivot 
            return A[ Right - 1 ];                //Return pivot 
        }
 void quicksort( int A[ ], int Left, int Right )
        {
            int i, j;
            int Pivot;

      if( Left + 30 <= Right )
            {
         Pivot = Median3( A, Left, Right );
         i = Left; j = Right - 1;
          for( ; ; )
                {
             while( A[ ++i ] < Pivot ){ }
              while( A[ --j ] > Pivot ){ }
              if( i < j )
                  Swap1( &A[ i ], &A[ j ] );
                    else
                  break;
                }
          Swap1( &A[ i ], &A[ Right - 1 ] );  // Restore pivot 

         quicksort( A, Left, i - 1 );
         quicksort( A, i + 1, Right );
            }
            else  
          insertionSort( A + Left, Right - Left + 1 );// Do an insertion sort on the subarray 
        }


//heap sort function

void Swap(int *xp2, int *yp2)
{
    int temp2 = *xp2;
    *xp2 = *yp2;
    *yp2 = temp2;
}

  #define LeftChild( i )  ( 2 * ( i ) + 1 )

 void PercDown( int A[ ], int i, int N )
 {
    int Child;
    int Tmp;

    for( Tmp = A[ i ]; LeftChild( i ) < N; i = Child )
      {
       Child = LeftChild( i );
        if( Child != N - 1 && A[ Child + 1 ] > A[ Child ] )
             Child++;
          if( Tmp < A[ Child ] )
              A[ i ] = A[ Child ];
                else
              break;
      }
    A[ i ] =Tmp;
 }

 void heapsort( int A[ ], int N )
 {
   int i;
   for( i = N / 2; i >= 0; i-- )  // BuildHeap 
     PercDown( A, i, N );
   for( i = N - 1; i > 0; i-- )
      {
       Swap( &A[ 0 ], &A[ i ] );  // DeleteMax 
       PercDown( A, 0, i );
       }
    }


//main function

int main()
{ int r=0;
double time_i,time_i1,time_i2,time_s,time_s1,time_s2;
double time_b,time_b1,time_b2,time_q,time_q1,time_q2;
double time_h,time_h1,time_h2,time_m,time_m1,time_m2;
   do {
   printf("Enter lenth of array: ");
   int size;
   scanf("%d",&size);
   int arr[size];
   {for (int p=0;p < size;++p)
	   {arr[p] = rand() % 10000000 + 1;}

   }
  
//time complexity of insertion sort
clock_t t;
t = clock(); 
insertionSort(arr, size);//calling insertion sort
t = clock() - t;
double time_taken = ((double)t)/CLOCKS_PER_SEC; // in seconds
 if(r==0)
  {time_i = time_taken;}
  else if (r==1)
  {time_i1 = time_taken;}
 else 
  {time_i2 = time_taken;}


//time complexity of selection sort
clock_t t1;
t1 = clock();
selectionSort(arr, size);//calling selection sort
t1 = clock() - t1;
double time_taken1 = ((double)t1)/CLOCKS_PER_SEC; // in seconds

 if(r==0)
  {time_s = time_taken1;}
 else if (r==1)
  {time_s1 = time_taken1;}
 else 
  {time_s2 = time_taken1;}


//time complexity of quick sort
 clock_t t4;
 t4 = clock();
 quicksort(arr,0,size-1);//calling quick sort
 t4 = clock() - t4;
 double time_taken4 = ((double)t4)/CLOCKS_PER_SEC; // in seconds

 if(r==0)
   { time_q = time_taken4;}
 else if (r==1)
   { time_q1 = time_taken4;}
 else 
    {time_q2 = time_taken4;}

//time complexity of heap sort

 clock_t t5;
 t5 = clock();
 heapsort(arr,size);//calling heap sort
 t5 = clock() - t5;
 double time_taken5 = ((double)t5)/CLOCKS_PER_SEC; // in seconds

 if(r==0)
  {time_h = time_taken5;}
 else if (r==1)
  { time_h1 = time_taken5;}
 else 
  {time_h2 = time_taken5;}
 
//time complexity of merge sort
clock_t t6;
    t6 = clock();
    mergeSort(arr,0,size - 1);//calling merge sort
   t6 = clock() - t6;
   double time_taken6 = ((double)t6)/CLOCKS_PER_SEC; // in seconds

 if(r==0)
  {time_m = time_taken6;}
 else if (r==1)
  { time_m1 = time_taken6;}
 else 
  {time_m2 = time_taken6;}

//time complexity of bubble sort
 clock_t t3;
 t3 = clock();
 bubbleSort(arr, size);//calling bubble sort
 t3 = clock() - t3;
 double time_taken3 = ((double)t3)/CLOCKS_PER_SEC; // in seconds
 if(r==0)
   {time_b = time_taken3;}
 else if (r==1)
   {time_b1 = time_taken3;}
 else 
    {time_b2 = time_taken3;

// writing the output to file     
  FILE *fp;
  fp = fopen("time of sorts.txt","w");
  if (fp == NULL)
    {
        fprintf(stderr, "\nError opening file\n");
        exit (1);
    }
fprintf(fp,"\t \t \t \t n=100 \t \t \t \t \t   n=1000 \t \t \t \t n=10000 \n");

fprintf(fp,"insrtn \t \t \t \t %f \t \t \t \t %f \t \t \t \t %f \n",time_i, time_i1,time_i2);
fprintf(fp,"bubble \t \t \t \t %f \t \t \t \t %f \t \t \t \t %f \n",time_b, time_b1,time_b2);
fprintf(fp,"slctn \t \t \t \t %f \t \t \t \t %f \t \t \t \t %f \n",time_s, time_s1,time_s2);
fprintf(fp,"merge \t \t \t \t %f \t \t \t \t %f \t \t \t \t %f \n",time_m, time_m1,time_m2);
fprintf(fp,"quick \t \t \t \t %f \t \t \t \t %f \t \t \t \t %f \n",time_q, time_q1,time_q2);
fprintf(fp,"heap \t \t \t \t %f \t \t \t \t %f \t \t \t \t %f \n",time_h, time_h1,time_h2);
fclose(fp);

}

r++;
}while (r<= 2);


return 0;

}




  
