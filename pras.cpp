#include<iostream>
#include<omp.h>

using namespace std;

int main()
{


	omp_set_num_threads(3);
	#pragma omp atomic

	#pragma omp barrier

	#pragma omp for
	{
		/* code */
	}
	int k=10;

	#pragma omp parallel shared(x,k) private(k) firstprivate(c)
	{
		code(1);
		k =10;
		k = omp_get_thread_num();
		omp_get_num_threads();
		#pragma omp master
		{

		}

		#pragma omp barrier
		code(2);
	}
	start = k*NUM_POINTS/NUM_THREADS;
	 , end =(k+1)*NUM_POINTS/NUM_THREADS;
	return 0;
};