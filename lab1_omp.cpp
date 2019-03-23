#include "lab1_sequential.h"
#include<bits/stdc++.h> 
#include <omp.h>

using namespace std;


int checkN = 0;
int* for_loop ;
int check_count = 0;
omp_lock_t writelock;

void initializing_centroids( float* &dup_centroid_ss,int*& datapoints,int k){ 
    
    dup_centroid_ss = (float * )malloc(k* 3 * sizeof(float)); 
    
    int kar;
    srand(0);

    for(int i =0;i<k;i++){
        kar = rand()%checkN;
        int i_kar = 3*kar;
        int i_3 = 3*i;
        dup_centroid_ss[i_3] = datapoints[i_kar];
        dup_centroid_ss[i_3+1] = datapoints[i_kar+1];
        dup_centroid_ss[i_3+2] = datapoints[i_kar+2];
    }
       
}

float dist_p1_p2(int* &p1,float *&p2){
 
    float dist =0;

    for(int i =0;i<3;i++){

        dist = dist + (p1[i]-p2[i])*(p1[i]-p2[i]); 
    }
   
    return sqrt(dist);
}


void updating_centroids(int k,float*& centroid_ss)
{
  #pragma omp parallel
{
   #pragma omp  for
    for(int i =0;i<k;i++){
    	    int i_3 = 3*i;
    	    int i_4 = 4*i;
        	centroid_ss[i_3+0] =    float(for_loop[i_4+0])/for_loop[i_4+3]; // has to be float recheck
        	centroid_ss[i_3+1] =    float(for_loop[i_4+1])/for_loop[i_4+3];
        	centroid_ss[i_3+2] =    float(for_loop[i_4+2])/for_loop[i_4+3];
       
    }
}

}

void clustering(int t,int*& datapoint_cluster,int N,int k,int*& datapoints,float*& centroid_ss){
    
    omp_set_num_threads(t);

    #pragma omp parallel shared(datapoint_cluster,N,k,datapoints,centroid_ss) 
    {     	

    	int ii;
        int ind ;
        float dis ;
        float min_dis ;
        int* p1 = new int[3];
        float* p2 = new float[3];
        int temp_check_count = 0;
       	int *shared_data =new int[4*k]();
       	int temp_thread_num = omp_get_thread_num();
       	

      #pragma omp for 
    	for(int i= 0;i<N;i++){ 
       		// cout<< "tid "<<temp_thread_num  << " vali "<<i<<" N is "<<N<<endl;
            int i_d = 3*i;

        	p1[0] = datapoints[i_d+0];

        	p1[1] = datapoints[i_d+1];

        	p1[2] = datapoints[i_d+2];

        	min_dis = INT_MAX;
        // dis = 0;
        // ind = -1;
        
        for(ii=0;ii<k;ii++)
        {
            int i1_d = 3*ii;
            p2[0] = centroid_ss[i1_d];
            p2[1] = centroid_ss[i1_d+1];
            p2[2] = centroid_ss[i1_d+2];

            dis = dist_p1_p2(p1,p2);
             
            if(dis<min_dis)
            {
                min_dis = dis;
                ind = ii;
            } 
        }

        if(datapoint_cluster[4*i+3] != ind){
            datapoint_cluster[4*i+3] = ind;
            temp_check_count++;
            // temp_check_count[id] ++;
        }
     
        int i_ind = 4*ind;
        shared_data[i_ind ] += p1[0];
        shared_data[i_ind +1] += p1[1];
        shared_data[i_ind +2] += p1[2];
        shared_data[i_ind +3] ++; 

    }
    omp_set_lock(&writelock);
    	check_count+=temp_check_count;
       	for(int i=0;i<4*k;i++){
       		for_loop[i] += shared_data[i];
        }        
    omp_unset_lock(&writelock);
    free(p1);
    free(p2);
    free(shared_data);
}
// exit(1);
}
void centroids_final(int t,int N,float*& centroid,vector<float> &cvector,int iteration,int K){

       int i1  = cvector.size();
       centroid = (float *)malloc((i1) * sizeof(float)); 
       // #pragma omp parallel
       // {
       // 	int id = omp_get_thread_num();
       // 	int len = (i1)/t;
       // 	int start = id*len;
       // for(int i =start;i<start+len;i++){ 
       for(int i =0;i<cvector.size();i++){
       
               centroid[i] = cvector[i];    
        }  
    // }        
}

void kmeans_omp(int num_threads,
				int N,
				int K,
				int* data_points,
				int** data_point_cluster,
				float** centroids,
				int* num_iterations
				)
{

    checkN = N;
	vector<float> cvector;

	omp_init_lock(&writelock);
	
	*data_point_cluster = (int*)malloc(N*4*sizeof(int));

    for_loop = (int*)malloc(4*K*sizeof(int));

#pragma omp parallel
    {
       #pragma omp for
        for(int i =0;i<N;i++){
       	int i_4 = 4*i;
       	int i_3 = 3*i;
        (*data_point_cluster)[i_4+0] = data_points[i_3+0];
        (*data_point_cluster)[i_4+1] = data_points[i_3+1];
        (*data_point_cluster)[i_4+2] = data_points[i_3+2];
        (*data_point_cluster)[i_4+3] = 0;
    }
}

    float* centroid_ss; 
                    
    initializing_centroids(centroid_ss,data_points,K);

    bool change = true;
    
    int iteration =0;
    
    while(change){   

     	iteration = iteration + 1;
         
          for(int i0=0;i0<K;i0++){
          	int i_3 = 3*i0;
          	int i_4 = 4*i0;
   			cvector.push_back( centroid_ss[i_3]);
   			cvector.push_back( centroid_ss[i_3+1]);
   			cvector.push_back( centroid_ss[i_3+2]);
            for_loop[ i_4 + 0] = 0;
            for_loop[ i_4 + 1] = 0;
            for_loop[ i_4 + 2] = 0;
            for_loop[ i_4 + 3] = 0;
  		  }
  		
     
        clustering(num_threads,*data_point_cluster,N,K,data_points,centroid_ss);
        // cout<< "iteration "<<iteration<<" " <<check_count<<endl; 
        updating_centroids( K,centroid_ss);
        // updating_centroids( N,K,centroid_ss,*data_point_cluster);
        if(check_count<1){
            change = false;
         }
        else{
           change = true;	
         }      
        // change = check_update(N,K,cvector,centroid_ss,iteration);
        check_count =0;       
    }

    num_iterations[0] = iteration;
    
    centroids_final(num_threads,N,*centroids,cvector,iteration,K);
    
}