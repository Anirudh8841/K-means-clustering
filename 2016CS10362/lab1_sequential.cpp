#include "lab1_sequential.h"
#include<bits/stdc++.h> 
using namespace std;

int checkN = 0;
int check_count = 0;
int* for_loop ;


void initializing_centroids( float* &dup_centroid_ss,int*& datapoints,int k){ 
    dup_centroid_ss = (float * )malloc(k* 3 * sizeof(float)); 
    srand(0);

    int kar;
    for(int i =0;i<k;i++){
        kar = rand()%checkN;
        dup_centroid_ss[3*i] = datapoints[3*kar];
        dup_centroid_ss[3*i+1] = datapoints[3*kar+1];
        dup_centroid_ss[3*i+2] = datapoints[3*kar+2];
    }
}


float dist_p1_p2(int* &p1,float *&p2){
 
    float dist =0;

    for(int i =0;i<3;i++){

        dist = dist + pow(p1[i]-p2[i],2); 
    }

    return sqrt(dist);
}

void print_cent(float *& centroid_ss,int k){

    for(int i =0;i<3*k;i++){
        
        cout<<"print_cen "<< i <<" "<< centroid_ss[i]<<endl;
    }
}

void updating_centroids(int N,int k,float*& centroid_ss,int*& datapoint_cluster)
{

    for(int i =0;i<4*k;i++){
       for_loop[ i ] = 0;
    }

    for(int i =0;i<N;i++){

        for_loop[4* (datapoint_cluster[4*i+3]) + 0] += datapoint_cluster[4*i+0];
        for_loop[4* (datapoint_cluster[4*i+3]) + 1] += datapoint_cluster[4*i+1];
        for_loop[4* (datapoint_cluster[4*i+3]) + 2] += datapoint_cluster[4*i+2];
        for_loop[4* (datapoint_cluster[4*i+3]) + 3] ++;
    }

// averaging 
    for(int i =0;i<k;i++){
    	
        centroid_ss[3*i+0] = float(for_loop[4*i+0])/for_loop[4*i+3]; // has to be float recheck
        centroid_ss[3*i+1] = float(for_loop[4*i+1])/for_loop[4*i+3];
        centroid_ss[3*i+2] = float(for_loop[4*i+2])/for_loop[4*i+3];
    }

}

void clustering(int*& datapoint_cluster,int N,int k,int*& datapoints,float*& centroid_ss){


    int* p1 ;
    float* p2 ;
    p1 = (int*)malloc(3*sizeof(int));
    p2 = (float*)malloc(3*sizeof(float));

    for(int i=0;i<N;i++){

        p1[0] = datapoints[3*i+0];
        p1[1] = datapoints[3*i+1];
        p1[2] = datapoints[3*i+2];

        float min_dis = INT_MAX;
        float dis = 0;
        int ind = -1;

        for(int ii=0;ii<k;ii++){
            p2[0] = centroid_ss[3*ii+0];
            p2[1] = centroid_ss[3*ii+1];
            p2[2] = centroid_ss[3*ii+2];
            dis = dist_p1_p2(p1,p2);

            if(dis<min_dis){
                min_dis = dis;
                ind = ii;
            } 
        }

        if(datapoint_cluster[4*i+3] != ind){
            datapoint_cluster[4*i+3] = ind;
            check_count++;
        }
    }
    free(p1);
    free(p2);
}

bool check_update(int N,int k,vector<float> &cvector,float*& centroid_ss,int num_iterations){

    if(check_count<1){
        return false;
    }
    return true;
/*
    bool updated = false;
    for(int i =0;i<k;i++){
        if(cvector[3* (num_iterations-1) + 3*i] != centroid_ss[3*i]){
           updated = true;
           break;
        }
         if(cvector[3* (num_iterations-1) + 3*i+1] != centroid_ss[3*i+1]){
           updated = true;
           break;
        }
         if(cvector[3* (num_iterations-1) + 3*i +2] != centroid_ss[3*i+2]){
           updated = true;
           break;
        }
    return updated;

    }*/

}

void centroids_final(float*& centroid,vector<float> &cvector,int iteration,int K){

    centroid = (float *)malloc((iteration+1)*K* 3 * sizeof(float)); 

    for(int i =0;i<cvector.size();i++){
        centroid[i] = cvector[i];    
    }      
}

void kmeans_sequential(int N, int K, int* data_points, int** data_point_cluster, float** centroids, int* num_iterations )
{
    checkN = N;
	vector<float> cvector;    
    for_loop = (int*)malloc(4*K*sizeof(int));
    *data_point_cluster = (int*)malloc(N*4*sizeof(int));

    for(int i =0;i<N;i++){
        (*data_point_cluster)[4*i+0] = data_points[3*i+0];
        (*data_point_cluster)[4*i+1] = data_points[3*i+1];
        (*data_point_cluster)[4*i+2] = data_points[3*i+2];
        (*data_point_cluster)[4*i+3] = 0;
    }

    float* centroid_ss; 
                    
    initializing_centroids(centroid_ss,data_points,K);

    bool change = true;
    
    int iteration =0;

    
    while(change){   

     	iteration = iteration + 1; 

     	clustering(*data_point_cluster,N,K,data_points,centroid_ss);

//  		            cout<< "iteration "<<iteration<<" " <<check_count<<endl; 
 
 		for(int i0=0;i0<3*K;i0++){
   			cvector.push_back( centroid_ss[i0]);  }

        updating_centroids( N,K,centroid_ss,*data_point_cluster);

        if(check_count<1){
           change = false;
        }
        else{
           change = true;   
         }
      
        check_count =0;
       
    }

    num_iterations[0] = iteration;
    
    centroids_final(*centroids,cvector,iteration,K);
    
}