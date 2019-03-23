#include "lab1_sequential.h"

#include<bits/stdc++.h> 
#include <iostream>
using namespace std;

float* initializing_centroids(int k){
 
    float* dup_centroid_ss;

    dup_centroid_ss = (float *)malloc(k* 3 * sizeof(float)); 

    // for(int i =0;i<3*k;i++){
    //     dup_centroid_ss[i] = (float*) malloc(sizeof(float));
    // } 

    for(int i =0;i<k;i++){

        dup_centroid_ss[3*i][0] = (float)(rand()%2000 - 1000);
        dup_centroid_ss[3*i+1][0] = (float)(rand()%2000 - 1000);
        dup_centroid_ss[3*i+2][0] = (float)(rand()%2000 - 1000);

    }

    return dup_centroid_ss;

}


float dist_p1_p2(int* p1,float *p2){
 
    float dist =0;
    for(int i =0;i<3;i++){

        dist = dist + pow(p1[i]-p2[i],2); 
    }

    return sqrt(dist);
}

int* num_elem(int N,int k,int** datapoint_cluster){

    int* n1;
    n1 = (int*)malloc(k*sizeof(int));

    for(int i =0;i<N;i++){
     
       n1[datapoint_cluster[4*i+3][0]] = n1[datapoint_cluster[4*i+3][0]]+1;
    }

    return n1;
}

float** updating_centroids(int N,int k,float** centroid_ss,int** datapoint_cluster,int* num_elements)
{

    float** dup_centroid_ss;

    dup_centroid_ss = (float **)malloc(k* 3 * sizeof(float*)); 

    for(int i =0;i<3*k;i++){
        dup_centroid_ss[i] = (float*) malloc(sizeof(float));
    } 

    for(int i =0;i<N;i++){
        dup_centroid_ss[3*datapoint_cluster[4*i+3][0] + 0][0] = dup_centroid_ss[3*datapoint_cluster[4*i+3][0]+0][0]+ datapoint_cluster[4*i+0][0];
        dup_centroid_ss[3*datapoint_cluster[4*i+3][0] + 1][0] = dup_centroid_ss[3*datapoint_cluster[4*i+3][0]+1][0]+ datapoint_cluster[4*i+1][0];
        dup_centroid_ss[3*datapoint_cluster[4*i+3][0] + 2][0] = dup_centroid_ss[3*datapoint_cluster[4*i+3][0]+2][0]+ datapoint_cluster[4*i+2][0];
    }

// averaging 
    for(int i =0;i<k;i++){
        dup_centroid_ss[3*i+0][0] = dup_centroid_ss[3*i+0][0]/num_elements[i]; // has to be float recheck
        dup_centroid_ss[3*i+1][0] = dup_centroid_ss[3*i+1][0]/num_elements[i];
        dup_centroid_ss[3*i+2][0] = dup_centroid_ss[3*i+2][0]/num_elements[i];
    }

    return dup_centroid_ss;
}

int** clustering(int N,int k,int* datapoints,float** centroid_ss){

   int** datapoint_cluster;
   datapoint_cluster = (int**)malloc(N*4*sizeof(int*));

    for(int i =0;i<4*N;i++){
        datapoint_cluster[i] = (int*) malloc(sizeof(int));
    } 

   for(int i=0;i<N;i++){

        int* p1 ;
        p1 = (int*)malloc(3*sizeof(int));

        p1[0] = datapoints[3*i+0];
        p1[1] = datapoints[3*i+1];
        p1[2] = datapoints[3*i+2];

        float min_dis = INT_MAX;
        float dis = 0;
        int ind = -1;

        for(int ii=0;ii<k;ii++){
            float* p2 ;
            p2 = (float*)malloc(3*sizeof(float));
            p2[0] = centroid_ss[3*ii+0][0];
            p2[1] = centroid_ss[3*ii+1][0];
            p2[2] = centroid_ss[3*ii+2][0];
            dis = dist_p1_p2(p1,p2);

            if(dis<min_dis){
                min_dis = dis;
                ind = ii;
            } 
        }
       

        datapoint_cluster[4*i+0][0] = datapoints[3*i+0];
        datapoint_cluster[4*i+1][0] = datapoints[3*i+1];
        datapoint_cluster[4*i+2][0] = datapoints[3*i+2];
        datapoint_cluster[4*i+3][0] = ind;


    }

    return datapoint_cluster;

}

bool check_update(int N,int k,float** centroids,float** centroid_ss,int num_iterations){

    bool updated = false;
    for(int i =0;i<k;i++){
        if(centroids[3* (num_iterations-1) + i][0] != centroid_ss[3*i][0]){
           updated = true;
           break;
        }
         if(centroids[3* (num_iterations-1) + i+1][0] != centroid_ss[3*i+1][0]){
           updated = true;
           break;
        }
         if(centroids[3* (num_iterations-1) + i +2][0] != centroid_ss[3*i+2][0]){
           updated = true;
           break;
        }

    }

    return updated;
}




void kmeans_sequential(int N, int K, int* data_points, int** data_point_cluster, float** centroids, int* num_iterations )

{
   
    int* cluster_sizes; // size of each cluster

    float** centroid_ss;//output for current iteration
    centroid_ss = initializing_centroids(K); // initialising k centroids
    
    centroids = centroid_ss;

    bool change = true;
    
    int iteration =0;

    while(change){
       cout<<"entered0"<<endl;

     	iteration = iteration + 1; 
     	data_point_cluster = clustering(N,K,data_points,centroid_ss);

     	cout<<"entered1"<<endl;
     	
     	cluster_sizes = num_elem(N,K,data_point_cluster);
     	
     	cout<<"entered2"<<endl;
        
        centroid_ss = updating_centroids(N,K,centroid_ss,data_point_cluster,cluster_sizes);
        
        cout<<"entered3"<<endl;

        change = check_update(N,K,centroids,centroid_ss,iteration);
       
       cout<<"entered4"<<endl;

		if(change){
        	for(int i =0;i<K;i++){

            	centroids[3*K*(iteration-1) +3*i+0][0] = centroid_ss[3*i+0][0];
             	centroids[3*K*(iteration-1) +3*i+1][0] = centroid_ss[3*i+1][0];
              	centroids[3*K*(iteration-1) +3*i+2][0] = centroid_ss[3*i+2][0];
        	}        
     	}
     	change = false;
    }


 
     for(int pp=0;pp<N;pp++){
        // cout<< "loop "<< "pp0  " << pp<<" _ " << data_point_cluster[4*pp + 0][0]<< endl;
        // cout<< "loop "<< "pp1 " << pp<<" _ "<< data_point_cluster[4*pp + 1][0]<< endl;
        // cout<< "loop "<< "pp2 " << pp<<" _ "<< data_point_cluster[4*pp + 2][0]<< endl;
        cout<< "loop "<< "pp3 " << pp<<" _ "<< data_point_cluster[4*pp + 3][0]<< endl;


     }
     cout<<"entered5"<<endl;
    num_iterations[0] = iteration;
    cout<<"entered6"<<endl;
   
   
}