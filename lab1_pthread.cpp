#include "lab1_sequential.h"

#include<bits/stdc++.h> 
#include <iostream>
#include <pthread.h>
using namespace std;

int n;
int k;
int t;
int* for_loop;
int checkN = 0;
int* datapoints;
float* centroid;
float* centroid_ss;
int check_count =0;
vector<float> cvector;
int* datapoint_cluster;
pthread_mutex_t lock;


void initializing_centroids( float* &dup_centroid_ss,int*& datapoints,int k){ 
    
    dup_centroid_ss = (float * )malloc(k* 3 * sizeof(float)); 
    
    srand(0);

    int kar;
  
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
        dist = dist + pow(p1[i]-p2[i],2); 
    }

    return sqrt(dist);
}

void updating_centroids(int N,int k,float*& centroid_ss,int*& datapoint_cluster)
{
// averaging 
    for(int i =0;i<k;i++){
             int i_4 = 4*i;
             int i_3 = 3*i;
            centroid_ss[i_3+0] =    float(for_loop[i_4+0])/for_loop[i_4+3]; // has to be float recheck
            centroid_ss[i_3+1] =    float(for_loop[i_4+1])/for_loop[i_4+3];
            centroid_ss[i_3+2] =    float(for_loop[i_4+2])/for_loop[i_4+3];
    }

}

void* clustering_pt(void* tid){

    int tem_chk =0;
    int *id = (int*) tid;
    int len_per_thread = n/t;
    int start = (*id)*len_per_thread;
    int *shared_data =new int[4*k]();
    
    for(int i =start;i< start + len_per_thread;i++)
    {
        int* p1 ;
        int i_3 = 3*i;
        p1 = (int*)malloc(3*sizeof(int));

        p1[0] = datapoints[i_3+0];
        p1[1] = datapoints[i_3+1];
        p1[2] = datapoints[i_3+2];

        int ind = -1;
        float dis = 0;
        float min_dis = INT_MAX;

        for(int ii=0;ii<k;ii++){
            float* p2 ;
            p2 = (float*)malloc(3*sizeof(float));
            int ii_3 = 3*ii;
            p2[0] = centroid_ss[ii_3+0];
            p2[1] = centroid_ss[ii_3+1];
            p2[2] = centroid_ss[ii_3+2];
            dis = dist_p1_p2(p1,p2);

            // cout<<"id "<<ii <<" " <<u*id<<"  ind "<< centroid_ss[ii_3+0] << dis<<endl;

            if(dis<min_dis){
                min_dis = dis;
                ind = ii;                             
            } 
            free(p2);
        }
        // datapoint_cluster[4*i+0] = datapoints[3*i+0];
        // datapoint_cluster[4*i+1] = datapoints[3*i+1];
        // datapoint_cluster[4*i+2] = datapoints[3*i+2];
        // datapoint_cluster[4*i+3] = ind;
        // cout<<"id "<<*id<<"  ind "<< ind<<endl;
        if(datapoint_cluster[4*i+3] != ind){
            datapoint_cluster[4*i+3] = ind;
            tem_chk++;
         }
        int i_ind = 4*ind;
        shared_data[i_ind ] += p1[0];
        shared_data[i_ind +1] += p1[1];
        shared_data[i_ind +2] += p1[2];
        shared_data[i_ind +3] ++;

        free(p1);

    }

    if(*id == t-1 && n%t!=0){
        for(int i =start;i< start + len_per_thread;i++)
        {
            int* p1 ;
            int i_3 = 3*i;
            p1 = (int*)malloc(3*sizeof(int));

            p1[0] = datapoints[i_3+0];
            p1[1] = datapoints[i_3+1];
            p1[2] = datapoints[i_3+2];

            float min_dis = INT_MAX;
            float dis = 0;
            int ind = -1;

            for(int ii=0;ii<k;ii++){
                float* p2 ;
                p2 = (float*)malloc(3*sizeof(float));
                int ii_3 = 3*ii;
                p2[0] = centroid_ss[ii_3+0];
                p2[1] = centroid_ss[ii_3+1];
                p2[2] = centroid_ss[ii_3+2];
                dis = dist_p1_p2(p1,p2);

                if(dis<min_dis){
                    min_dis = dis;
                    ind = ii;
                } 
                free(p2);
            }

            if(datapoint_cluster[4*i+3] != ind){
                datapoint_cluster[4*i+3] = ind;
                tem_chk++;
            }
            int i_ind = 4*ind;
            shared_data[i_ind ] += p1[0];
            shared_data[i_ind +1] += p1[1];
            shared_data[i_ind +2] += p1[2];
            shared_data[i_ind +3] ++;
            free(p1);
        }
   }
    
    pthread_mutex_lock(&lock);

    check_count += tem_chk; 

     for(int i=0;i<4*k;i++){        
        for_loop[i] += shared_data[i];
       }
    // cout<<"entered "<<*id << " ck "<<check_count<<endl;// TODO
    pthread_mutex_unlock(&lock);

    return NULL;
} 
void centroids_final(float*& centroid,vector<float> cvector,int iteration,int K){

    for(int i =0;i<cvector.size();i++){
       
       centroid[i] = cvector[i];    

    }
}

void print_cent(float *& centroid_ss,int k){

    for(int i =0;i<3*k;i++){
        
        cout<<"print_cen "<< i <<" "<< centroid_ss[i]<<endl;
    }
}

void print_dpc(int *& datapoint_cluster,int N){

    for(int i =0;i<4;i++){
        
        cout<<"print_db "<< 4*i+3 <<" "<< datapoint_cluster[4*i+3]<<endl;
    }
}
  // float* centroid_ss;
void kmeans_pthread(int num_threads,
                    int N,
                    int K,
                    int* data_points,
                    int** data_point_cluster,
                    float** centroids,
                    int* num_iterations
                    )
{

    // vector<float> cvector;
    n=N;
    k=K;
    checkN = N;
    t=num_threads;
    *data_point_cluster = (int*)malloc(N*4*sizeof(int));

    for_loop = (int*)malloc(4*K*sizeof(int));

    for(int i =0;i<N;i++){
        (*data_point_cluster)[4*i+0] = data_points[3*i+0];
        (*data_point_cluster)[4*i+1] = data_points[3*i+1];
        (*data_point_cluster)[4*i+2] = data_points[3*i+2];
        (*data_point_cluster)[4*i+3] = 0;
    }
     
    datapoints = data_points;
    // float* centroid_ss; 

    centroid_ss = (float * )malloc(K* 3 * sizeof(float)); 
                    
    pthread_t some_thread[num_threads];

    int *tid = (int *) malloc (sizeof (int) * num_threads);

    pthread_mutex_init(&lock, NULL);

    initializing_centroids(centroid_ss,data_points,K);

    datapoint_cluster = *data_point_cluster;
 // print_cent(centroid_ss,K);
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
        

        for (int i = 0; i < num_threads; i++) {// i on heap?
            tid[i] = i;
            // args[i].id = &tid[i];
            pthread_create(&some_thread[i], NULL, clustering_pt, (void *)&tid[i]);
            }

        for (int i = 0; i < num_threads; i++){
            pthread_join(some_thread[i], NULL);
        }
                // cout<< "iteration "<<iteration<<" " <<check_count<<endl; 
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

    int i_v = cvector.size();

    *centroids = (float *)malloc(i_v* sizeof(float)); // thinkk abt size again 

    centroids_final(*centroids,cvector,iteration,K);
       
}