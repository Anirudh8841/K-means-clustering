#include "lab1_sequential.h"

#include<bits/stdc++.h> 
#include <iostream>
#include <pthread.h>
using namespace std;

pthread_mutex_t lock;

struct arg_struct {
int N;
int k;
int t;
float* centroid_ss;
float* centroid;
int* datapoint_cluster;
int* datapoints;
int* num_elements_t;
int* num_elements;
int* id;
};



vector<float> cvector;


void initializing_centroids( float* &dup_centroid_ss,int k){
 
    dup_centroid_ss = (float * )malloc(k* 3 * sizeof(float)); 

    for(int i =0;i<3*k;i++){

        dup_centroid_ss[i] = (rand()%2000 - 1000);

    }
}

float dist_p1_p2(int* &p1,float *&p2){
 
    float dist =0;
    for(int i =0;i<3;i++){

        dist = dist + pow(p1[i]-p2[i],2); 
    }

    return sqrt(dist);
}

void* num_elem_pt(void* arguements){
 

    struct arg_struct *args = (struct arg_struct *) arguements;

    float* & centroid_ss = args->centroid_ss;
    int*& datapoints = args->datapoints;
    int*& datapoint_cluster = args->datapoint_cluster;
    int N = args->N;
    int k = args->k;
    int t = args->t;
    int* n1 = args->num_elements_t;
    int* orig_size = args->num_elements;
    int *id = args -> id;


    int len_per_thread = N/t;
    int start = (*id)*len_per_thread;

for(int i =0;i< k;i++)
    {
      n1[i] =0;
    }

for(int i =start;i< start + len_per_thread;i++)
    {
      n1[datapoint_cluster[4*i+3]] += 1;
    }

 if((*id == t-1)&&(N%t!=0)){

    for(int i =start + len_per_thread;i<N;i++)
    {
      n1[datapoint_cluster[4*i+3]] += 1;
    }
 }
        

pthread_mutex_lock(&lock);
for(int i =0;i<k;i++){
    orig_size[i] += n1[i];
}
pthread_mutex_unlock(&lock);

  return NULL;
}

void num_elem(int* &n1,int N,int k,int* &datapoint_cluster){

    for(int i =0;i<k;i++){
     
       n1[i] = 0;
    }

     
    for(int i =0;i<N;i++){
     
       n1[datapoint_cluster[4*i+3]] = n1[datapoint_cluster[4*i+3]]+1;
    }

}





// void* updating_centroids_pt(void* arguements){

//  struct arg_struct *args = (struct arg_struct *) arguements;

//     float* & centroid_ss = args->centroid_ss;
//     int*& datapoints = args->datapoints;
//     int*& datapoint_cluster = args->datapoint_cluster;
//     int N = args->N;
//     int k = args->k;
//     int t = args->t;
//     int *id = args -> id;


//     int len_per_thread = N/t;
//     int start = (*id)*len_per_thread;
    

  
//   for(int i =start;i<start + len_per_thread ;i++){
//         centroid_ss[3* (datapoint_cluster[4*i+3]) + 0] += datapoint_cluster[4*i+0];
//         centroid_ss[3* (datapoint_cluster[4*i+3]) + 1] += datapoint_cluster[4*i+1];
//         centroid_ss[3* (datapoint_cluster[4*i+3]) + 2] += datapoint_cluster[4*i+2];
//     }



//  return NULL;
// }

void updating_centroids(int N,int k,float*& centroid_ss,int*& datapoint_cluster,int*& num_elements)
{

    for(int i =0;i<k;i++){
        centroid_ss[3* i + 0] = 0;
        centroid_ss[3* i + 1] = 0;
        centroid_ss[3* i + 2] = 0;
    }

    for(int i =0;i<N;i++){
        centroid_ss[3* (datapoint_cluster[4*i+3]) + 0] += datapoint_cluster[4*i+0];
        centroid_ss[3* (datapoint_cluster[4*i+3]) + 1] += datapoint_cluster[4*i+1];
        centroid_ss[3* (datapoint_cluster[4*i+3]) + 2] += datapoint_cluster[4*i+2];
    }
    
// averaging 
    for(int i =0;i<k;i++){
        if(num_elements[i]!=0){
            centroid_ss[3*i+0] =    centroid_ss[3*i+0]/num_elements[i]; // has to be float recheck
            centroid_ss[3*i+1] =    centroid_ss[3*i+1]/num_elements[i];
            centroid_ss[3*i+2] =    centroid_ss[3*i+2]/num_elements[i];
        }
        else{
            centroid_ss[3*i+0] = 0;
            centroid_ss[3*i+1] = 0;
            centroid_ss[3*i+2] = 0;
        }
    }

}
 
void* clustering_pt(void* arguements){

 struct arg_struct *args = (struct arg_struct *) arguements;

    float* & centroid_ss = args->centroid_ss;
    int*& datapoints = args->datapoints;
    int*& datapoint_cluster = args->datapoint_cluster;
    int N = args->N;
    int k = args->k;
    int t = args->t;
    int *id = args -> id;


    int len_per_thread = N/t;
    int start = (*id)*len_per_thread;


    
    for(int i =start;i< start + len_per_thread;i++)
    {
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
            p2[0] = centroid_ss[3*ii+0];
            p2[1] = centroid_ss[3*ii+1];
            p2[2] = centroid_ss[3*ii+2];
            dis = dist_p1_p2(p1,p2);

            if(dis<min_dis){


                min_dis = dis;
                ind = ii;
            } 
        }
       

        datapoint_cluster[4*i+0] = datapoints[3*i+0];
        datapoint_cluster[4*i+1] = datapoints[3*i+1];
        datapoint_cluster[4*i+2] = datapoints[3*i+2];
        datapoint_cluster[4*i+3] = ind;


    }

     if(*id == t-1){
        if(N%t!=0){
              for(int i =start + len_per_thread; i< N ;i++)
               {
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
                        p2[0] = centroid_ss[3*ii+0];
                        p2[1] = centroid_ss[3*ii+1];
                        p2[2] = centroid_ss[3*ii+2];
                        dis = dist_p1_p2(p1,p2);

                        if(dis<min_dis){
                            min_dis = dis;
                            ind = ii;
                        } 
                    }
       
                    datapoint_cluster[4*i+0] = datapoints[3*i+0];
                    datapoint_cluster[4*i+1] = datapoints[3*i+1];
                    datapoint_cluster[4*i+2] = datapoints[3*i+2];
                    datapoint_cluster[4*i+3] = ind;
                        }
                }
    }
    
    return NULL;
} 


void clustering(int*& datapoint_cluster,int N,int k,int*& datapoints,float*& centroid_ss){

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
            p2[0] = centroid_ss[3*ii+0];
            p2[1] = centroid_ss[3*ii+1];
            p2[2] = centroid_ss[3*ii+2];
            dis = dist_p1_p2(p1,p2);

            if(dis<min_dis){


                min_dis = dis;
                ind = ii;
            } 
        }
       

        datapoint_cluster[4*i+0] = datapoints[3*i+0];
        datapoint_cluster[4*i+1] = datapoints[3*i+1];
        datapoint_cluster[4*i+2] = datapoints[3*i+2];
        datapoint_cluster[4*i+3] = ind;


    }

}




bool check_update(int N,int k,vector<float> cvector,float*& centroid_ss,int num_iterations){

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

    }

    return updated;
}




void* centroids_final_pt(void* arguements){

 struct arg_struct *args = (struct arg_struct *) arguements;

    float* & centroid = args->centroid;
    
    int N = args->N;
    int k = args->k;
    int t = args->t;
    int *id = args -> id;


    int len_per_thread = cvector.size()/t;
    int start = (*id)*len_per_thread;

    for(int i =start;i<start+len_per_thread;i++){
       
       centroid[i] = cvector[i];    

    }
    if((cvector.size())%t!=0 && *id ==t-1 ){
        for(int i =start+len_per_thread;i<cvector.size();i++){
       
          centroid[i] = cvector[i];     

            }
    }

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

    for(int i =0;i<4*N;i++){
        
        cout<<"print_db "<< i <<" "<< datapoint_cluster[i]<<endl;
    }
}

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

    *data_point_cluster = (int*)malloc(N*4*sizeof(int));

    int* cluster_sizes; // size of each cluster
    int* cluster_sizes_t[num_threads];

    cluster_sizes = (int*)malloc(K*sizeof(int));
    for(int i =0;i<num_threads;i++){
            cluster_sizes_t[i] = (int*)malloc(K*sizeof(int));
    }

    float* centroid_ss; 

    centroid_ss = (float * )malloc(K* 3 * sizeof(float)); 
                    
    pthread_t some_thread[num_threads];

    int *tid = (int *) malloc (sizeof (int) * num_threads);
    // for(int i =0;i<K;i++){
    //  tem_n1[i] = (int *) malloc (sizeof (int) * num_threads);
    // }
    
    pthread_mutex_init(&lock, NULL);

    struct arg_struct args[num_threads];

  for(int i =0;i<num_threads;i++){


    args[i].N = N;

    args[i].k = K;
    args[i].t = num_threads;
    args[i].centroid = *centroids;
    args[i].centroid_ss = centroid_ss;
    args[i].datapoint_cluster = *data_point_cluster;
    args[i].datapoints = data_points;
    args[i].num_elements_t = cluster_sizes_t[i] ;
    
}

    initializing_centroids(centroid_ss,K);

    bool change = true;
    
    int iteration =0;
  
    while(change){
        
       // cout<<"entered0"<<endl;

        iteration = iteration + 1; 
           
    for (int i = 0; i < num_threads; i++) {// i on heap?
            tid[i] = i;
            args[i].id = &tid[i];
            pthread_create(&some_thread[i], NULL, clustering_pt, (void *)&args[i]);
            }

    for (int i = 0; i < num_threads; i++){
        pthread_join(some_thread[i], NULL);}

        // clustering(*data_point_cluster,N,K,data_points,centroid_ss);

        // print_cent(centroid_ss,K);
        // print_dpc(*data_point_cluster,N); 
        if(change){       
            for(int i0=0;i0<3*K;i0++){
            cvector.push_back( centroid_ss[i0]);
         }
        }

    for(int i=0;i<K;i++){cluster_sizes[i] =0;}

    for (int i = 0; i < num_threads; i++) {// i on heap?
            tid[i] = i;
            args[i].id = &tid[i];
            args[i].num_elements = cluster_sizes;
            pthread_create(&some_thread[i], NULL, num_elem_pt, (void *)&args[i]);
            }

    for (int i = 0; i < num_threads; i++){
        pthread_join(some_thread[i], NULL);}

        updating_centroids( N,K,centroid_ss,*data_point_cluster,cluster_sizes);


        change = check_update(N,K,cvector,centroid_ss,iteration);
        // change = false;
       

    }

    num_iterations[0] = iteration;
    *centroids = (float *)malloc((iteration+1)*K* 3 * sizeof(float)); 

  for (int i = 0; i < num_threads; i++) {// i on heap?
            tid[i] = i;
            args[i].id = &tid[i];
            args[i].centroid = *centroids;
            pthread_create(&some_thread[i], NULL, clustering_pt, (void *)&args[i]);
            }

    for (int i = 0; i < num_threads; i++){
        pthread_join(some_thread[i], NULL);}

    // centroids_final(*centroids,cvector,iteration,K);
    
   
   
}