#include <algorithm>
#include <vector>
#include "mpi.h" 
#include <iostream>
using namespace std;
 

void merge(double *col,double *col_res, int inicio, int mid, int final)
{
    int i = inicio, j = mid + 1, k = inicio;
    while (i <= mid && j <= final)
    {
        if (col[i] > col[j])
        {
            col_res[k] = col[i];
            i ++;
        }else
        {
            col_res[k] = col[j];
            j ++;
        }
        k ++;
    }
    if (j > final){
        while (i<=mid){
            col_res[k] = col[i];
            i ++;
            k ++;
        }
    } else if (i > mid){
        while (j<=final){
            col_res[k] = col[j];
            j ++;
            k ++;
        }
    }
    for (i=inicio;i<=final;i++){
        col[i]=col_res[i];
    }
}

int main(int argc, char *argv[]) 
{ 
    int rank, size, N;
    int *A;
    int *lcl;
 
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
 
    if( size % 2 != 0 ){
    	printf("El numero de procesos debe ser par");
    }
 
    if(argc < 2){
        printf("Se debe indicar el tamaño del vector");
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    else {
        N = atoi(argv[1]);
    }
 
    if(rank == 0){//el proceso 0 genera un vector desordenado.
    	for(i=0;i<N;i++){
            A[i] = rand()%10+1;
        }
    }

    lcl = (int*)malloc(sizeof(double)*(N/size)*2);
    bool ok = true;
    int tag, step, out, stage,max_stage;
    out=0; step=0;
    max_stage = //logaritmo en base 2 del size
    //Comienza el ordenado
    while(ok)
    {
        if(rank != 0){ // Cada uno de los workers
            MPI_Request request;
            MPI_Irecv(&(lcl[0]),N/size*2,MPI_INT,0,tag,MPI_COMM_WORLD,&request); 
            MPI_Send(&(lcl[0]),0,MPI_INT,0,-1,MPI_COMM_WORLD);
            MPI_Status status; 
            int MPI_Wait(&request, &status); 
            int stage = tag / size;
            int pos = tag % size;
            if (stage == 0){
                // stage cero es sort
                sortArray(&lcl, 0, N/size);
                MPI_Send(&(lcl[0]),N/size,MPI_INT,0,tag,MPI_COMM_WORLD);
            }
            else if (stage > 0){
                // stage mayor a cero es merge
                int half = N / size * (stage - 1);
                merge(&lcl, 0, half - 1, half*2 - 1);
                MPI_Send(&(lcl[0]),N/size*2,MPI_INT,0,tag,MPI_COMM_WORLD);
            }
            else {
                 ok = false;
            }
        }
        else // master
        {
            MPI_Status status;
            MPI_Recv(&(lcl[0]),N/size*2,MPI_INT,MPI_ANY_SOURCE,tag,MPI_COMM_WORLD,&status);
            if (tag == 0){
                if (stage <= max_stage){
                    if (stage == 0){
                        MPI_Isend(&(A[step*(N/size)]),N/size,MPI_INT,status.MPI_SOURCE,step,MPI_COMM_WORLD);
                        step = step + (N/size);
                    }
                    else {
                        MPI_Isend(&(A[step_2*(N/size*2)]),N/size*2,MPI_INT,status.MPI_SOURCE,step,MPI_COMM_WORLD);
                        step = step + 1;
                    }
                    if (step == N){
                        step = 0;
                        stage++;
                    }
                }
                else {
                    MPI_Isend(&(lcl[0]),0,MPI_INT,status.MPI_SOURCE,-1,MPI_COMM_WORLD);
                    out++;
                    if (out >= size){
                        ok = false;
                    }
                }
            }
            else {
                if (tag < size){
                    
                }
                else if (tag >= size && tag < (size + size/2)){

                }
            }
        
        }
    }
 
    MPI_Finalize();
    return 0; 
}