#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <mpi.h>
 

void merge(int *col,int *col_res, int inicio, int mid, int final)
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

void sort(int *C, int N){
    int l,j;
    for (l=0;l<N-1;l++){
        for(j=0;j<N-l;j++){
            if (C[j] < C[j+1]){
                int temp = C[j];
                C[j] = C[j+1];
                C[j+1] = temp;
            }
        }
    }
}

void save(int *B, int from, int to, int *A){
    int i;
    for (i=from; i<= to; i++){
        B[i] = A[i];
    }
}

int main(int argc, char *argv[]) 
{ 
    int rank, size, N, i, max_step;
    int *A, *B;
    int *lcl, *tmp;
 
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

    A = (int*)malloc(sizeof(double)*N);
    B = (int*)malloc(sizeof(double)*N);
    lcl = (int*)malloc(sizeof(double)*N);
    tmp = (int*)malloc(sizeof(double)*N);
    int ok = 1;
    int receiving = 1;
    int tag, step, out, stage,max_stage,done;
    max_stage = //logaritmo en base 2 del size
    max_step = size;
    tmp = (int*)malloc(sizeof(double)*(N/size)*2);
    out=0; step=0; done=0;
    tag = 0;

    //Comienza el ordenado
    while(ok == 1)
    {
        printf("Proceso %i inicia",rank);
        if(rank != 0){ // Cada uno de los workers
            MPI_Request request;
            MPI_Irecv(&(lcl[0]),N/size*2,MPI_INT,0,tag,MPI_COMM_WORLD,&request); 
            MPI_Send(&(lcl[0]),0,MPI_INT,0,0,MPI_COMM_WORLD);
            printf("Proceso %i pide datos",rank);
            MPI_Status status; 
            MPI_Wait(&request, &status);
            int stage = tag / size;
            int pos = tag % size;
            if (tag != N){
                if (stage == 0){
                    // stage cero es sort
                    sort(lcl, N/size);
                    MPI_Send(lcl,N/size,MPI_INT,0,tag,MPI_COMM_WORLD);
                }
                else if (stage > 0){
                    // stage mayor a cero es merge
                    int half = N / size * (stage - 1);
                    merge(lcl,tmp, 0, half - 1, half*2 - 1);
                    MPI_Send(lcl,N/size*2,MPI_INT,0,tag,MPI_COMM_WORLD);
                }
            }
            else {
                 ok = 0;
            }
        }
        else // master
        {
            MPI_Status status;
            MPI_Recv(lcl,N/size*2,MPI_INT,MPI_ANY_SOURCE,tag,MPI_COMM_WORLD,&status);
            if (tag == 0){
                MPI_Request request;
                if (receiving == 1){
                    if (stage == 0){
                        MPI_Isend(&(A[step*(N/size)]),N/size,MPI_INT,status.MPI_SOURCE,step,MPI_COMM_WORLD,&request);
                        step = step + 1;
                    }
                    else {
                        int tag = stage * size + step;
                        MPI_Isend(&(A[step*(N/size*2)]),N/size*2,MPI_INT,status.MPI_SOURCE,tag,MPI_COMM_WORLD,&request);
                        step = step + 1;
                    }
                    if (step == max_step){
                        if (max_step > 1){
                            step = 0;
                            max_step = max_step / 2;
                            stage++;
                        }
                        else {
                            receiving = 0;
                        }
                    }

                }
                else {
                    MPI_Isend(&(lcl[0]),0,MPI_INT,status.MPI_SOURCE,N,MPI_COMM_WORLD,&request);
                    out++;
                    if (out >= size){
                        ok = 0;
                    }
                }
            }
            else {
                int stage = tag / size;
                int pos = tag % size;
                if (stage == 0){
                    int temp = pos * (N/size);
                    save(A,temp,temp + (N/size),lcl);
                    done++;
                }
                else {
                    int aux = tag - size;
                    int temp = aux * (N/size) * 2;
                    save(A,temp, temp + (N/size) * 2, lcl);
                }
            }
        
        }
    }
 
    MPI_Finalize();
    return 0; 
}