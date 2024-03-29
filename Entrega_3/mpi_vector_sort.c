#include <malloc.h>
#include <math.h>
#include <stdio.h>
#include <mpi.h>
#include <limits.h>
#include <stdlib.h>


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

void mergeSort(int *col,int *col_res, int inicio, int final)
{
    if (inicio == final) return;
    int mid = (inicio + final)/2;
    mergeSort(col,col_res, inicio, mid);
    mergeSort(col,col_res,mid+1, final);
    merge(col,col_res, inicio, mid, final);
}


double dwalltime(){
    double sec;
    struct timeval tv;
    gettimeofday(&tv,NULL);
    sec = tv.tv_sec + tv.tv_usec/1000000.0;
    return sec;
}

void save(int *B, int from, int to, int *A){
    int i;
    for (i=from; i<= to; i++){
        B[i] = A[i-from];
    }
}

void imprimeVector(int *S,int from, int N,char comment[]){
    int j;
    printf("%s",comment);
    printf("Contenido del vector: \n" );
    for(j=from;j<N;j++){
        printf("%i ",S[j]);
    }
    printf("\n\n");
}

int main(int argc, char *argv[]) 
{ 
    int rank, size, N, i, max_step, M;
    int tag, step, out, stage, done, count, ok, receiving;
    int *A, *U, *lcl, *tmp;
    double timetick;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    //size = size -1;
    out=0;step=0; done=0;stage = 0;ok=1;receiving=1;
    MPI_Status status;
    MPI_Request request;
    
    if(argc < 2){
        printf("Se debe indicar el tamaño del vector\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    else {
        N = atoi(argv[1]);
        if (N % size != 0){
            printf("Tamaño de vector inadecuado\n");
            MPI_Abort(MPI_COMM_WORLD,1);
        }
    }
    if(argc < 3){
        printf("Se debe indicar el tamaño de la porcion\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    else {
        M = atoi(argv[2]);
        if (N % M != 0){
            printf("Tamaño de porcion inadecuado\n");
            MPI_Abort(MPI_COMM_WORLD,1);
        }
        max_step = N/M; 
    }

    if(rank == 0){
        printf("Cantidad de workers: %i\n",size);
        A = (int*)malloc(sizeof(double)*N);
        U = (int*)malloc(sizeof(int)*N);
        srand(time(0));
    	for(i=0;i<N;i++){
            A[i] = rand()%1000+1;
        }
        timetick = dwalltime();
    }

    //MPI_Barrier(MPI_COMM_WORLD);
    
    lcl = (int*)malloc(sizeof(int)*N);
    tmp = (int*)malloc(sizeof(int)*N);
    int stage_size = M;
    while(ok == 1)
    {
        if(rank != 0){ 
            MPI_Send(&(lcl[0]),0,MPI_INT,0,0,MPI_COMM_WORLD);
            MPI_Recv(&(lcl[0]),N,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status); 
            MPI_Get_count(&status,MPI_INT,&count);
            tag = status.MPI_TAG;
            if (tag == 1){
                mergeSort(lcl,tmp,0, M);
                MPI_Send(&(lcl[0]),M,MPI_INT,0,1,MPI_COMM_WORLD);
            }
            else if (tag == 2){
                int half = count/2;
                merge(lcl,tmp, 0, half - 1, count - 1);
                MPI_Send(&(lcl[0]),count,MPI_INT,0,1,MPI_COMM_WORLD);
            }
            else {
                 break;
            }
        }
        else
        {
            if (receiving == 1 || receiving == 0){
                MPI_Recv(&(lcl[0]),N,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
            }
            else if (receiving == 2){
                MPI_Recv(&(lcl[0]),N,MPI_INT,MPI_ANY_SOURCE,1,MPI_COMM_WORLD,&status);
            }
            tag = status.MPI_TAG;
            if (tag == 0){
                MPI_Request request;
                if (receiving == 1){
                    if (stage == 0){
                        U[status.MPI_SOURCE] = step*M;
                        MPI_Send(&(A[step*M]),M,MPI_INT,status.MPI_SOURCE,1,MPI_COMM_WORLD);
                    }
                    else {
                        U[status.MPI_SOURCE] = step*stage_size;
                        MPI_Send(&(A[step*stage_size]),stage_size,MPI_INT,status.MPI_SOURCE,2,MPI_COMM_WORLD);
                    }
                    step++;
                    if (step == max_step){
                        receiving = 2;
                    }
                }
                else if (receiving == 0){
                    MPI_Send(&(lcl[0]),0,MPI_INT,status.MPI_SOURCE,3,MPI_COMM_WORLD);
                    out++;
                    if (out == size - 1){
                        //imprimeVector(A,0,N,"");
                        printf("Tiempo en segundos total: %f\n\n", dwalltime() - timetick);  
                        break;
                    }
                }
            }
            else {
                MPI_Get_count(&status,MPI_INT,&count);
                int temp = U[status.MPI_SOURCE];
                save(A,temp,temp + count - 1, lcl);
                done++;
                if ((receiving == 2) && (done >= max_step)){
                    if (max_step > 1){
                        step = 0;done = 0;receiving = 1;
                        max_step = max_step / 2;
                        stage++;
                        stage_size = stage_size * 2;
                    }
                    else {
                        receiving = 0;
                    }
                }
            }
            if ((step > 0) && (step % (size-1) == 0)){
                if (stage == 0){   
                    mergeSort(A,tmp,step*M, step*M + M -1);
                    
                }
                else {
                    int half = stage_size/2;
                    merge(A,tmp, step*stage_size, step*stage_size + half - 1, step*stage_size + stage_size - 1);
                }
                step++;
                done++;
                if (step == max_step){
                    receiving = 2;
                     if (done >= max_step){
                        if (max_step > 1){
                            step = 0;done = 0;receiving = 1;
                            max_step = max_step / 2;
                            stage++;
                            stage_size = stage_size * 2;
                        }
                        else {
                            receiving = 0;
                        }
                    }
                }
            }       
        }
    }
    MPI_Finalize();
    return 0; 
}