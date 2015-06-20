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
    int i,j,I,J,despB;
    printf("%s",comment);
    printf("Contenido del vector: \n" );
    for(j=from;j<N;j++){
        printf("%i ",S[j]);
    }
    printf("\n\n");
}

int main(int argc, char *argv[]) 
{ 
    int rank, size, N, i, max_step;
    int *A, *B, *U;
    int *lcl, *tmp;
    double timetick;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    size = size -1;
    if(argc < 2){
        printf("Se debe indicar el tamaño del vector\n");
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    else {
        N = atoi(argv[1]);
    }

    if(rank == 0){//el proceso 0 genera un vector desordenado.
        printf("Cantidad de workers: %i\n",size);
        printf("tamaño de bloque por worker: %i\n", N/size);
        A = (int*)malloc(sizeof(double)*N);
        B = (int*)malloc(sizeof(double)*N);
        srand(time(0));
    	for(i=0;i<N;i++){
            A[i] = rand()%1000+1;
        }
        //imprimeVector(A,0,N,"");
        timetick = dwalltime();
    }
    U = (int*)malloc(sizeof(double)*size);
    lcl = (int*)malloc(sizeof(double)*N);
    tmp = (int*)malloc(sizeof(double)*N);
    int ok = 1;
    int receiving = 1;
    int tag, step, out, stage,done;
    max_step = size;
    stage = 0;
    tmp = (int*)malloc(sizeof(double)*(N/size)*2);
    out=0; step=0; done=0;
    tag = 0;
    int count = 0;
    int stage_size = (N/size);
    while(ok == 1)
    {
        if(rank != 0){ // Cada uno de los workers
            MPI_Request request;
            MPI_Send(&(lcl[0]),0,MPI_INT,0,0,MPI_COMM_WORLD);
            MPI_Status status; 
            MPI_Recv(&(lcl[0]),N,MPI_INT,0,MPI_ANY_TAG,MPI_COMM_WORLD,&status); 
            MPI_Get_count(&status,MPI_INT,&count);
            tag = status.MPI_TAG;
            if (tag != 0){
                if (tag == 3){
                    // wait recibido
                }
                if (tag == 1){
                    // stage cero es sort
                    // printf("Proceso %i recibe sort\n",rank);
                    // imprimeVector(lcl,0,count,"");
                    sort(lcl, N/size);
                    // printf("Proceso %i devuelve sort\n",rank);
                    // imprimeVector(lcl,0,count,"");
                    MPI_Send(&(lcl[0]),N/size,MPI_INT,0,1,MPI_COMM_WORLD);
                    //printf("Proceso %i devuelve sort\n",rank);
                }
                else if (tag == 2){
                    // stage mayor a cero es merge
                    // printf("Proceso %i recibe merge\n",rank);
                    // imprimeVector(lcl,0,count,"");
                    int half = count / 2;
                    merge(lcl,tmp, 0, half - 1, count - 1);
                    // printf("Proceso %i devuelve merge\n",rank);
                    // imprimeVector(lcl,0,count,"");
                    MPI_Send(&(lcl[0]),count,MPI_INT,0,2,MPI_COMM_WORLD);
                }
            }
            else {
                 ok = 0;
            }
        }
        else // master
        {
            MPI_Status status;
            MPI_Recv(&(lcl[0]),N,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
            tag = status.MPI_TAG;
            if (tag == 0){
            MPI_Get_count(&status,MPI_INT,&count);
                //printf("MASTER recibe request de %i\n", status.MPI_SOURCE);
                MPI_Request request;
                if (receiving == 1){
                    if (stage == 0){
                        U[status.MPI_SOURCE] = step*(N/size);
                        MPI_Send(&(A[step*(N/size)]),N/size,MPI_INT,status.MPI_SOURCE,1,MPI_COMM_WORLD);
                        //printf("MASTER manda sort a %i\n", status.MPI_SOURCE);
                        //imprimeVector(A,step*(N/size),step*(N/size) + (N/size),"");
                        step = step + 1;
                    }
                    else {
                        U[status.MPI_SOURCE] = step*stage_size;
                        MPI_Send(&(A[step*stage_size]),stage_size,MPI_INT,status.MPI_SOURCE,2,MPI_COMM_WORLD);
                        //printf("MASTER manda merge a %i\n", status.MPI_SOURCE);
                        step = step + 1;
                    }
                    if (step == max_step){
                        receiving = 2;
                    }
                }
                else if (receiving == 2){
                    //printf("MASTER manda wait a %i\n", status.MPI_SOURCE);
                    MPI_Send(&(lcl[0]),0,MPI_INT,status.MPI_SOURCE,3,MPI_COMM_WORLD);
                }
                else if (receiving == 0){
                    //printf("MASTER manda close a %i\n", status.MPI_SOURCE);
                    MPI_Send(&(lcl[0]),0,MPI_INT,status.MPI_SOURCE,0,MPI_COMM_WORLD);
                    out++;
                    if (out >= size){
                        ok = 0;
                        printf("Tiempo en segundos total: %f\n\n", dwalltime() - timetick);  
                        //imprimeVector(A,0,N,"");
                    }
                }
            }
            else {
                if (tag == 1){
                    //printf("MASTER recibe sort de %i\n",status.MPI_SOURCE);
                    int temp = U[status.MPI_SOURCE];
                    save(A,temp,temp + (N/size) - 1,lcl);
                    //imprimeVector(A,temp,temp + (N/size),"");
                }
                else {
                    //printf("MASTER recibe merge de %i\n",status.MPI_SOURCE);
                    int temp = U[status.MPI_SOURCE];
                    //imprimeVector(lcl,0,stage_size,"");
                    //printf("MASTER guarda el merge en el indice: %i\n",temp);
                    save(A,temp,stage_size - 1, lcl);
                }
                done++;
                if ((receiving == 2) && (done == max_step)){
                    if (max_step > 1){
                        step = 0;
                        done = 0;
                        receiving = 1;
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
 
    MPI_Finalize();
    return 0; 
}