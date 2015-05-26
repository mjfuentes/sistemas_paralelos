#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <pthread.h>

double *A,*B,*C,*AB,*D;
int N,T;
int P = 0;
double res[24];
double pow (double base , double exponent);
float powf (float base  , float exponent);
long double powl (long double base, long double exponent);
double maxA = (double) -1;
double maxB = (double) -1;
double minB = (double) INT_MAX;
double minA = (double) INT_MAX;
double totA = 0;
double totB = 0;
double avgA,avgB,totA,totB,factor;
pthread_barrier_t barrier;
double *col;
int *pos;
double *col_res;
int *pos_res;
double timetick;

//Para calcular tiempo
double dwalltime(){
	double sec;
	struct timeval tv;
	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}

void merge(double *col, int *pos,double *col_res, int *pos_res, int inicio, int mid, int final)
{
	int i = inicio, j = mid + 1, k = inicio;
	while (i <= mid && j <= final)
	{
		if (col[i] > col[j])
		{
			col_res[k] = col[i];
			pos_res[k] = pos[i];
			i ++;
		}else
		{
			col_res[k] = col[j];
			pos_res[k] = pos[j];
			j ++;
		}
		k ++;
	}
	if (j > final){
		while (i<=mid){
			col_res[k] = col[i];
			pos_res[k] = pos[i];
			i ++;
			k ++;
		}
	} else if (i > mid){
		while (j<=final){
			col_res[k] = col[j];
			pos_res[k] = pos[j];
			j ++;
			k ++;
		}
	}
}

void mergeSort(double *col, int *pos,double *col_res, int *pos_res, int inicio, int final)
{
	if (inicio == final) return;
	int mid = (inicio + final)/2;
	mergeSort(col,pos,col_res,pos_res, inicio, mid);
	mergeSort(col,pos,col_res,pos_res,mid+1, final);
	merge(col,pos,col_res,pos_res, inicio, mid, final);
}

void imprimeMatriz(double *S,int N, int order,char comment[]){
	if (P){
		int i,j,I,J,despB;
		printf("%s",comment);
		printf("Contenido de la matriz: \n" );
		for (i=0; i<N; i++){
			for(j=0;j<N;j++){
				if (order == 1){
					printf("%f ",S[i*N+j]);
				}
				else {
					printf("%f ",S[j*N+i]);
				}
			}
			printf("\n ");
		};
		printf(" \n\n");
	}
}

void imprimeVector(double *S,int N,char comment[]){
	if (P){
		int i,j,I,J,despB;
		printf("%s",comment);
		printf("Contenido del vector: \n" );
		for(j=0;j<N;j++){
			printf("%f ",S[j]);
		}
		printf("\n\n");
	}
}


void *calcular(void *s){
	int principio,final,i,id,j,k,l;
	int *id_pointer;
	id_pointer=(int*) s;
	id=*id_pointer;
	principio=id*(N/T);
	final=principio+(N/T);
	double temp;int pos_temp;
	double *backup=(double*)malloc(sizeof(double)*N*N);

	// ********************************************
	// *************** ETAPA 1 ********************
	// ********************************************

	for (i=principio; i<final; i++){
		for(j=0;j<N;j++){
			for(k=0;k<N;k++){
				C[i*N+j] += A[i*N+k] * B[j*N+k];
			}
			C[i*N+j] = C[i*N+j] * D[j]; 

			if (A[i*N+j] > res[id*6]){
				res[id*6] = A[i*N+j];
			} 
			if (A[i*N+j] < res[id*6+1]){
				res[id*6+1] = A[i*N+j];
			}
			if (B[j*N+i] > res[id*6+2]){
				res[id*6+2] = B[j*N+i];
			} 
			if (B[j*N+i] < res[id*6+3]){
				res[id*6+3] = B[j*N+i];
			}
			res[id*6+4] += A[i*N+j];
			res[id*6+5] += B[j*N+i];
		}
	}

	pthread_barrier_wait (&barrier);
	
	// ********************************************
	// *************** ETAPA 2 ********************
	// ********************************************

	if (id == 0){
		imprimeMatriz(C,N,1,"Matriz resultante de A*B*D\n");
		for(i=0;i<T;i++){
			if (res[i*6] > maxA){
				maxA = res[i*6];
			} 
			if (res[i*6+1] < minA){
				minA = res[i*6+1];
			}
			if (res[i*6+2] > maxB){
				maxB = res[i*6+2];
			} 
			if (res[i*6+3] < minB){
				minB = res[i*6+3];
			}
			totA += res[i*6+4];
			totB += res[i*6+5];
		}

		avgA = (double) totA / (N * N);
		avgB = (double) totB / (N * N);
		double a = maxA - minA;
		double b = maxB - minB;
		factor = ((a*a)/avgA) * ((b*b)/avgB);
		if (P){
			printf("Factor: %f\n\n", factor);  
		}
	}

	pthread_barrier_wait (&barrier);

	for(i=principio;i<final;i++){
		for(j=0;j<N;j++){
			C[i*N+j] = C[i*N+j] * factor;
		}
	}
	
	pthread_barrier_wait (&barrier);

	if (id == 0){
		printf("Etapa 1 y 2 terminada en: %f\n",dwalltime() - timetick);
		timetick = dwalltime();  
		imprimeMatriz(C,N,1,"Matriz con factor aplicado\n");
	}

	// ********************************************
	// *************** ETAPA 3 ********************
	// ********************************************

	for(i=0;i<N;i++){
		for (l=principio;l<final;l++){
			col[l]=C[i+l*N];
			pos[l]=l;
		}

		mergeSort(col,pos,col_res,pos_res, principio, final-1);

		pthread_barrier_wait (&barrier);

		if (id % 2 == 0){
  
			int mid = principio + (N/T) - 1;
			merge(col_res,pos_res,col,pos,principio,mid, principio + 2*(N/T)-1);
		}

		pthread_barrier_wait (&barrier);

		if (id == 0){
			int mid = (N/2) - 1;
			merge(col,pos,col_res,pos_res,0,mid, N-1);
		}

		pthread_barrier_wait (&barrier);

		for (k=principio;k<final;k++){
			for (j=i;j<N;j++){
				backup[k*N+j] = C[pos_res[k]*N+j];
			}
		}

		pthread_barrier_wait (&barrier);

		for (k=principio;k<final;k++){
			for (j=i;j<N;j++){
				C[k*N+j] = backup[k*N+j];
			}
		}

		pthread_barrier_wait (&barrier);

	}

}

int main(int argc,char*argv[]){
	int id[T];
	int i,j,k;
	int check=1;
	T = 4;
	double timetick_total;
	pthread_t p_threads[T];
	pthread_attr_t attr;
	pthread_attr_init (&attr);

	if (argc < 3){
		printf("\n Faltan argumentos: Dimension de matriz e Impresion\n");
		return 0;
	}

	N=atoi(argv[1]);
	if (atoi(argv[2]) == 1){
		P = 1;
	} 
	P=atoi(argv[2]);

   //Aloca memoria para las matrices
	A=(double*)malloc(sizeof(double)*N*N);
	B=(double*)malloc(sizeof(double)*N*N);
	C=(double*)malloc(sizeof(double)*N*N);
	D=(double*)malloc(sizeof(double)*N);
	col=(double*)malloc(sizeof(double)*N);
	pos=(int*)malloc(sizeof(int)*N);
	col_res=(double*)malloc(sizeof(double)*N);
	pos_res=(int*)malloc(sizeof(int)*N);

   //Inicializa las matrices A y B en 1, C en diagonal
	srand(time(0));
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			A[i*N+j] = rand()%10+1;
			B[j*N+i] = rand()%10+1;
			C[i*N+j] = 0;
		}
		D[i] = rand()%10+1;
	}

	for (i=0;i<T;i++){
		res[(i*6)] = (double) INT_MIN;
		res[(i*6)+1] = (double) INT_MAX;
		res[(i*6)+2] = (double) INT_MIN;
		res[(i*6)+3] = (double) INT_MAX;
		res[(i*6)+4] = 0;
		res[(i*6)+5] = 0;
	}


	imprimeMatriz(A,N,1,"Matriz A\n");
	imprimeMatriz(B,N,0,"Matriz B\n");
	imprimeVector(D,N,"Vector D\n");

	pthread_barrier_init (&barrier, NULL, 4);

	timetick_total = dwalltime();  
	timetick = dwalltime();  

	for(i=0;i<T;i++){ 
		id[i] = i;
		pthread_create(&p_threads[i], &attr, calcular, (void *) &id[i]);
	}

	for (i=0; i< T; i++)
	{ 
		pthread_join(p_threads[i], NULL);
	}

	printf("Etapa 3 terminada en: %f\n",dwalltime() - timetick);
	imprimeMatriz(C,N,1,"Matriz ordenada\n");

	printf("Tiempo en segundos total: %f\n\n", dwalltime() - timetick_total);  

	free(A);
	free(B);
	free(C);
	free(D);
	free(pos);
	free(col);
	free(pos_res);
	free(col_res);

	return(0);
}