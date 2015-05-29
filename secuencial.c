#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>

//Para calcular tiempo
double dwalltime(){
	double sec;
	struct timeval tv;
	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}

double pow (double base , double exponent);
float powf (float base  , float exponent);
long double powl (long double base, long double exponent);

void imprimeMatriz(double *S,int N, int order){
	int i,j,I,J,despB;

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

void imprimeVector(double *S,int N,char comment[]){
		int i,j,I,J,despB;
		printf("%s",comment);
		printf("Contenido del vector: \n" );
		for(j=0;j<N;j++){
			printf("%f ",S[j]);
		}
		printf("\n\n");
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
	for (i=inicio;i<=final;i++){
		col[i]=col_res[i];
		pos[i]=pos_res[i];
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

int main(int argc,char*argv[]){
	double *A,*B,*C,*AB,*D, *CO;
	int i,j,k,N;
	int check=1;
	double timetick;

	if (argc < 2){
		printf("\n Falta un argumento:: N dimension de la matriz \n");
		return 0;
	}

	N=atoi(argv[1]);

   //Aloca memoria para las matrices
	A=(double*)malloc(sizeof(double)*N*N);
	B=(double*)malloc(sizeof(double)*N*N);
	C=(double*)malloc(sizeof(double)*N*N);
	D=(double*)malloc(sizeof(double)*N);
	CO=(double*)malloc(sizeof(double)*N*N);

   //Inicializa las matrices A y B en 1, D en diagonal
	srand(time(0));
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			A[i*N+j] = rand()%10+1;
			B[j*N+i] = rand()%10+1;
		}
		D[i] = rand()%10+1;
	}


	// ********************************************
	// *************** ETAPA 1 ********************
	// ********************************************
	timetick = dwalltime();

	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			C[i*N+j] = 0;
			for(k=0;k<N;k++){
				C[i*N+j] += A[i*N+k] * B[j*N+k];
			}
			C[i*N+j] = C[i*N+j] * D[j]; 
		}
	}   


	// ********************************************
	// *************** ETAPA 2 ********************
	// ********************************************

	int maxA,minA,maxB,minB,totA,totB;
	double factor, avgA, avgB;
	maxA = INT_MIN;
	maxB = INT_MIN;
	minB = INT_MAX;
	minA = INT_MAX;
	totA = 0;
	totB = 0;


	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			if (A[i*N+j] > maxA){
				maxA = A[i*N+j];
			} 
			if (A[i*N+j] < minA){
				minA = A[i*N+j];
			}
			if (B[j*N+i] > maxB){
				maxB = B[j*N+i];
			} 
			if (B[j*N+i] < minB){
				minB = B[j*N+i];
			}
			totA += A[i*N+j];
			totB += B[j*N+i];
		}
	}

	avgA = (double) totA / (N * N);
	avgB = (double) totB / (N * N);
	double a = maxA - minA;
	double b = maxB - minB;
	factor = ((a*a)/avgA) * ((b*b)/avgB);

	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			C[i*N+j] = C[i*N+j] * factor;
		}
	}

	// ********************************************
	// *************** ETAPA 3 ********************
	// ********************************************

	int l;
	double *col, *col_res;
	int *pos, *pos_res;
	col=(double*)malloc(sizeof(double)*N);
	pos=(int*)malloc(sizeof(int)*N);
	col_res=(double*)malloc(sizeof(double)*N);
	pos_res=(int*)malloc(sizeof(int)*N);
	double temp;int pos_temp;


	///////////////////////////////////////////////////////
	// Utilizando vectores para disminuir fallo en cache //
	///////////////////////////////////////////////////////


	// Itera por cada una de las columnas
	for(i=0;i<N;i++){

		// Transforma la columna en un vector
		for (l=0;l<N;l++){
			col[l]=C[i+l*N];
			pos[l]=l;
		}

		// Merge sort al vector
		mergeSort(col,pos,col_res,pos_res,0,N-1);


		////////////////////////////////////////////////////////
		// Ordena la matriz a partir del vector de posiciones //
		////////////////////////////////////////////////////////

		//Utiliza el vector pos[] para ordenar las filas de la matriz hacia la derecha
		for (k=0;k<N;k++){
			for (l=i;l<N;l++){
				CO[l+k*N] = C[l+(pos[k]*N)];
			}
		} 
	}

	printf("Tiempo en segundos total: %f\n\n", dwalltime() - timetick);  

free(A);
free(B);
free(C);
free(D);
free(CO);
return(0);
}