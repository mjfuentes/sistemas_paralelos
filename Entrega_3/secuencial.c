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

void imprimeVector(double *S,int N,char comment[]){
		int i,j,I,J,despB;
		printf("%s",comment);
		printf("Contenido del vector: \n" );
		for(j=0;j<N;j++){
			printf("%f ",S[j]);
		}
		printf("\n\n");
}

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

void mergeSort(double *col,double *col_res, int inicio, int final)
{
	if (inicio == final) return;
	int mid = (inicio + final)/2;
	mergeSort(col,col_res, inicio, mid);
	mergeSort(col,col_res,mid+1, final);
	merge(col,col_res, inicio, mid, final);
}

int main(int argc,char*argv[]){
	double *D;
	int i,j,k,N;
	int check=1;
	double timetick;

	if (argc < 2){
		printf("\n Falta un argumento:: N dimension del vector \n");
		return 0;
	}
	N=atoi(argv[1]);

   //Aloca memoria para el vector
	D=(double*)malloc(sizeof(double)*N);

   //Inicializa el vector con numeros random
	srand(time(0));
	for(i=0;i<N;i++){
		D[i] = rand()%100+1;
	}

	int l;
	double *col_res;
	col_res=(double*)malloc(sizeof(double)*N);

	// Merge sort al vector
	mergeSort(D,col_res,0,N-1);
	
	imprimeVector(D,N,"vector ordenado");
	printf("Tiempo en segundos total: %f\n\n", dwalltime() - timetick);  

free(D);
free(col_res);
return(0);
}