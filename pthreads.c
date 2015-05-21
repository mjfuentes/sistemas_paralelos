#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <pthread.h>

double *A,*B,*C,*AB,*D;
int N,T;
double res[6*T];
double pow (double base , double exponent);
float powf (float base  , float exponent);
long double powl (long double base, long double exponent);
double maxA = (double) INT_MIN;
double maxB = (double) INT_MIN;
double minB = (double) INT_MAX;
double minA = (double) INT_MAX;
double totA = 0;
double totB = 0;
double avgA,avgB,totA,totB,factor;
int pthread_barrier_t barrier;

//Para calcular tiempo
double dwalltime(){
	double sec;
	struct timeval tv;
	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}

void *calcular(void *s){
	int principio,final,i,id,j,k;
	int *id_pointer;
	id_pointer=(int*) s;
	id=*id_pointer;
	principio=id*(N/T);
	final=principio+(N/T);

	for (i=principio; i<final; i++){
		for(j=0;j<N;j++){
			for(k=0;k<N;k++){
				C[i*N+j] += A[i*N+k] * B[j*N+k];
			}
			C[i*N+j] = C[i*N+j] * D[j]; 
		}
	}

	for (i=principio; i<final; i++){
		for(j=0;j<N;j++){
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

	int pthread_barrier_wait (pthread_barrier_t *barrier);
	
	if (id == 0){
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

	}

	int pthread_barrier_wait (pthread_barrier_t *barrier);

	for(i=principio;i<final;i++){
		for(j=0;j<N;j++){
			C[i*N+j] = C[i*N+j] * factor;
		}
	}

	int pthread_barrier_wait (pthread_barrier_t *barrier);



}

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

int main(int argc,char*argv[]){
	int id[T];
	int i,j,k;
	int check=1;
	double timetick;
	T = 4;
	pthread_t p_threads[T];
	pthread_attr_t attr;
	pthread_attr_init (&attr);

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

   //Inicializa las matrices A y B en 1, C en diagonal
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			A[i*N+j] = rand()%10;
			B[j*N+i] = rand()%10;
			C[i*N+j] = 0;
		}
		D[i] = rand()%10;
	}

	pthread_barrier_init (&barrier, NULL, 4);

	// ********************************************
	// *************** ETAPA 1 ********************
	// ********************************************

	timetick = dwalltime();

	// for(i=0;i<N;i++){
	// 	for(j=0;j<N;j++){
	// 		for(k=0;k<N;k++){
	// 			C[i*N+j] += A[i*N+k] * B[j*N+k];
	// 		}
	// 		C[i*N+j] = C[i*N+j] * D[j]; 
	// 	}
	// }   

	for(i=0;i<T;i++){ 
		id[i] = i;
		pthread_create(&p_threads[i], &attr, calcular, (void *) &id[i]);
	}

	for (i=0; i< T; i++)
	{ 
		pthread_join(p_threads[i], NULL);
	}

	printf("Tiempo en segundos total: %f\n\n", dwalltime() - timetick);  

	// ********************************************
	// *************** ETAPA 3 ********************
	// ********************************************

	int l;double *col;int *pos;
	col=(double*)malloc(sizeof(double)*N);
	pos=(int*)malloc(sizeof(int)*N);
	double temp;int pos_temp;
	timetick = dwalltime();

	// Itera por cada una de las columnas
	for(i=0;i<N;i++){

		// Transforma la columna en un vector
		for (l=0;l<N;l++){
			col[l]=C[i+l*N];
			pos[l]=l;
		}

		// Utiliza bubble sort para ordenar el vector
		// Guarda en Col[] los valores del vector ordenados, y en pos[] los indices iniciales del vector ordenados
		for (l=0;l<N;l++){
			for(j=0;j<N-l;j++){
				if (col[j] < col[j+1]){
					temp = col[j];
					col[j] = col[j+1];
					col[j+1] = temp;
					pos_temp = pos[j];
					pos[j] = pos[j+1];
					pos[j+1] = pos_temp;
				}
			}
		}

		//Utiliza el vector pos[] para ordenar las filas de la matriz hacia la derecha
		for (k=0;k<N;k++){
			if (k != pos[k]){
				for (l=i;l<N;l++){
					temp = C[l+k*N];
					C[l+k*N] = C[l+(pos[k]*N)];
					C[l+(pos[k]*N)] = temp;
					for (j=k;j<N;j++){
						if (pos[j]==k){
							pos[j]=pos[k];
						}
					}
				}
			}
		}

	}
	printf("Tiempo de la ordenacion: %f\n\n", dwalltime() - timetick);  


// free(A);
// free(B);
// free(C);
// free(D);
// free(pos);
// free(col);

return(0);
}