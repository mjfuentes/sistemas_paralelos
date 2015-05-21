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

int main(int argc,char*argv[]){
	double *A,*B,*C,*AB,*D;
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

   //Inicializa las matrices A y B en 1, C en diagonal
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			A[i*N+j] = rand()%10;
			B[j*N+i] = rand()%10;
		}
		D[i] = rand()%10;
	}

	// imprimeMatriz(A,N,1);
	// imprimeMatriz(B,N,0);
	// imprimeMatriz(D,N,1);

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

	// imprimeMatriz(C,N,1);
	printf("Tiempo en segundos de multiplicacion de matrices: %f\n\n", dwalltime() - timetick);  

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

	timetick = dwalltime();

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

	printf("Tiempo en segundos de multiplicacion por factor: %f\n\n", dwalltime() - timetick);  

	printf("Antes de la ordenacion: \n");
  	imprimeMatriz(C,N,1);

	// ********************************************
	// *************** ETAPA 3 ********************
	// ********************************************

	int l;double *col;int *pos;
	col=(double*)malloc(sizeof(double)*N);
	pos=(int*)malloc(sizeof(int)*N);
	double temp;int pos_temp;

	// for(i=0;i<N;i++){
	// 	for (l=0;l<N-1;l++){
	// 		for(j=0;j<N-l;j++){
	// 			if (C[i+j*N] < C[i+(j+1)*N]){
	// 				for(k=i;k<N;k++){
	// 					temp = C[k+j*N];
	// 					C[k+j*N] = C[k+(j+1)*N];
	// 					C[k+(j+1)*N] = temp;
	// 				}
	// 			}
	// 		}
	// 	}
	// }

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

		for (l=0;l<N;l++){
			C[i+l*N] = col[l];
		}



		// Utiliza el vector pos[] para ordenar las filas de la matriz hacia la derecha
		// for (k=0;k<N;k++){
		// 	if (k != pos[k]){
		// 		for (l=i;l<N;l++){
		// 			temp = C[l+k*N];
		// 			C[l+k*N] = C[l+(pos[k]*N)];
		// 			C[l+(pos[k]*N)] = temp;
		// 			pos[k]=pos[k];
		// 		}
		// 	}
		// }
			// printf("pos: %i\n", pos[k]);  
	}

	printf("Despues de la ordenacion: \n");
  	imprimeMatriz(C,N,1);




	// FRANCO  -- > burbuja para el secuencial
	// ENZO    -- > merge sort, cada cachito un hilo. MEJOR OPCION, RECURSIVO. MAS FACIL PARA HACERLO PARALELO



free(A);
free(B);
free(C);
free(D);
return(0);
}