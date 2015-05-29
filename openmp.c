#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include <omp.h>

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

int main(int argc,char*argv[]){	double *A,*B,*C,*AB,*D, *CO;
	int i,j,k,N;
	int check=1;
	double timetick_total;
	double timetick;
	omp_set_num_threads(4);
	double etapa1,etapa2; 

	
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

    // Variables para etapa 2	
	int tid;
	double totA,totB,factor, avgA, avgB,maxA,minA,maxB,minB;
	maxA = INT_MIN;
	maxB = INT_MIN;
	minB = INT_MAX;
	minA = INT_MAX;
	totA = 0;
	totB = 0;
	
	
  // Varibles para la etapa 3
	int l;
	double *col, *col_res;
	int *pos, *pos_res;
	col=(double*)malloc(sizeof(double)*N);
	pos=(int*)malloc(sizeof(int)*N);
	col_res=(double*)malloc(sizeof(double)*N);
	pos_res=(int*)malloc(sizeof(int)*N);
	double temp;int pos_temp;
	
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
	timetick_total = dwalltime();
	timetick = dwalltime();

	// i por defecto va al private porq es el indice del for
        #pragma omp parallel default(none) private(tid,i,j,k,l) shared(CO,factor,A,B,C,D,N,totA,totB,avgA,avgB,maxA,maxB,minA,minB,timetick,col,pos,col_res,pos_res,etapa1,etapa2)	
	{ // abro directiva de region 
	    #pragma omp for private(j,k)
		for(i=0;i<N;i++){
			for(j=0;j<N;j++){
				C[i*N+j] = 0;
				for(k=0;k<N;k++){
					C[i*N+j] += A[i*N+k] * B[j*N+k];
				}
				C[i*N+j] = C[i*N+j] * D[j]; 
			}
		}   

		tid = omp_get_thread_num();


#pragma omp single

		{ 

			printf("Matriz A:\n");
			imprimeMatriz(A,N,1);

			printf("Matriz B: \n");
			imprimeMatriz(B,N,1);

			printf("Matriz C: \n");
			imprimeMatriz(C,N,1);

		}



	// ********************************************
	// *************** ETAPA 2 ********************
	// ********************************************


		double localMaxA,localMaxB,localMinA,localMinB;

	#pragma omp for private(j) reduction(+: totA,totB)
		for(i=0;i<N;i++){
			for(j=0;j<N;j++){

				if (A[i*N+j] > localMaxA){
					localMaxA = A[i*N+j];
				} 
				if (A[i*N+j] < localMinA){
					localMinA = A[i*N+j];
				}
				if (B[j*N+i] > localMaxB){
					localMaxB = B[j*N+i];
				} 
				if (B[j*N+i] < localMinB){
					localMinB = B[j*N+i];
				}
				totA += A[i*N+j];
				totB += B[j*N+i];
			}


		}

			#pragma omp critical
		if (localMaxA > maxA){
			maxA = localMaxA;
		}
		#pragma omp critical
		if (localMinA < minA){
			minA = localMinA;
		}
		#pragma omp critical
		if (localMaxB > maxB){
			maxB = localMaxB;
		}
		#pragma omp critical
		if (localMinB < minB){
			minB = localMinB;
		}

		tid = omp_get_thread_num();

	#pragma omp single
		{ 


			avgA = (double) totA / (N * N);
			avgB = (double) totB / (N * N);
			double a = maxA - minA;
			double b = maxB - minB;
			factor = ((a*a)/avgA) * ((b*b)/avgB);

			printf("Factor: %f\n",factor);	
			printf("Max A: %f\n",maxA);	
			printf("Max B: %f\n",maxB);	
			printf("Total A: %f\n",totA);
			printf("Total B: %f\n",totB);

		}

	#pragma omp for private(j) firstprivate (factor)
		for(i=0;i<N;i++){
			for(j=0;j<N;j++){
				C[i*N+j] = C[i*N+j] * factor;
			}
		}

	#pragma omp single
		{
			printf("Matriz C * factor: \n");
			imprimeMatriz(C,N,1);
		}



#pragma omp single
		{
			etapa2 = dwalltime() - timetick;
			timetick = dwalltime(); 


		}


	// ********************************************
	// *************** ETAPA 3 ********************
	// ********************************************


	///////////////////////////////////////////////////////
	// Utilizando vectores para disminuir fallo en cache //
	///////////////////////////////////////////////////////


	// Itera por cada una de las columnas
		for(i=0;i<N;i++){

	      // Transforma la columna en un vector
		#pragma omp for
			for (l=0;l<N;l++){
				col[l]=C[i+l*N];
				pos[l]=l;
			}

		// Merge sort al vector
			mergeSort(col,pos,col_res,pos_res,tid*(N/4),(tid+1)*(N/4)-1);
	//	if (tid == 0) mergeSort(col,pos,col_res,pos_res,0,N-1);

		#pragma omp barrier

			if (tid % 2 == 0){

				int mid = tid*(N/4) + (N/4) - 1;
				merge(col,pos,col_res,pos_res,tid*(N/4),mid, tid*(N/4) + 2*(N/4)-1);
			}

		#pragma omp barrier

			if (tid == 0){
				int mid = (N/2) - 1;
				merge(col,pos,col_res,pos_res,0,mid, N-1);
			}

		#pragma omp barrier

		////////////////////////////////////////////////////////
		// Ordena la matriz a partir del vector de posiciones //
		////////////////////////////////////////////////////////


		#pragma omp for
			for (k=0;k<N;k++){
				for (j=i;j<N;j++){
					CO[k*N+j] = C[pos[k]*N+j];
				}
			}

		#pragma omp single
			{
				double * tmp = C;
				C = CO;
				CO = tmp;
			}

		}

		#pragma omp for
		for (k=0;k<N;k++){
			for (j=0;j<N/2;j++){
				C[k*N+(j*2)] = CO[k*N+(j*2)];
				C[k*N+(j*2)+1] = C[k*N+(j*2)+1];
			}
		}
	}
	// CIERRO EL DEL OMP PARALLEL
	printf("Etapa 1 y 2 terminada en: %f\n", etapa2);
	printf("Etapa 3 terminada en: %f\n",dwalltime() - timetick);
	printf("Tiempo en segundos total: %f\n\n", dwalltime() - timetick_total);  
		//printf("Despues de la ordenacion: \n");
		// printf("Despues del acomodo: \n");
 	// 	imprimeMatriz(CO,N,1);
		// printf("Despues del acomodo: \n");
 	// 	imprimeMatriz(C,N,1);




	free(A);
	free(B);
	free(C);
	free(D);
	free(CO);
	return(0);
}