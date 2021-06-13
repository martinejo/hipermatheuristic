#include "generador.h"

/** @brief Genera una matriz de numeros aleatorios.
 *  @param n numero de Columnas.
 *  @param m numero de filas.
 *  @param max valor maximo a generar.
 *
 */
int** generar_matriz (unsigned int n, unsigned int m, unsigned int max)
{
	clock_t ticks;
	int **matriz;
	int i = 0, j = 0;

	ticks = clock();

	srand(ticks);

	matriz = (int**)malloc(m*sizeof(int*));

	for ( i = 0 ; i < m ; i++ )
	{
		*(matriz+i) = (int*)malloc(n*sizeof(int));
	}

	//m filas y n columnas
	for ( i = 0 ; i < m ; i++ )
	{
		for ( j = 0 ; j < n ; j++ )
			{
			*(*(matriz+i)+j) = (int)rand()%max;//(int)rand()%max;
			}
	}
	return matriz;
}
/**
 * Genera una matriz m x n de valores double hasta max.
 * n Columnas
 * m Filas
 */
double** generar_matriz_double (unsigned int m, unsigned int n, unsigned int max) {
    //clock_t ticks;
    //time_t seconds;

	double **matriz;
	int i = 0, j = 0;




    //Reserva de memoria
	matriz = (double**)malloc(m*sizeof(double*));
	for ( i = 0 ; i < m ; i++ )
	{
		matriz[i] = (double*)malloc(n*sizeof(double));
	}
	//m filas y n columnas
	for ( i = 0 ; i < m ; i++ )
	{
		for ( j = 0 ; j < n ; j++ )
        {
			matriz[i][j] = ((double)(rand()/(double)RAND_MAX)+rand()%max);
        }

	}
	return matriz;

}

double* generar_vector_double (unsigned int m, unsigned int s, unsigned int n, unsigned int max) {
    //clock_t ticks;
    //time_t seconds;

FILE *fentradas, *fsalidas;

	double *matriz, *matriz2;
	int i = 0, j = 0, k=0;

    //Reserva de memoria
	matriz = (double*)malloc((m*n)*sizeof(double));
	matriz2 = (double*)malloc((s*n)*sizeof(double));
	//m filas y n columnas
	if(max==1){
        fentradas=fopen("entradas.txt","r");
        for(i=0;i<m;i++){
            for(j=0;j<n;j++){
                fscanf(fentradas,"%lf",&matriz[(i*n)+j]);
            }
        }
        fclose(fentradas);
        return matriz;
    }
    if(max==2){
        fsalidas=fopen("salidas.txt","r");
        for(i=0;i<s;i++){
            for(j=0;j<n;j++){
                fscanf(fsalidas,"%lf",&matriz2[(i*n)+j]);
            }
        }
        fclose(fsalidas);
        return matriz2;
    }
    return matriz;
}
/** @brief Genera un cromosoma aleatorio de longitud n.
 *  @param numero_genes Longitud del cromosoma.
 *
 * El cromosoma será un vector de enteros de longitud 'numero_genes'.
 * En cada posición habrá aleatoriamente un valor 0 o 1.
 *
 * @return cromosoma
 */

int* generate_random_int_chain(unsigned const int length) {

	int i = 0;
	int* cromosoma;

	cromosoma = (int*)malloc(length*sizeof(int));
	srand (time (NULL));

	for ( i = 0 ; i < length ; i++)
    {
		cromosoma[i] = (int)rand()%2;
    }
	return cromosoma;
}

/**
 * Genera una cadena de doubles.
 */
double* generate_random_double_chain (unsigned int length) {
	int i = 0;
	double* cromosoma;

	cromosoma = (double*)malloc(length*sizeof(double));

	for ( i = 0 ; i < length ; i++)
    {
		*(cromosoma+i) = (double)Random();
    }

	return cromosoma;
}

void imprimir_matriz (double **matriz, unsigned int m, unsigned int n)
{
	int i = 0, j = 0;

		for ( i = 0 ; i < m ; i++ )
		{
			for ( j = 0 ; j < n ; j++ ) {
				printf("%lf\t",matriz[i][j]);
            }
			printf("\n");
		}

}

/** @brief Resume la matriz quedandose con los elementos que son 1s en el cromosoma dado
 *  @param matriz
 *  @param cromosoma
 *  @param longitud Numero de columnas de la matriz y longitud del cromosoma
 *  @param numero_unos Longitud de la matriz resumida (al terminar la función).
 *
 *  Se da por supuesto que la longitud del cromosoma corresponderá con el numero de columnas de la matriz
 */


int** resumir_matriz (int **matriz, int *cromosoma, unsigned int longitud)
{
	int i = 0, j = 0;
	int** matriz_resumida;
    int numero_unos;

	numero_unos = count_ones_chain(cromosoma,longitud);

	matriz_resumida = (int**)malloc(numero_unos*sizeof(int*));

	for ( i = 0 ; i < longitud ; i++ )
	{
		if ( *(cromosoma+i) == 1 ) // Si vale 1 hay que copiar el puntero a la nueva matriz
		{
			*(matriz_resumida+j) = *(matriz+i);
			j++;
		}
	}

	return matriz_resumida;
}
/*
int* generar_cromosoma(unsigned const int numero_genes) {
    return generate_random_int_chain(numero_genes);
}
*/
