//
//  chains_tools.c
//  proyecto_fin_carrera
//
//  Created by Raúl Martínez Moreno on 13/06/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include "chains_tools.h"

/*
 * Devuelve el numero de valores(0s o 1s) que varian en una cadena con respecto a la otra.
 * Por ejemplo 0101 y 0100 se diferencia en un valor, pero 0000 y 1111 en 4.
 * @param Cadena 1.
 * @param Cadena 2.
 * @param Longitud cadenas.
 */
int count_ones_chain (int* cromosoma, unsigned int longitud)
{
	int i = 0, cont = 0;

	for ( i = 0 ; i < longitud ; i++)
		if (*(cromosoma+i) == 1 )
			cont++;

	return cont;
}
/**
 * Imprime una cadena por pantalla.
 * @param *c cadena.
 * @param l Longitud de la cadena.
 */
void print_chain(int *c, int l) {
    int i = 0;

    for ( i = 0 ; i < l ; i++ )
        printf("%d ",*(c+i));
}
/**
 * Compara dos candenas dadas.
 * @param *c1 Primera cadena.
 * @param *c2 Segunda cadena.
 * @param l Longitud cadenas.
 * @return 0 si son iguales, 1 si no.
 */
int compare_chain (int *c1, int *c2, int l) {
    int i = 0, equal = 0;

    while (equal == 0 && i < l) {
        if ( *(c1+i) != *(c2+i) )
            equal = 1;
        i++;
    }
    return equal;
}
/**
 * @brief Cuenta el numero de valores diferentes entre dos cadenas.
 * @param *c1 Primera cadena.
 * @param *c2 Segunda cadena.
 * @param l Longitud cadenas.
 * @return Numero de elementos diferentes entre dos cadenas.
 */
int number_different_values (int *c1, int *c2, int l) {
    int count = 0, i = 0;

    for ( i = 0 ; i < l ; i++ ) {
        if ( *(c1+i) != *(c2+i) )
            count++;
    }

    return count;
}
/**
 * @brief Devuelve una copia de una cadena dada.
 * @param *c Cadena a copiar.
 * @param l Longitud cadena.
 * @return Cadena copia.
 */
int* copy_chain(int *c, int l) {
    int *aux,
    i = 0;

    aux = (int*)malloc(l*sizeof(int));

    for ( i = 0 ; i < l ; i++ )
        *(aux+i) = *(c+i);

    return aux;
}
/**
 * Devuelve una nueva cadena que es el resultado de realizar la operacion XOR sobre las dos cadenas dadas.
 * @param c1 Cadena 1.
 * @param c2 Cadena 2.
 * @param l longitud cadenas.
 * @return Resultado XOR sobre las dos cadenas dadas.
 */
int* combine_chain_xor(int *c1, int *c2, int l) {
    int *aux,
    i = 0;

    aux = (int*)malloc(l*sizeof(int));

    for ( i = 0 ; i < l ; i++ ) {
        if ( *(c1+i) != *(c2+i) )
            *(aux+i) = 1;
        else
            *(aux+i) = 0;
    }

    return aux;
}

/**
 * Combina dos cadenas dadas. Intercambia la primera mitad de las cadenas dadas.
 * @param *c1 Cadena 1.
 * @param *c2 Cadena 2.
 */
void combine_chain_head(int *c1, int *c2, int l) {
    int *aux,
    i = 0;

    aux = (int*)malloc((l/2)*sizeof(int));
    for ( i  = 0 ; i < l/2 ; i++) {
        *(aux+i) = *(c1+i);
        *(c1+i) = *(c2+i);
        *(c2+i) = *(aux+i);
    }

    free(aux);
}
/**
 * Combina dos cadenas dadas. Intercambia la segunda mitad de las cadenas dadas.
 * @param *c1 Cadena 1.
 * @param *c2 Cadena 2.
 */
void combine_chain_tail(int *c1, int *c2, int l) {
    int *aux,
    i = 0;

    aux = (int*)malloc(((l+1)/2)*sizeof(int));
    for ( i  = l/2 ; i < l ; i++) {
        *(aux+i) = *(c1+i);
        *(c1+i) = *(c2+i);
        *(c2+i) = *(aux+i);
    }

    free(aux);
}
/**
 * @brief cambia un valor aleatorio en una cadena
 * @param c Cadena.
 * @param l Longitud cadena.
 */
void change_random_value(int *c,int l) {
    int posicion = 0;
    clock_t ticks;

    ticks = clock();
    srand(ticks);
    posicion = rand()%l;

    if ( *(c+posicion) == 0 )
        *(c+posicion) = 1;
    else {
        *(c+posicion) = 0;
    }
}
/**
 * @brief Cambia el valor de una posicion dada en una cadena.
 * @param c Cadena.
 * @param p Posicion.
 * @param l Longitud cadena.
 */
void change_value(int *c, int p, int l) {
    if ( p < l ){
        if ( *(c+p) == 0 )
            *(c+p) = 1;
        else {
            *(c+p) = 0;
    }
    }
}
/**
 * @brief Intercambia los valores de dos cadenas en intervalos de un valor dado.
 * Ejemplo: Para una cadena de longitud 10, si se da el intervalo 2, se intercambiaran los valores 0,1,4,5,8,9.
 * Para una cadena de longitud 10, si se da el intervalo 5, se intercambiaran los valores 0,1,2,3,4.
 * @param *c1 Primera cadena.
 * @param *c2 Segunda cadena.
 * @param lenght Longitud de las cadenas.
 * @param intervalo Intervalo
 */
void exchange_chains_values (int *c1, int *c2, int lenght, int interval) {
    int i = 0,
        aux = 0,
        aux_int = 0,
        exchange = 0;

    if ( interval < lenght ) {
        for ( i = 0 ; i < lenght ; i++ ) {
            if ( aux == 0 ) {
                aux = interval;
                if ( exchange == 0 )
                    exchange = 1;
                else
                    exchange = 0;
            }
            if ( exchange == 1 ) {
                aux_int = *(c1+i);
                *(c1+i) = *(c2+i);
                *(c2+i) = aux_int;
            }
            aux--;
        }
    }
}
double difference (double *c, int starts, int ends, double difference)
{
    double aux = 0.0;
    int i = 0;

    for ( i = starts ; i < ends ; i++ ) {
        aux+=c[i];  //Acumula en aux el sumatorio
    }
    aux = difference-aux;   //Diferencia menos sumatorio

    return aux;
}

void change_chain (double *chain, int starts, int ends)
{
    double diff, cant;
    int i = 0, num_no_ceros = 0, minor, mayor;

    for ( i = starts ; i < ends ; i++ )
        if ( chain[i] != 0.0 )num_no_ceros++;
        diff = difference(chain, starts, ends, 1.0);  // Cantidad que hay que sumar (si es positiva) o restar (si es negativa)

        if ( num_no_ceros == 0 )
            diff = 0.0;
        else
            cant = diff/num_no_ceros;

            i = starts;
            while (diff < -0.2 || diff > 0.2 ) {
                if ( chain[i] != 0.0 && chain[i]+cant > 0.0 ) {
                chain[i]+=cant;
                diff-=cant;
            }
            i++;
            if ( i == ends ) {
                i = starts;
                cant = diff/num_no_ceros;
            }
        }

    // Busca el primer elemento.
    i = starts;
    while ( i < ends && chain[i] == 0.0 )
        i++;

    if ( i < ends ) {
        minor = i;
        mayor = i;
    }

    for ( i = starts ; i < ends ; i++ ) {
        if ( chain[i] < chain[minor] && chain[i] != 0.0 )
            minor = i;
        if ( chain[i] > chain[mayor] )
            mayor = i;
    }

    if ( minor != mayor ) {
        if ( diff < 0.0 )
            chain[mayor]+=diff;
        else if ( diff > 0.0 )
            chain[minor]+=diff;
    }
}

