

#ifndef proyecto_fin_carrera_generador_h
#define proyecto_fin_carrera_generador_h

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include "chains_tools.h"
#include "rngs.h"


int** generar_matriz (unsigned int, unsigned int, unsigned int);
double** generar_matriz_double (unsigned int, unsigned int, unsigned int);
double* generar_vector_double(unsigned int, unsigned int, unsigned int, unsigned int);
int* generate_random_int_chain(unsigned const int);
double* generate_random_double_chain (unsigned int);

int** resumir_matriz (int**, int*, unsigned int);
void imprimir_matriz (double**, unsigned int, unsigned int);

#endif
