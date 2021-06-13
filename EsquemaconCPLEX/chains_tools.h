//
//  chains_tools.h
//  proyecto_fin_carrera
//
//  Created by Raúl Martínez Moreno on 13/06/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef proyecto_fin_carrera_chain_tools_h
#define proyecto_fin_carrera_chain_tools_h

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "rngs.h"

int* copy_chain(int*, int);

int count_ones_chain (int*, unsigned int);
int number_different_values (int*, int*, int);
int compare_chain (int*,int*, int);

void print_chain(int*,int);

// Cambia un valor aleatorio en una cadena.
void change_random_value(int*,int);
// Cambia un dado concreto en una cadena.
void change_value(int*,int,int);
// Combinaciones de cadenas.
int* combine_chain_xor(int*, int*, int);
void combine_chain_head(int*, int*, int);
void combine_chain_tail(int*, int*, int);

double difference (double *, int, int, double);
void change_chain (double *chain, int starts, int ends);

#endif
