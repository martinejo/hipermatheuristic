//#include <iostream>
//#include <windows.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/time.h>
#include <time.h>
//#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <string.h>
//#include <omp.h>
#include "generador.h"
#include "chains_tools.h"


/*
struct Chromosome {
	int etiq;
	int evaluated;		//0 no evaluado, 1 si evaludado
	double evaluation;	//tiene la puntuacion del cromosoma, al crearlo es BAD_CHROM
	int *model;	//es el modelo, la tarea y el procesador
	int longitud;
};*/

struct Chromosome {
    int etiq;                   // Etiqueta del cromosoma.

    double  evaluation,         // Valor del cromosoma
            pre_evaluation,
            beta,
            *vik,
            *urk,
            *tneg,
            *tpos,
            *djk,
            *alfas;
    int     *discrete_model,
            metodo,
            fi,
            valido,
            contador,
            terminar;
};

struct hiperheuristica{

    int NIM, //Cantidad de metaheuristicas generadas inicialmente
        NFM, //Cantidad de metaheuristicas generadas por combinacion
        NMS; //Numero de metaheuristicas seleccionadas
};

struct metaheuristica{

     int INEIni,
        FNEIni,
        CANTIDAD_CRUCE,
        NBESel,
        NWESel,
        PBBCom,
        PBWCom,
        PWWCom,
        MNIEnd,
        PEIIni,
        IIEIni,
        NIREnd,
        PEIImp,
        IIEImp,
        PEDImp,
        IIDImp,
        parametros_ini,
        parametros_sel,
        parametros_com,
        parametros_imp,
        parametros_end,
        parametros;

    double fitness;
    };

typedef struct Chromosome * Ptr_Chromosome;     //declaramos un ptr a chromosome
typedef struct hiperheuristica * Ptr_hiperheuristica;
typedef struct metaheuristica * Ptr_metaheuristica;

/*Constantes*/
#define BAD_CHROM  999999
#define MAX_REP_MEJOR 10      // Maximo de veces que se repite el mejor
#define MAX_ITER 100000         // Maximo de iteraciones del mejor

#define NUMERO_FILAS 15
#define NUMERO_COLUMNAS 10
#define VALOR_MAX 99
#define BASADA_PUNTO 1
#define MAX_CHROM 100

#define CANTIDAD_M 2               //ENTRADAS
#define CANTIDAD_N 30
#define CANTIDAD_S 1             //SALIDAS
#define POSITION_V CANTIDAD_M
#define POSITION_T_POS CANTIDAD_M*2
#define POSITION_NUN POSITION_T_POS+CANTIDAD_S
#define POSITION_ALFA POSITION_NUN+CANTIDAD_S
#define POSITION_D POSITION_ALFA+CANTIDAD_N
#define LONG_CHROMOSOME_CONTINOUS 2*(CANTIDAD_M+CANTIDAD_N+CANTIDAD_S)

/*--------------------------- FUNCIONES ----------------------------------*/

Ptr_Chromosome create_chromosome(int, double *, double *, int, double*, double*, double*, double*, double*, double*, int);
void generate_b(Ptr_Chromosome, int, double*, double*);
void generate_alfas(Ptr_Chromosome, double *, double*);
void generate_beta (Ptr_Chromosome, double*,int);
int modify_alfas(Ptr_Chromosome, double *, double *, int, int, double*, double*, double*);
void modify_alfas2(Ptr_Chromosome, double *, double *, int);

double pre_evaluate_chromosome(Ptr_Chromosome, double*, double*,int,double*,double*);
double difference (double *, int,  int , double );
void free_chromosome(Ptr_Chromosome);
void show_chromosome(Ptr_Chromosome, double*, double*, int, double*, double*);
int good_chromosome(Ptr_Chromosome, double *, double *, int, double*, double*);
double evaluate_chromosome(Ptr_Chromosome, double *, double*, int, double*, double*);
void classify_chromosome(Ptr_Chromosome*, int);
void pre_classify_chromosome(Ptr_Chromosome*, int);

void crossover(Ptr_Chromosome*,Ptr_Chromosome, double*, double*, int, double*, double*,int);
void mutate (Ptr_Chromosome, double*, double*, int, double*, double*);
double genetico(double *, double *, int, FILE *, Ptr_metaheuristica);

// Restricciones
double check_constraint_1 (Ptr_Chromosome, double*, int, double, double* );
double check_constraint_2 (Ptr_Chromosome, double*, int, double );
double check_constraint_3 (Ptr_Chromosome, double*, int, double );
double check_constraint_4 (Ptr_Chromosome, double *, double*, int, double);
int check_constraint_5 (Ptr_Chromosome);
int check_constraint_6 (Ptr_Chromosome);
int check_constraint_11 (Ptr_Chromosome);
int check_constraint_11 (Ptr_Chromosome);
int check_constraint_12 (Ptr_Chromosome);
int check_constraint_13 (Ptr_Chromosome);
int check_constraint_14 (Ptr_Chromosome);
int check_constraint_15 (Ptr_Chromosome);

void modify_constraint_1(Ptr_Chromosome,double*,int);
void modify_constraint_2(Ptr_Chromosome,double*,int);
void modify_beta(Ptr_Chromosome, double*, double*, int, int, double*, double*,double*);
void modify_beta_cruce(Ptr_Chromosome, double*, double*, int);
void mejora(Ptr_Chromosome,double*,double*,int,double*,double*, int);

void deducir_valores_alfa(Ptr_Chromosome, double*, int, int, double*, double*, double*);
void deducir_valores_uv(Ptr_Chromosome, double*, double*, int);
void calculo_D(Ptr_Chromosome, double*, double*);
void calculo_alfas (Ptr_Chromosome, double*, int);

void sumatorios(Ptr_Chromosome, double*, double*, double*, double*);
int mejor_caso(Ptr_Chromosome, double*, double*);
void suma_filas(Ptr_Chromosome, double*, double*);
void suma_columnas(double*, double*, double*, double*);
int termino_inout(Ptr_Chromosome, double*, double*, int i);
void inversos(double*,double*,double*,double*);
void recalcular(Ptr_Chromosome,double*,double*,int,double*,double*,int);

double maximizar(Ptr_Chromosome,double*,int);
void classify_chromosome_maximo(Ptr_Chromosome*, double*,int,int);
void mejorar_u(Ptr_Chromosome, double*, double*);
void suma_columnas2(double*, double*, double*, double*);
void generate_beta2 (Ptr_Chromosome, double*, int);
void deducir_valores_uv2(Ptr_Chromosome, double *, double *, int);
Ptr_Chromosome cromosoma_vacio(int);

int eficiencia(int*, double*, double*);
void generar_meta(Ptr_metaheuristica*, Ptr_hiperheuristica);

