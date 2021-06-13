//
//  main.c
//  proyecto_fin_carrera
//
//  Created by Raúl Martínez Moreno on 10/07/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.


#include <stdio.h>
#include "rngs.h"
#include "genetico.h"
#include <time.h>



int main(int argc, const char * argv[])
{
    int empresa=0,i,j,k,numero=0;
    int *empresas;
    double *inputs,
        *outputs,
        resultado=0.0;

    double t1,t2,t3,t4;

/*---------------------------- PARAMETROS METAHEURISTICAS -------------------------------*/

    int bucle=0;


/*------------------------------------------------------------------------*/

/* ------------------- Parametros Hiperheuristica ----------------------- */

    Ptr_hiperheuristica hiper = NULL;
    hiper = (Ptr_hiperheuristica)malloc(sizeof(struct hiperheuristica));

/*------------------------------------------------------------------------*/

    //Con PlantSeeds -1 el estado se recoge del reloj
    PlantSeeds(-1);
    srand(time(NULL));

    FILE *fsolucion;

    t3=time(0);

    inputs = (double*)malloc((CANTIDAD_M*CANTIDAD_N+1)*sizeof(double));
    outputs = (double*)malloc((CANTIDAD_M*CANTIDAD_N+1)*sizeof(double));
    empresas = (int*)malloc(CANTIDAD_N*sizeof(int));


    //printf("\n-----------------Entradas--------------");
    inputs = generar_vector_double(CANTIDAD_M, CANTIDAD_S, CANTIDAD_N, 1);
    //printf("\n-----------------Salidas--------------");
    outputs = generar_vector_double(CANTIDAD_M, CANTIDAD_S, CANTIDAD_N, 2);

    fsolucion=fopen("solucion.txt","a+");


    /*AQUI LANZAMOS EL PROGRAMA DE EFICIENCIA Y RECOGEMOS LAS DMU A CALCULAR*/
    //empresas = (int*)calloc((CANTIDAD_N),sizeof(int));

    numero = eficiencia(empresas,inputs,outputs);

    i=0;

    for(i=0;i<numero;i++){
            fprintf(fsolucion,"\nEmpresa[%d] = %d",i,empresas[i]);
            printf("\nEmpresa[%d] = %d",i,empresas[i]);
    }
    fprintf(fsolucion,"\nNUMERO = %d",numero);
    printf("\nNUMERO = %d",numero);

    fclose(fsolucion);

    //parametros = parametros_ini+parametros_sel+parametros_com+parametros_imp+parametros_end;

    hiper->NIM = 1;
    hiper->NFM = 0;
    hiper->NMS = 0;

    /* CREAR VECTOR DONDE ALMACENAR TODOS LOS ESQUEMAS PARAMETRIZADOS */

    Ptr_metaheuristica *List_meta = NULL;
    List_meta=(Ptr_metaheuristica*)malloc((hiper->NIM+hiper->NFM)*sizeof(Ptr_metaheuristica));

    /*------------------------------------------------------------------*/

    while(bucle<hiper->NIM){

    Ptr_metaheuristica meta = NULL;
    meta = (Ptr_metaheuristica)malloc(sizeof(struct metaheuristica));

    fsolucion=fopen("solucion.txt","a+");
    //INICIALIZACION DE CROMOSOMAS
/*
    meta->INEIni = (rand()%200)+1;  //INEIni
    meta->FNEIni = (rand()%100)+1;
    if(meta->FNEIni>meta->INEIni){meta->FNEIni=meta->INEIni;}
    meta->CANTIDAD_CRUCE = 2;
    meta->NBESel = (rand()%100)+1;
    meta->NWESel = (rand()%10)+1;
    if(meta->NBESel > meta->FNEIni){meta->NBESel = meta->FNEIni;}
    if(meta->NWESel > meta->FNEIni){meta->NWESel = meta->FNEIni;}
    meta->PBBCom = (rand()%90)+51;
    meta->PBWCom = 0;
    meta->PWWCom = (rand()%90)+1;
    meta->MNIEnd = (rand()%100)+1;
    meta->PEIIni = (rand()%100)+1;
    meta->IIEIni = (rand()%50)+1;
    meta->NIREnd = (rand()%25)+1;
    meta->PEIImp = (rand()%100)+1;
    meta->IIEImp = (rand()%5)+1;
    meta->PEDImp = (rand()%10)+1;
    meta->IIDImp = (rand()%5)+1;
*/

    meta->INEIni = 185;  //INEIni
    meta->FNEIni = 67;
    if(meta->FNEIni>meta->INEIni){meta->FNEIni=meta->INEIni;}
    meta->CANTIDAD_CRUCE = 2;
    meta->NBESel = 72;
    meta->NWESel = 67;
    if(meta->NBESel > meta->FNEIni){meta->NBESel = meta->FNEIni;}
    if(meta->NWESel > meta->FNEIni){meta->NWESel = meta->FNEIni;}
    meta->PBBCom = 84;
    meta->PBWCom = 0;
    meta->PWWCom = 81;
    meta->MNIEnd = 91;
    meta->PEIIni = 92;
    meta->IIEIni = 18;
    meta->NIREnd = 46;
    meta->PEIImp = 60;
    meta->IIEImp = 10;
    meta->PEDImp = 59;
    meta->IIDImp = 9;


    printf("\n-------------------------------------------------------------------------------");
    printf("\nINEIni = %d || FNEIni = %d || PEIIni = %d || IIEIni = %d", meta->INEIni, meta->FNEIni, meta->PEIIni, meta->IIEIni);
    printf("\nMNIEnd = %d || NIREnd = %d", meta->MNIEnd, meta->NIREnd);
    printf("\nNBESel = %d || NWESel = %d", meta->NBESel, meta->NWESel);
    printf("\nPBBCom = %d || PBWCom = %d || PWWCom = %d", meta->PBBCom, meta->PBWCom, meta->PWWCom);
    printf("\nPEIImp = %d || IIEImp = %d || PEDImp = %d || IIDImp = %d\n", meta->PEIImp, meta->IIEImp, meta->PEDImp, meta->IIDImp);
    printf("-------------------------------------------------------------------------------\n");
    //fprintf(fsolucion,"--------------------------------------------------------------------------------");
    //fprintf(fsolucion,"\nINEIni = %d || FNEIni = %d || PEIIni = %d || IIEIni = %d", meta->INEIni, meta->FNEIni, meta->PEIIni, meta->IIEIni);
    //fprintf(fsolucion,"\nMNIEnd = %d || NIREnd = %d", meta->MNIEnd, meta->NIREnd);
    //fprintf(fsolucion,"\nNBESel = %d || NWESel = %d", meta->NBESel, meta->NWESel);
    //fprintf(fsolucion,"\nPBBCom = %d || PBWCom = %d || PWWCom = %d", meta->PBBCom, meta->PBWCom, meta->PWWCom);
    //fprintf(fsolucion,"\nPEIImp = %d || IIEImp = %d || PEDImp = %d || IIDImp = %d\n", meta->PEIImp, meta->IIEImp, meta->PEDImp, meta->IIDImp);
    //fprintf(fsolucion,"-------------------------------------------------------------------------------\n");
    //fprintf(fsolucion,"TAMAÑO DEL PROBLEMA --> m = %d  s = %d  n = %d\n",CANTIDAD_M,CANTIDAD_S,CANTIDAD_N);
    printf("TAMANO DEL PROBLEMA --> m = %d  s = %d  n = %d\n",CANTIDAD_M,CANTIDAD_S,CANTIDAD_N);
    empresa=0;
    fflush(stdout);

    while(empresa<numero){
    printf ("\n\t\t   ========== EMPRESA %d ==========\n", empresas[empresa]);
    resultado += genetico(inputs, outputs, empresas[empresa]-1, fsolucion, meta);
    empresa++;
    }
    meta->fitness = resultado;
    List_meta[bucle]= meta;
    fprintf(fsolucion, "\n--------------------- RESULTADO FINAL = %f-------------------------\n",resultado);
    resultado=0;
    bucle++;
    t4=time(0);
    printf ("\n\n\n  El tiempo (reloj de pared) de ejecucion del programa ha sido (%f) segundos  \n\n\n", t4-t3);
    //fprintf (fsolucion,"\n\n\n  El tiempo (reloj de pared) de ejecucion del programa ha sido (%f) segundos  \n\n\n", t4-t3);

    fflush(stdout);
    }

   /*for(i=0;i<hiper->NIM;i++){
    printf("\nFitness[%d] = %f || %d %d %d",i,List_meta[i]->fitness, List_meta[i]->INEIni, List_meta[i]->FNEIni, List_meta[i]->PBBCom);
   }*/

    generar_meta(List_meta,hiper);

    for(i=0;i<hiper->NFM;i++){

    empresa=0;
    fflush(stdout);
    while(empresa<numero){
    printf ("\n\t\t   ========== EMPRESA %d ==========\n", empresas[empresa]);
    resultado += genetico(inputs, outputs, empresas[empresa]-1, fsolucion, List_meta[bucle]);
    empresa++;
    }
    List_meta[bucle]->fitness = resultado;
    fprintf(fsolucion, "\n--------------------- RESULTADO FINAL = %f-------------------------\n",resultado);
    resultado=0;
    bucle++;
    resultado=0;
    fflush(stdout);
    }


    fflush(stdout);
    for(i=0;i<(hiper->NIM+hiper->NFM);i++){
        fprintf(fsolucion,"\nFitness_ordenado[%d] = %f || %d %d %d %d %d %d %d %d %d %d %d %d %d %d",i,List_meta[i]->fitness, List_meta[i]->INEIni,List_meta[i]->FNEIni,List_meta[i]->PEIIni,List_meta[i]->IIEIni,List_meta[i]->MNIEnd,List_meta[i]->NIREnd,List_meta[i]->NBESel,List_meta[i]->NWESel,List_meta[i]->PBBCom,List_meta[i]->PWWCom,List_meta[i]->PEIImp,List_meta[i]->IIEImp,List_meta[i]->PEDImp,List_meta[i]->IIDImp);
    }

    fclose(fsolucion);
    free(empresas);
    free(inputs);
    free(outputs);

    return 0;
}


