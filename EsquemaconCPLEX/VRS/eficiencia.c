
#include <stdio.h>
#include <stdlib.h>
#include "genetico.h"
#include <string.h>
#include <math.h>
#include <ilcplex/cplex.h>
#include <time.h>


int eficiencia(int *empresas, double *inputs, double *outputs){


int i,j,k=0,status=0,M=10000,estado=0,eficaz=0,eficientes=1,cantidad=0,contador=0;

double **entrada,
       **salida;

double objval;

double t1,t2;

FILE *fentradas, *fsalidas, *fsolucion, *fmodelo;

CPXENVptr env=NULL;
CPXLPptr lp=NULL;


/* --------------------------------------------------------------------------------------------------*/
/* --------------------------------------------------------------------------------------------------*/
/* --------------------------------------------------------------------------------------------------*/
/* ----------------------------------------- INICIALIZACION------------------------------------------*/
/* --------------------------------------------------------------------------------------------------*/
/* --------------------------------------------------------------------------------------------------*/
/* --------------------------------------------------------------------------------------------------*/


t1=time(0);

fentradas=fopen("entradas.txt","r");
fsalidas=fopen("salidas.txt","r");

entrada = (double **)malloc((CANTIDAD_N+1)*sizeof(double*));
for (i=0;i<CANTIDAD_N+1;i++) {
		entrada[i] = (double*)malloc((CANTIDAD_M+1)*sizeof(double));
}

salida = (double **)malloc((CANTIDAD_N+1)*sizeof(double*));
for (i=0;i<CANTIDAD_N+1;i++) {
		salida[i] = (double*)malloc((CANTIDAD_S+1)*sizeof(double));
}

//empresas = (int*)malloc(CANTIDAD_N*sizeof(int));

//////////////////// GUARDAR DATOS //////////////////////////

for(i=1;i<CANTIDAD_M+1;i++){
    for(j=1;j<CANTIDAD_N+1;j++){
        fscanf(fentradas,"%lf",&entrada[j][i]);
    }
}
for(i=1;i<CANTIDAD_S+1;i++){
    for(j=1;j<CANTIDAD_N+1;j++){
        fscanf(fsalidas,"%lf",&salida[j][i]);
    }
}
fclose(fentradas);
fclose(fsalidas);
fflush(stdout);


/* --------------------------------------------------------------------------------------------------*/
/* --------------------------------------------------------------------------------------------------*/
/* --------------------------------------------------------------------------------------------------*/
/* ----------------------------------------- CALCULO DE EFICIENCIA ----------------------------------*/
/* --------------------------------------------------------------------------------------------------*/
/* --------------------------------------------------------------------------------------------------*/
/* --------------------------------------------------------------------------------------------------*/

fsolucion=fopen("solucion.txt","a+");
fprintf(fsolucion,"\n\n ============ EFICIENTES ============= \n\n");
fclose(fsolucion);

for(k=1;k<CANTIDAD_N+1;k++){

fmodelo=fopen("modelo.lp","w");

fflush(stdout);
/* CREAR MODELO A OPTIMIZAR */

fprintf(fmodelo,"maximize\nfuncion: ");
fflush(stdout);

fprintf(fmodelo,"sn1");
for(i=2;i<CANTIDAD_M+1;i++){
    fprintf(fmodelo,"+sn%d",i);
}
fflush(stdout);

for(i=1;i<CANTIDAD_S+1;i++){
    fprintf(fmodelo,"+sp%d",i);
}
fflush(stdout);
/* CREAR RESTRICCIONES UTILIZADAS */

fprintf(fmodelo,"\nst\n");

for(i=1;i<CANTIDAD_M+1;i++){
        fprintf(fmodelo,"sumatorio_entradas%d: ",i);
        for(j=1;j<CANTIDAD_N+1;j++){
            fprintf(fmodelo,"+%fd%d",entrada[j][i],j);
        }
       fprintf(fmodelo,"+sn%d=%f\n",i,entrada[k][i]);
}
fflush(stdout);

for(i=1;i<CANTIDAD_S+1;i++){
        fprintf(fmodelo,"sumatorio_salidas%d: ",i);
        for(j=1;j<CANTIDAD_N+1;j++){
            fprintf(fmodelo,"+%fd%d",salida[j][i],j);
        }
       fprintf(fmodelo,"-sp%d=%f\n",i,salida[k][i]);
}
fflush(stdout);

fprintf(fmodelo,"VRS: ");
//#pragma omp parallel for schedule (static) private (j)
for(i=1;i<CANTIDAD_N+1;i++){
    fprintf(fmodelo,"+d%d",i);
}
fprintf(fmodelo," = 1\n");
/* CREAR COTAS DE VARIABLES */

fprintf(fmodelo,"bounds\n");

	for(i=1;i<CANTIDAD_N+1;i++){
        fprintf(fmodelo,"d%d>=0\n",i);
	}
	for(i=1;i<CANTIDAD_M+1;i++){
        fprintf(fmodelo,"sn%d>=0\n",i);
	}
	for(i=1;i<CANTIDAD_S+1;i++){
        fprintf(fmodelo,"sp%d>=0\n",i);
	}
fprintf(fmodelo,"end");
fclose(fmodelo);
fflush(stdout);

fsolucion=fopen("solucion.txt","a+");

/*Me conecto a CPLEX*/

env=CPXopenCPLEX(&status);

/*status=CPXsetintparam(env,CPX_PARAM_SCRIND,CPX_ON);*/

status=CPXsetdblparam(env,CPX_PARAM_EPRHS,1e-09);

/*Creo el problema lp*/

lp=CPXcreateprob(env,&status,"problema");

/*Leo el modelo lp desde el fichero*/

status=CPXreadcopyprob(env,lp,"modelo.lp",NULL);

/*Resuelvo el problema de programación -lineal-*/

status=CPXprimopt(env,lp);

/*accedo al valor óptimo*/

status=CPXgetobjval(env,lp,&objval);

/*Escribo el valor óptimo del problema */
if(fabs(objval) < 0.00000001)
    {
        eficaz=1;
        fprintf(fsolucion,"DMU%d= %.10lf\t|||| EFICIENTE = %d\n",k,objval, eficaz);
        cantidad++;
    }
    else
        {
            fprintf(fsolucion,"DMU%d= %.10lf\t|||| EFICIENTE = %d\n",k,objval, eficaz);
            empresas[contador]=k;
            contador++;
        }

//fprintf(fsolucion,"DMU%d= %.10lf\t|||| EFICIENTE = %d\n",k,objval, eficaz);

eficaz=0;
fflush(stdout);


/*Libero memoria de CPLEX*/

status = CPXfreeprob(env, &lp);

fclose(fsolucion);

}

status = CPXcloseCPLEX (&env);

t2=time(0);

fsolucion=fopen("solucion.txt","a+");
fprintf(fsolucion,"\n\nEficientes = %d",cantidad);
fprintf(fsolucion,"\n-----------Tiempo calculo eficiencia = %f--------------\n",t2-t1);
fclose(fsolucion);

return contador;

}
