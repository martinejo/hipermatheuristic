#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ilcplex/cplex.h>
#include <time.h>

double timeval_diff(struct timeval *a, struct timeval *b)
{
  return
    (double)(a->tv_sec + (double)a->tv_usec/1000000) -
    (double)(b->tv_sec + (double)b->tv_usec/1000000);
}

void main(){

int i,j,k=0,status=0,CANTIDAD_N, CANTIDAD_M, CANTIDAD_S, M=10000000000000000, estado=0,eficaz=0, efi=1, noefi=1;

double **entrada,
    **salida;

double objval;

double t1,t2,t3,t4;

    struct timeval t_ini, t_fin;
	double secs;

FILE *fentradas, *fsalidas, *fsolucion;
CPXENVptr env=NULL;
CPXLPptr lp=NULL;


/* --------------------------------------------------------------------------------------------------*/
/* --------------------------------------------------------------------------------------------------*/
/* --------------------------------------------------------------------------------------------------*/
/* ----------------------------------------- INICIALIZACION------------------------------------------*/
/* --------------------------------------------------------------------------------------------------*/
/* --------------------------------------------------------------------------------------------------*/
/* --------------------------------------------------------------------------------------------------*/


fentradas=fopen("entradas.txt","r");
fsalidas=fopen("salidas.txt","r");
fsolucion=fopen("solucion.txt","a+");

/*
printf("\nNUMERO DE DMU's: ");
scanf("%d", &CANTIDAD_N);
printf("\nNUMERO DE ENTRADAS: ");
scanf("%d", &CANTIDAD_M);
printf("\nNUMERO DE SALIDAS: ");
scanf("%d", &CANTIDAD_S);
fflush(stdout);
*/

CANTIDAD_N=2500;
CANTIDAD_S=1;
CANTIDAD_M=2;

int *eficientes, *noeficiente;

eficientes = (int*)malloc((CANTIDAD_N+1)*sizeof(int));
noeficiente = (int*)malloc((CANTIDAD_N+1)*sizeof(int));

entrada = (double **)malloc((CANTIDAD_N+1)*sizeof(double*));
for (i=0;i<CANTIDAD_N+1;i++) {
		entrada[i] = (double*)malloc((CANTIDAD_M+1)*sizeof(double));
}

salida = (double **)malloc((CANTIDAD_N+1)*sizeof(double*));
for (i=0;i<CANTIDAD_N+1;i++) {
		salida[i] = (double*)malloc((CANTIDAD_S+1)*sizeof(double));
}

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

/////////////////// ESCRIBIR DATOS //////////////////////
/*
printf("\n --------------ENTRADAS------------------");
fprintf(fsolucion,"\n --------------ENTRADAS------------------");
for(i=1;i<CANTIDAD_N+1;i++){
        printf("\n DMU%d",i);
        fprintf(fsolucion,"\n DMU%d",i);
    for(j=1;j<CANTIDAD_M+1;j++){
        printf(" %f ", entrada[i][j]);
        fprintf(fsolucion," %f ", entrada[i][j]);
    }
}
printf("\n --------------SALIDAS------------------");
fprintf(fsolucion,"\n --------------SALIDAS------------------");
for(i=1;i<CANTIDAD_N+1;i++){
        printf("\n DMU%d",i);
        fprintf(fsolucion,"\n DMU%d",i);
    for(j=1;j<CANTIDAD_S+1;j++){
        printf(" %lf ", salida[i][j]);
        fprintf(fsolucion," %lf ", salida[i][j]);
    }
}
*/
fflush(stdout);

// --------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
// ----------------------------------------- CALCULO DE EFICIENCIA ----------------------------------
// --------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------


//fprintf(fsolucion,"\n\n ============ EFICIENTES ============= \n\n");

for(k=1;k<CANTIDAD_N+1;k++){

FILE*fmodelo;

fmodelo=fopen("modelo.lp","w");
fflush(stdout);
// CREAR MODELO A OPTIMIZAR

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
// CREAR RESTRICCIONES UTILIZADAS

fprintf(fmodelo,"\nst\n");

for(i=1;i<CANTIDAD_M+1;i++){
        fprintf(fmodelo,"sumatorio%d: ",i);
        for(j=1;j<CANTIDAD_N+1;j++){
            fprintf(fmodelo,"+%fd%d",entrada[j][i],j);
        }
       fprintf(fmodelo,"+sn%d=%f\n",i,entrada[k][i]);
}
fflush(stdout);

for(i=1;i<CANTIDAD_S+1;i++){
        fprintf(fmodelo,"sumatorio%d: ",i+CANTIDAD_M);
        for(j=1;j<CANTIDAD_N+1;j++){
            fprintf(fmodelo,"+%fd%d",salida[j][i],j);
        }
       fprintf(fmodelo,"-sp%d=%f\n",i,salida[k][i]);
}
fflush(stdout);
// CREAR COTAS DE VARIABLES

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


//Me conecto a CPLEX

env=CPXopenCPLEX(&status);

//status=CPXsetintparam(env,CPX_PARAM_SCRIND,CPX_ON);

status=CPXsetdblparam(env,CPX_PARAM_EPRHS,1e-09);
status=CPXsetintparam(env,CPX_PARAM_PREIND,1);

//Creo el problema lp

lp=CPXcreateprob(env,&status,"problema");

//Leo el modelo lp desde el fichero

status=CPXreadcopyprob(env,lp,"modelo.lp",NULL);

//Resuelvo el problema de programación -lineal-

status=CPXprimopt(env,lp);

//accedo al valor óptimo

status=CPXgetobjval(env,lp,&objval);

//Escribo el valor óptimo del problema
if(fabs(objval) < 0.00000001)
    {
        eficaz=1;
        //fprintf(fsolucion,"DMU%d= %.10lf\t|||| EFICIENTE = %d\n",k,objval, eficaz);
		eficientes[efi]=k;
		efi++;
    }
    else
        {
            //fprintf(fsolucion,"DMU%d= %.10lf\t|||| EFICIENTE = %d\n",k,objval, eficaz);
			noeficiente[noefi]=k;
			noefi++;
        }

//fprintf(fsolucion,"DMU%d= %.10lf\t|||| EFICIENTE = %d\n",k,objval, eficaz);

eficaz=0;
fflush(stdout);


//Libero memoria de CPLEX

status = CPXfreeprob(env, &lp);
status = CPXcloseCPLEX (&env);


}
k=0;

// --------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------
// --------------------------------------- CALCULO DEL MODELO----------------------------------------
// --------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------

t1=time(0);
gettimeofday(&t_ini, NULL);

FILE*fmodelo;

int calcula = 1;

fprintf(fsolucion,"\n================%d %d %d===================",CANTIDAD_N,CANTIDAD_M,CANTIDAD_S);

for(k=1;k<CANTIDAD_N+1;k++)
{
	
	if(noeficiente[calcula]==k)
	{
		//printf("\nEficientes[%d] = %d",calcula,eficientes[calcula]);
		fmodelo=fopen("modelo.lp","w");
		fflush(stdout);
		// CREAR MODELO A OPTIMIZAR

		fprintf(fmodelo,"\nmaximize\nfuncion: ");
		fflush(stdout);


		fprintf(fmodelo,"beta",(1.0/CANTIDAD_M));
		for(i=1;i<CANTIDAD_M+1;i++){
			fprintf(fmodelo,"-%.10lftn%d",(1.0/CANTIDAD_M)*(1.0/entrada[k][i]),i);
		}
		fflush(stdout);

		fprintf(fmodelo,"\nst");

		fprintf(fmodelo,"\nr1:");
		fprintf(fmodelo,"beta");
		for(i=1;i<CANTIDAD_S+1;i++){
			fprintf(fmodelo,"+%.10lftp%d",(1.0/CANTIDAD_S)*(1.0/salida[k][i]),i);
		}
		fprintf(fmodelo," = 1");
		fflush(stdout);

		for(i=1;i<CANTIDAD_M+1;i++){
			fprintf(fmodelo,"\nr2_%d:",i);
			for(j=1;j<efi;j++){
				fprintf(fmodelo,"+%falfa%d",entrada[eficientes[j]][i],eficientes[j]);
			}
			fprintf(fmodelo,"+tn%d-%fbeta = 0",i,entrada[k][i]);
		}
		fflush(stdout);

		for(i=1;i<CANTIDAD_S+1;i++){
			fprintf(fmodelo,"\nr3_%d:",i);
			for(j=1;j<efi;j++){
				fprintf(fmodelo,"+%falfa%d",salida[eficientes[j]][i],eficientes[j]);
			}
			fprintf(fmodelo,"-tp%d-%fbeta = 0",i,salida[k][i]);
		}
		fflush(stdout);

		for(i=1;i<efi;i++){
			fprintf(fmodelo,"\nr4_%d:",i);
			for(j=1;j<CANTIDAD_S+1;j++){
				fprintf(fmodelo,"+%fu%d",salida[eficientes[i]][j],j);
			}
			for(j=1;j<CANTIDAD_M+1;j++){
				fprintf(fmodelo,"-%fv%d",entrada[eficientes[i]][j],j);
			}
			fprintf(fmodelo,"+d%d = 0",eficientes[i]);
		}
		fflush(stdout);

		for(i=1;i<efi;i++)
			{	
			fprintf(fmodelo,"\nr7_%d:",eficientes[i]);		
			fprintf(fmodelo,"+d%d -100000000000000000b%d <= 0",eficientes[i],eficientes[i]);
			}
		fflush(stdout);

		for(i=1;i<efi;i++)
			{	
			fprintf(fmodelo,"\nr8_%d:",eficientes[i]);		
			fprintf(fmodelo,"+alfa%d +100000000000000000b%d <= 100000000000000000",eficientes[i],eficientes[i]);
			}
		fflush(stdout);



		fprintf(fmodelo,"\nbounds\n");

		for(i=1;i<CANTIDAD_M+1;i++){
			fprintf(fmodelo,"v%d>=1\n",i);
		}
		for(i=1;i<CANTIDAD_S+1;i++){
			fprintf(fmodelo,"u%d>=1\n",i);
		}


		fprintf(fmodelo,"beta>=0\n");
		for(i=1;i<CANTIDAD_M+1;i++){
			fprintf(fmodelo,"tn%d>=0\n",i);
		}
		for(i=1;i<CANTIDAD_S+1;i++){
			fprintf(fmodelo,"tp%d>=0\n",i);
		}
		for(i=1;i<efi;i++){
			fprintf(fmodelo,"d%d>=0\n",eficientes[i]);
		}
		for(i=1;i<efi;i++){
			fprintf(fmodelo,"alfa%d>=0\n",eficientes[i]);
		}

		fprintf(fmodelo,"\nBINARY\n");
		for(i=1;i<efi;i++)
			{		
			fprintf(fmodelo,"b%d\n",eficientes[i]);
			}
		fflush(stdout);
		//fprintf(fmodelo,"SOS\n");
		//for(i=1;i<CANTIDAD_N+1;i++)
		//{
		//    fprintf(fmodelo,"set%d: S1:: d%d:1 alfa%d:2\n",i-1,i,i);
		//}
		fprintf(fmodelo,"end");
		fclose(fmodelo);



		//Me conecto a CPLEX

		env=CPXopenCPLEX(&status);

		//status=CPXsetintparam(env,CPX_PARAM_SCRIND,CPX_ON);
		//status = CPXsetintparam (env, CPX_PARAM_SCRIND, CPX_OFF);
		status=CPXsetdblparam(env,CPX_PARAM_NODELIM,1);
		status=CPXsetintparam(env,CPX_PARAM_PREIND,1);
		//QUITAR PRESOLVER EN LA EJECUCION

		//Creo el problema lp

		lp=CPXcreateprob(env,&status,"problema");

		//Leo el modelo lp desde el fichero

		status=CPXreadcopyprob(env,lp,"modelo.lp",NULL);

		//Resuelvo el problema de programación -lineal-

		status=CPXmipopt(env,lp);

		//accedo al valor óptimo

		status=CPXgetmipobjval(env,lp,&objval);

		estado = CPXgetstat(env,lp);

			if(objval<1.000 && estado<=102){
				//fprintf(fsolucion,"\n");
				fprintf(fsolucion,"\nEmpresa\t%d:\t%.10lf\t||\tEstado\t=\t%d",k,objval,estado);
				printf("\nEmpresa %d: %.10lf || Estado = %d",k,objval,estado);
			}else
			{
				//fprintf(fsolucion,"\n");
				fprintf(fsolucion,"\nEmpresa\t%d:\t0.0\t||\tEstado\t=\t%d",k,estado);
				printf("\nEmpresa %d: 0.0 || Estado = %d",k,estado);
			}

		estado=0;
		status = CPXfreeprob(env, &lp);
		status = CPXcloseCPLEX (&env);
		
		calcula++;
	}
	
	

}



	for(i=0;i<CANTIDAD_N+1;i++){
	free(entrada[i]);
	}
	free(entrada);
	for(i=0;i<CANTIDAD_N+1;i++){
	free(salida[i]);
	}
	free(salida);

	

	t2=time(0);
	gettimeofday(&t_fin, NULL);
	secs = timeval_diff(&t_fin, &t_ini);
	printf("\n\nEl proceso a tardado %f\n",t2-t1);
	
	fprintf(fsolucion,"\nTiempo = %f",t2-t1);
	fclose(fsolucion);



}



