#include "genetico.h"
#include <stdlib.h>
#include <stdio.h>

extern void dgesv_( int* n, int* nrhs, double* a, int* lda, int* ipiv,
                double* b, int* ldb, int* info );

/**
 * Crea un cromosoma. El cromosoma cumplirá las restricciones
 * @param etiq, etiqueta del cromosoma.
 * @param *inputs, vector con las entradas
 * @param *outputs, vector con las salidas
 * @param col, variable auxiliar para fijar alguna columna
 * @param *Vx_columna, vector con los sumatorios en columnas de las entradas
 * @param *Vy_columna, vector con los sumatorios en columnas de las salidas
 */

Ptr_Chromosome create_chromosome( int etiq, double *inputs, double *outputs, int col, double *Vx_columna, double *Vy_columna, double *inputs_inv, double* outputs_inv, double* sumatorio_entradas, double* sumatorio_salidas, int rellenar )
{
    int maxiteraciones = 0,
        iteraciones = 0,
        metodo = 0,
        indicador_metodo=0,
        indicador_valido=1,
        alfacero=0,
        opcion=CANTIDAD_M,
        alfamenor=0;
    double aux = 1.0,
//         tiempo,
           precision = 0.000001,
           *Vx_fila;
//    struct timeval ti, tf;

    Ptr_Chromosome chrom = NULL;

    chrom = (Ptr_Chromosome)malloc(sizeof(struct Chromosome));
    Vx_fila = (double*)malloc(CANTIDAD_M*sizeof(double));

    //gettimeofday(&ti, NULL);
    chrom->etiq = etiq;
    chrom->metodo = 0;
    chrom->valido = 0;
    chrom->fi=0;
    chrom->evaluation = BAD_CHROM;
    chrom->vik = (double*)malloc(CANTIDAD_M*sizeof(double));
    chrom->urk = (double*)malloc(CANTIDAD_S*sizeof(double));
    chrom->tneg = (double*)malloc(CANTIDAD_M*sizeof(double));
    chrom->tpos = (double*)malloc(CANTIDAD_S*sizeof(double));
    chrom->djk = (double*)malloc(CANTIDAD_N*sizeof(double));
    chrom->alfas = (double*)malloc(CANTIDAD_N*sizeof(double));
    chrom->discrete_model = (int*)malloc(CANTIDAD_N*sizeof(int));

    chrom->contador=0;
    chrom->terminar=0;

    suma_filas(chrom,inputs,Vx_fila);

    while(aux!=0.0&&maxiteraciones<((25*rellenar)+1))
    {
        //Primer método
        if(metodo==0&&chrom->valido==0){
            indicador_metodo = 1;
            opcion=(rand()%(CANTIDAD_M-CANTIDAD_S+1))+CANTIDAD_S;
            do
            {
            generate_b(chrom,opcion, sumatorio_entradas, sumatorio_salidas);
            generate_beta(chrom, outputs,col);
            calculo_alfas(chrom,outputs,col);
            modify_constraint_1(chrom, inputs, col);
            alfacero++;
            }
            while(check_constraint_14(chrom)&&alfacero<50);
            alfacero=0;

            deducir_valores_uv(chrom, inputs, outputs, col);
            //if(chrom->contador>=2000){maxiteraciones=101;}

            //Comprobamos si cumple con las restricciones
            aux=check_constraint_1(chrom, outputs, col, precision,outputs_inv);
            aux+=check_constraint_11(chrom);
            aux+=check_constraint_12(chrom);
            aux+=check_constraint_2(chrom, inputs, col, precision);
            aux+=check_constraint_3(chrom, outputs, col, precision);
            aux+=check_constraint_4(chrom, inputs, outputs, col, precision);
            aux+=check_constraint_5(chrom);
            aux+=check_constraint_6(chrom);
            aux+=check_constraint_13(chrom);
            aux+=check_constraint_14(chrom);
            aux+=check_constraint_15(chrom);
            if(aux==0.0)
            {
                chrom->valido=1;
                chrom->metodo = 1;
                alfacero=0;
            }
            iteraciones++;

            if(iteraciones==((2999*rellenar)+1)){
                /*opcion-=1;
                if(opcion<CANTIDAD_S)
                    {
                        //maxiteraciones = 101;
                        metodo=1;
                        opcion=CANTIDAD_M;
                }
                iteraciones=0;*/
                //show_chromosome(chrom,inputs,outputs,col,inputs_inv,outputs_inv);
                metodo=1;
            }
        }else{ if(metodo==1&&chrom->valido==0){
            indicador_metodo = 2;
            //Aquí calculamos las U y las V
            while(indicador_valido!=0&&alfacero<10){
            chrom->beta = Random();
            //Cálculo de los valores discretos
            generate_b(chrom,opcion, sumatorio_entradas, sumatorio_salidas);
            //Generamos las alfas llamando a la función
            generate_alfas(chrom, inputs, Vx_fila);
            alfamenor= mejor_caso(chrom, inputs, outputs);
            indicador_valido=modify_alfas(chrom, inputs, outputs, col, alfamenor, Vx_columna, Vy_columna,outputs_inv);
            alfacero++;
            chrom->contador++;
            //printf("\nOPCION = %d",opcion);
            if(chrom->contador == 250)
                {
                    chrom->terminar=1;
                }
            }
            alfacero=0;
            deducir_valores_uv(chrom, inputs, outputs, col);

            //Comprobamos si cumple con las restricciones
            aux=check_constraint_1(chrom, outputs, col, precision,outputs_inv);
            aux+=check_constraint_11(chrom);
            aux+=check_constraint_12(chrom);
            aux+=check_constraint_2(chrom, inputs, col, precision);
            aux+=check_constraint_3(chrom, outputs, col, precision);
            aux+=check_constraint_4(chrom, inputs, outputs, col, precision);
            aux+=check_constraint_5(chrom);
            aux+=check_constraint_6(chrom);
            aux+=check_constraint_13(chrom);
            aux+=check_constraint_14(chrom);
            aux+=check_constraint_15(chrom);
            maxiteraciones++;
            /*if(maxiteraciones==25){
                opcion-=1;
                maxiteraciones=0;
                if(opcion<CANTIDAD_S){maxiteraciones=101;}
            }*/
            if(aux==0.0)
            {
                chrom->valido=1;
                chrom->metodo=2;
            }
            }
        }
    }
    //gettimeofday(&tf, NULL);   // Instante final
    //tiempo= (tf.tv_sec - ti.tv_sec)*1000 + (tf.tv_usec - ti.tv_usec)/1000.0;
    printf(") CONTADOR = %d",chrom->contador);
    printf("||| Solucion: %d METODO %d b = %d |||", chrom->valido, indicador_metodo, opcion);
    //if(chrom->contador >= 2500){chrom->terminar=1;}
    fflush(stdout);
    //show_chromosome(chrom,inputs,outputs,col,inputs_inv,outputs_inv);
    //system("pause");
    free(Vx_fila);

    return chrom;
}

/**
 * Genera las b del cromosoma.
 */
void generate_b(Ptr_Chromosome chrom, int opcion, double *sumatorio_entradas, double *sumatorio_salidas)
{
    int i = 0,
    w = 0,
    repeticion=0,
    aux=opcion;
    //aux =(rand()%(CANTIDAD_M-CANTIDAD_S+1))+CANTIDAD_S;

    for(i=0;i<CANTIDAD_N;i++)
    {
         chrom->discrete_model[i]=1;
    }

    for(i=0;i<aux;i++)
    {
        w=rand()%CANTIDAD_N;
        while(chrom->discrete_model[w]==0&&repeticion<CANTIDAD_N)
            {
            w=rand()%CANTIDAD_N;
            repeticion++;
            }
        chrom->discrete_model[w]=0;
        repeticion=0;
    }
}

/**
 * @brief Función para generar las alfas
 * @param Chrom , es el cromosoma
 * @param *inputs, vector con las entradas
 * @param *Vx_fila, vector con los sumatorios en filas de las entradas
 */
void generate_alfas (Ptr_Chromosome chrom, double *inputs, double *Vx_fila)
{
    int i = 0,
        j = 0,
        alfasnocero = 0;
    double aux = 0.0;

    //Sacamos la posición de la fila del sumatorio mayor que es el peor caso
    //Además sacamos el número de alfas distintas de 0
    j=0;
    for(i=0;i<CANTIDAD_M;i++)
    {
        if(Vx_fila[j]<Vx_fila[i])j=i;
    }

    for(i=0;i<CANTIDAD_N;i++)
    {
        if(chrom->discrete_model[i]==0)alfasnocero += 1;
    }

    //Saber cuanto tiene que valer cada termino dintinto de 0 para que el sumatorio valga 0.5
    aux=1.0/alfasnocero;

    //Calculo de número de alfas distintas de 0 para forzar el peor caso que el sumatorio sea 0.5
    for(i=0;i<CANTIDAD_N;i++)
    {
        if(chrom->discrete_model[i]==1){
            chrom->alfas[i] = 0.0;
            chrom->djk[i] = Random();
        }else{
            chrom->djk[i] = 0.0;
            chrom->alfas[i]=aux/inputs[i+(j*CANTIDAD_N)];
        }
    }
    /*for(i=0;i<CANTIDAD_N;i++){
        printf("\nAlfa[%d] = %f",i,chrom->alfas[i]);
    }
    system("pause");*/
}


/**
 * @brief Función para modificar los alfas para poder cumplir las restricciones 1.2 y 1.3
 * Se limita el valor de los sumatorios a los valores de entrada y salida (ya que beta <=1)
 * @param Chrom , es el cromosoma
 * @param *inputs, vector con las entradas
 * @param *outputs, vector con las salidas
 * @param col, variable auxiliar para fijar alguna columna
 * @param alfamenor, indice del vector que contiene las alfas el cual apunta al alfa mas eficiente de modificar para satisfacer la restriccion 1
 * @param *Vx_columna, vector con los sumatorios en columnas de las entradas
 * @param *Vy_columna, vector con los sumatorios en columnas de las salidas
 */
int modify_alfas (Ptr_Chromosome chrom, double *inputs, double *outputs, int col,
                    int alfamenor,double *Vx_columna, double *Vy_columna, double *outputs_inv) {
    int i = 0,
        columna,
        d=0;

    double *Vx,
           *Vy;
    double valido = 1.0,
           aux1 = 0.01,
           precision = 0.000001,
           aux2 = 0.01;


    //Sumatorios alfa*entrada y alfa*salida
    Vx = (double*)malloc(CANTIDAD_M*sizeof(double));
    Vy = (double*)malloc(CANTIDAD_S*sizeof(double));

    sumatorios(chrom,inputs,outputs,Vx,Vy);

    //Forzar que Xik de la restricción 2 no sea menor nunca que el sumatorio (ya que el rango de beta es {0,1})
    while(valido!=0&&d<20)
    {
        valido=0;

        for(i=0;i<CANTIDAD_M;i++)
        {
            while(Vx[i]>(inputs[i*CANTIDAD_N+col]-(inputs[i*CANTIDAD_N+col]*aux1)))
            {
                columna=termino_inout(chrom,inputs,Vy_columna,i);

                chrom->alfas[columna]=chrom->alfas[columna]*0.95;
                if(chrom->alfas[columna]<0.000003)chrom->alfas[columna]=0.0;

                sumatorios(chrom,inputs,outputs,Vx,Vy);
            }
        }
        aux1+=0.03;

        for(i=0;i<CANTIDAD_S;i++)
        {
            while(Vy[i]<(outputs[i*CANTIDAD_N+col]+(outputs[i*CANTIDAD_N+col]*aux2)))
            {
                columna=termino_inout(chrom,outputs,Vx_columna,i);

                chrom->alfas[columna]=chrom->alfas[columna]*1.05;

                sumatorios(chrom,inputs,outputs,Vx,Vy);
            }
        }
        aux2+=0.03;
        deducir_valores_alfa(chrom, outputs, col, alfamenor,outputs_inv, Vx_columna, Vy_columna);
        modify_beta(chrom, inputs, outputs, col, alfamenor, Vx_columna, Vy_columna,outputs_inv);
        modify_constraint_1(chrom, inputs, col);
        modify_constraint_2(chrom, outputs, col);
        valido+=check_constraint_1(chrom, outputs, col, precision,outputs_inv);
        //printf("\n -------------- Valido 1 = %f ------------", valido);
        valido+=check_constraint_2(chrom, inputs, col, precision);
        //printf("\n -------------- Valido 2 = %f ------------", valido);
        valido+=check_constraint_3(chrom, outputs, col, precision);
        //printf("\n -------------- Valido 3 = %f ------------", valido);
        valido+=check_constraint_11(chrom);
        //printf("\n -------------- Valido 11 = %f ------------", valido);
        valido+=check_constraint_12(chrom);
        //printf("\n -------------- Valido 12 = %f ------------", valido);
        valido+=check_constraint_14(chrom);
        //printf("\n -------------- Valido 14 = %f ------------", valido);
        //system("pause");
        d++;
    }
    //if(valido==0)    printf("\nSalgo???");

    free(Vy);
    free(Vx);
    return (valido);
}

/**
 * @brief Función para modificar los alfas para poder cumplir las restricciones 1.2 y 1.3
 * Se limita el valor de los sumatorios a los valores de entrada y salida (ya que beta <=1)
 * @param Chrom , es el cromosoma
 * @param *inputs, vector con las entradas
 * @param *outputs, vector con las salidas
 * @param col, variable auxiliar para fijar alguna columna
 * @param *Vx_columna, vector con los sumatorios en columnas de las entradas
 * @param *Vy_columna, vector con los sumatorios en columnas de las salidas
 */
void modify_alfas_2 (Ptr_Chromosome chrom, double *inputs, double *outputs, int col, double *Vx_columna, double *Vy_columna)
{
    int i = 0,
        columna;
    double *Vx,
           *Vy;
    double aux1=0.01,
            aux2=0.01;

    //Sumatorios alfa*entrada y alfa*salida
    Vx = (double*)malloc(CANTIDAD_M*sizeof(double));
    Vy = (double*)malloc(CANTIDAD_S*sizeof(double));

    sumatorios(chrom,inputs,outputs,Vx,Vy);

    //Forzar que Xik de la restricción 2 no sea menor nunca que el sumatorio (ya que el rango de beta es {0,1})
    for(i=0;i<CANTIDAD_M;i++)
    {
        while(Vx[i]>(inputs[i*CANTIDAD_N]-(inputs[i*CANTIDAD_N]*aux1)))
        {
            columna=termino_inout(chrom,inputs,Vy_columna,i);

            chrom->alfas[columna]=chrom->alfas[columna]*0.97;
            if(chrom->alfas[columna]<0.000003)chrom->alfas[columna]=0.0;

            sumatorios(chrom,inputs,outputs,Vx,Vy);
        }
    }
    aux1+=0.03;

    for(i=0;i<CANTIDAD_S;i++)
    {
        while(Vy[i]<(outputs[i*CANTIDAD_N]+(outputs[i*CANTIDAD_N]*aux2)))
        {

            columna=termino_inout(chrom,outputs,Vx_columna,i);

            chrom->alfas[columna]=chrom->alfas[columna]*1.03;

            sumatorios(chrom,inputs,outputs,Vx,Vy);
        }
    }
    aux2+=0.03;

    free(Vx);
    free(Vy);
}


/////////////////////////////////////////////////////////////////////////////////////////////
//                                      RESTRICCIONES                                      //
/////////////////////////////////////////////////////////////////////////////////////////////

/**
 * @brief Comprueba que se cumple con la ecuación de la restricción 1
 * @param Chrom , es el cromosoma
 * @param *outputs, matriz con las salidas
 * @param col, variable auxiliar para fijar alguna columna
 * @param precision, variable con la que eliminaremos errores en los decimales residuales
 */

double check_constraint_1 (Ptr_Chromosome chrom, double *outputs, int col, double precision, double *outputs_inv)
{
    int i;
    double aux = 0.0,
    c_1 = 0.0;

    // Error de la restriccion c.1
    aux=0;
    for ( i = 0 ; i < CANTIDAD_S ; i++ )
    {
        aux += chrom->tpos[i]/outputs[(i*CANTIDAD_N)+col];
    }
    aux /= CANTIDAD_S;
    aux += chrom->beta;
    c_1 += 1-aux;
    if (((c_1 - precision)< 0.0) && ((c_1 + precision) > 0.0))
    {
            c_1=0.0;
    }
    return c_1;
}

/**
 * @brief Comprueba que se cumple con la ecuación de la restricción 2
 * @param Chrom , es el cromosoma
 * @param *inputs, matriz con las entradas
 * @param col, variable auxiliar para fijar alguna columna
 * @param precision, variable con la que eliminaremos errores en los decimales residuales
 */

double check_constraint_2 (Ptr_Chromosome chrom, double *inputs, int col, double precision)
{
    int i,k;
    double *aux,
           c_2 = 0.0;

    aux = (double*)malloc(CANTIDAD_M*sizeof(double));

    for(i=0;i<CANTIDAD_M;i++){
            aux[i]=0.0;
        for(k=0;k<CANTIDAD_N;k++){
            aux[i]+=chrom->alfas[k]*inputs[(i*CANTIDAD_N)+k];
        }
    }

    // Error de la restriccion c.2
    for(i=0;i<CANTIDAD_M;i++){
        aux[i]-=chrom->beta*inputs[(i*CANTIDAD_N)+col]-chrom->tneg[i];
    }

    for(i=0;i<CANTIDAD_M;i++){
        c_2+=aux[i];
    }

    if (((c_2 - precision)< 0.0) && ((c_2 + precision) > 0.0))
    {
            c_2=0.0;
    }
    free(aux);
    return c_2;
}

/**
 * @brief Comprueba que se cumple con la ecuación de la restricción 3
 * @param Chrom , es el cromosoma
 * @param *outputs, matriz con las salidas
 * @param col, variable auxiliar para fijar alguna columna
 * @param precision, variable con la que eliminaremos errores en los decimales residuales
 */

double check_constraint_3 (Ptr_Chromosome chrom, double *outputs, int col, double precision)
{
    int i,k;
    double *aux,
    c_3 = 0.0;

    aux = (double*)malloc(CANTIDAD_S*sizeof(double));

    // Error asociado a c.3
    for(i=0;i<CANTIDAD_S;i++){
            aux[i]=0.0;
        for(k=0;k<CANTIDAD_N;k++){
            aux[i]+=chrom->alfas[k]*outputs[(i*CANTIDAD_N)+k];
        }
    }

    for(i=0;i<CANTIDAD_S;i++){
        aux[i]-=chrom->beta*outputs[(i*CANTIDAD_N)+col]+chrom->tpos[i];
    }

    for(i=0;i<CANTIDAD_S;i++){
        c_3+=aux[i];
    }

    if (((c_3- precision)< 0.0) && ((c_3 + precision) > 0.0))
    {
            c_3=0.0;
    }

    free(aux);
    return c_3;
}


/**
 * @brief Comprueba que se cumple con la ecuación de la restricción 4
 * @param Chrom , es el cromosoma
 * @param *inputs, vector con las entradas
 * @param *outputs, vector con las salidas
 * @param col, variable auxiliar para fijar alguna columna
 * @param precision, variable con la que eliminaremos errores en los decimales residuales
 */
double check_constraint_4(Ptr_Chromosome chrom, double *inputs, double* outputs, int col, double precision)
{
    int i,k;
    double *aux,
           *aux2,
           c_4 = 0.0;

    aux = (double*)malloc(CANTIDAD_N*sizeof(double));
    aux2 = (double*)malloc(CANTIDAD_N*sizeof(double));

    // Error asociado a c.4
    for(i=0;i<CANTIDAD_N;i++){
            aux[i]=0.0;
        for(k=0;k<CANTIDAD_M;k++){
            aux[i]+=chrom->vik[k]*inputs[(k*CANTIDAD_N)+i];
        }
    }

    for(i=0;i<CANTIDAD_N;i++){
            aux2[i]=0.0;
        for(k=0;k<CANTIDAD_S;k++){
            aux2[i]+=chrom->urk[k]*outputs[(k*CANTIDAD_N)+i];
        }
    }

    for(i=0;i<CANTIDAD_N;i++){
        aux2[i]-=aux[i];
        aux2[i]+=chrom->djk[i];
    }

    for(i=0;i<CANTIDAD_N;i++){
        c_4+=aux2[i];
    }

    if (((c_4- precision)< 0.0) && ((c_4 + precision) > 0.0))
    {
            c_4=0.0;
    }

    free(aux);
    free(aux2);
    return c_4;
}


/**
 * @brief Comprueba que las 'v' son todas mayores o igual a 1
 * @param Chrom , es el cromosoma
 */

int check_constraint_5 (Ptr_Chromosome chrom)
{
    int valido = 0,
        i = 0;

    while (i < CANTIDAD_M && !valido)
    {
        if ( chrom->vik[i] < 1 )
            valido = 1;
        i++;
    }

    return valido;
}


/**
 * @brief Comprueba que las 'u' son todas mayores o igual a 1
 * @param Chrom , es el cromosoma
 */
int check_constraint_6 (Ptr_Chromosome chrom)
{
    int valido = 0,
    i = 0;

    while (i < CANTIDAD_S && !valido)
    {
        if ( chrom->urk[i] < 1 )
            valido = 1;
        i++;
    }

    return valido;
}


/**
 * @brief Comprueba que las t- son todas mayores o igual a 0
 * @param Chrom , es el cromosoma
 */
int check_constraint_11 (Ptr_Chromosome chrom)
{
    int valido = 0,
    i = 0;

    while (i < CANTIDAD_M && !valido)
    {
        if ( chrom->tneg[i] < 0 )
            valido = 1;
        //printf("\nTneg[%d] = %f",i,chrom->tneg[i]);
        i++;

    }

    return valido;
}


/**
 * @brief Comprueba que las t+ son todas mayores o igual a 0
 * @param Chrom , es el cromosoma
 */
int check_constraint_12 (Ptr_Chromosome chrom)
{
    int valido = 0,
    i = 0;

    while (i < CANTIDAD_S && !valido)
    {
        if ( chrom->tpos[i] < 0 )
            valido = 1;
        i++;
    }

    return valido;
}

/**
 * @brief Comprueba que las d son todas mayores o igual a 0
 * @param Chrom , es el cromosoma
 */
int check_constraint_13 (Ptr_Chromosome chrom)
{
    int valido = 0,
    i = 0;

    while (i < CANTIDAD_N && !valido)
    {
        if ( chrom->djk[i] < 0 )
            valido = 1;
        i++;
    }

    return valido;
}
/*Alfas todas positivas!!!!!*/

int check_constraint_14 (Ptr_Chromosome chrom)
{
    int valido = 0,
    i = 0;

    while (i < CANTIDAD_N && !valido)
    {
        if ( chrom->alfas[i] < 0 )
            valido = 1;
        i++;
    }

    return valido;
}

int check_constraint_15 (Ptr_Chromosome c)
{
    int valido = 0,
    i = 0;

     for(i=0;i<CANTIDAD_N;i++){
        if(c->discrete_model[i]==0){
            if(c->djk[i]!=0){
                valido=1;
            }
        }
    }

    return valido;
}
/**
 * Hace una pre-evaluacion del cromosoma para saber como de lejos está de ser válido.
 * @param Chrom , es el cromosoma
 * @param *inputs, vector con las entradas
 * @param *outputs, vector con las salidas
 * @param col, variable auxiliar para fijar alguna columna
 */

double pre_evaluate_chromosome (Ptr_Chromosome chrom, double *inputs, double *outputs, int col, double *inputs_inv, double *outputs_inv)
{
    int i;
    double aux = 0.0, precision=0.000001;

    //aux += check_constraint_1(chrom, outputs, col, precision,outputs_inv);
    //aux += check_constraint_2(chrom, inputs, col, precision);
    //aux += check_constraint_3(chrom, outputs, col, precision);
    //aux += check_constraint_4(chrom, inputs, outputs, col, precision);
    //aux += check_constraint_5(chrom);
    //aux += check_constraint_6(chrom);
    //aux += check_constraint_11(chrom);
    //aux += check_constraint_12(chrom);
    //aux += check_constraint_13(chrom);

    //if ( aux < 0 )
        //aux *= -1;
    for(i=0;i<CANTIDAD_M;i++){
        if(chrom->tneg[i]<0){
            aux+=fabs(chrom->tneg[i]);
        }
        if(chrom->vik[i]<1){
            aux+=fabs(chrom->vik[i]);
        }
    }
    for(i=0;i<CANTIDAD_S;i++){
        if(chrom->tpos[i]<0){
            aux+=fabs(chrom->tpos[i]);
        }
        if(chrom->urk[i]<1){
            aux+=fabs(chrom->urk[i]);
        }
    }
    for(i=0;i<CANTIDAD_N;i++){
        if(chrom->djk[i]<0){
            aux+=fabs(chrom->djk[i]);
        }
    }
    aux+=fabs(check_constraint_1(chrom,outputs,col,precision,outputs_inv));
    aux+=fabs(check_constraint_2(chrom,inputs,col,precision));
    aux+=fabs(check_constraint_3(chrom,outputs,col,precision));
    aux+=check_constraint_4(chrom,inputs,outputs,col,precision);

    return aux;
}


/**
 * Libera la memoria ocupada por el cromosoma.
 * @param c Cromosoma.
 */

void free_chromosome(Ptr_Chromosome c)
{
	free(c->discrete_model);
    free(c->tpos);
    free(c->tneg);
    free(c->vik);
    free(c->urk);
    free(c->djk);
    free(c->alfas);
	free(c);
	return;
}



/**
 * Muestra un cromosoma por pantalla
 * @param c Cromosoma.
 * @param *inputs, vector con las entradas
 * @param *outputs, vector con las salidas
 * @param col, variable auxiliar para fijar alguna columna
 */

void show_chromosome( Ptr_Chromosome c, double *inputs, double *outputs,int col,double *inputs_inv, double *outputs_inv )
{
    double precision=0.000001;
    int i,aux=0;

    printf("\n\t\t=================================");
    fflush(stdout);
    printf("\n\t\t---------Cromosoma---------\n");
    fflush(stdout);

    printf("\n\n\t\t====> Valor funcion maximizar = %f\n\n",maximizar(c,inputs,col));
    fflush(stdout);
    printf("\t\tEmpresa : %d",col+1);
    fflush(stdout);
    printf("\n\t\tMejorado mediante METODO = %d",c->fi);
    fflush(stdout);
    printf("\n\t\tVALIDO = %d",c->valido);
    printf("\n");
    fflush(stdout);

    for(i=0;i<CANTIDAD_N;i++){
        if(c->discrete_model[i]==0){
            /*if(c->alfas[i]==0){
                aux+=10;
            }*/
            if(c->djk[i]!=0){
                aux+=10;
            }
        }
    }

    for(i=0;i<CANTIDAD_N;i++){
        if(c->discrete_model[i]==0){
        printf("\n\t\tb[%d] = %d",i,c->discrete_model[i]);
        fflush(stdout);
        /*contador++;
        if(contador>=3){system("pause");}*/
        }
    }
           /*for(i=0;i<CANTIDAD_S;i++){
           printf("\nURK[%d] = %f",i,c->urk[i]);
           }
           printf("\n");
           for(i=0;i<CANTIDAD_M;i++){
           printf("\nVIK[%d] = %f",i,c->vik[i]);
           }
           printf("\n");*/
           /*for(i=0;i<CANTIDAD_N;i++){
           printf("\nDJK[%d] = %f",i,c->djk[i]);
           }
            printf("\n");*/
            for(i=0;i<CANTIDAD_S;i++){
            printf("\n\t\tTpos[%d] = %f",i,c->tpos[i]);
            }
            fflush(stdout);
            for(i=0;i<CANTIDAD_N;i++){
                if(c->alfas[i]!=0.0){
                    printf("\n\t\tAlfa[%d] = %f",i,c->alfas[i]);
                    printf("\n\t\tD[%d] = %f",i,c->djk[i]);
                }
            }
            fflush(stdout);
            /*
            printf("\n");
            for(i=0;i<CANTIDAD_M;i++){
            printf("\nTneg[%d] = %f",i,c->tneg[i]);
            }
            printf("\n");*/
            /*for(i=0;i<CANTIDAD_N;i++){
            printf("\nAlfa[%d] = %f",i,c->alfas[i]);
            }*/
/*
    printf("ETIQ: %4d\t", c->etiq);*/
    printf("\n\t\tEVAL:: %lf\t", c->evaluation);
    fflush(stdout);
/*  printf("PRE-EVAL: %lf\t", c->pre_evaluation);

    printf("DIF.ALFA: %f\t",difference(c->alfas, 0, CANTIDAD_N, 1.0));
    printf("\nBETA: %f\n",c->beta);*/
    fflush(stdout);

    printf("\n\t\tRestriccion 1: %f",check_constraint_1(c,outputs,col,precision,outputs_inv));
    fflush(stdout);
    printf("\n\t\tRestriccion 2: %f",check_constraint_2(c,inputs,col,precision));
    fflush(stdout);
    printf("\n\t\tRestriccion 3: %f",check_constraint_3(c,outputs,col,precision));
    fflush(stdout);
    printf("\n\t\tRestriccion 4: %f",check_constraint_4(c,inputs,outputs,col,precision));
    fflush(stdout);
    printf("\n\t\tRestriccion 5: %d",check_constraint_5(c));
    fflush(stdout);
    printf("\n\t\tRestriccion 6: %d",check_constraint_6(c));
    fflush(stdout);
    printf("\n\t\tRestriccion 11: %d",check_constraint_11(c));
    fflush(stdout);
    printf("\n\t\tRestriccion 12: %d",check_constraint_12(c));
    fflush(stdout);
    printf("\n\t\tRestriccion 13: %d",check_constraint_13(c));
    fflush(stdout);
    printf("\n\t\tRestriccion 14: %d",check_constraint_14(c));
    fflush(stdout);
    printf("\n\t\tRestriccion EXTRA: %d",aux);
    fflush(stdout);
    printf("\n\t\t=================================");
    printf("\n");
    fflush(stdout);
    return;
}



/**
 * Determina si un cromosoma es valido
 * @param c Cromosoma
 * @param *inputs, vector con las entradas
 * @param *outputs, vector con las salidas
 * @param col, variable auxiliar para fijar alguna columna
 */

int good_chromosome(Ptr_Chromosome c, double *inputs, double *outputs, int col, double *inputs_inv, double *outputs_inv)//, int NumTareas, Ptr_Tarea t, Ptr_Procesador)
{
    double precision=0.000001;
    int i;

    for(i=0;i<CANTIDAD_N;i++){
        if(c->discrete_model[i]==0){
            /*if(c->alfas[i]==0){
                return 0;
            }*/
            if(c->djk[i]!=0){
                return 0;
            }
        }
    }
    // Los cromosomas ya cumplen ciertas restriccione al generarlos y cuando la preevaluacion es igual a 0. Aquí se comprueba el resto.
    if (check_constraint_1(c,outputs,col,precision,outputs_inv)== 0.0&&
        check_constraint_2(c,inputs,col,precision)==0.0&&
        check_constraint_3(c,outputs,col,precision)==0.0&&
        check_constraint_4(c,inputs,outputs,col,precision)==0.0&&
        check_constraint_5(c) == 0 &&
        check_constraint_6(c) == 0 &&
        check_constraint_11(c) == 0 &&
        check_constraint_12(c) == 0 &&
        check_constraint_14(c) == 0 &&
        check_constraint_13(c) == 0)
        {
        return 1;
        }
    return 0;
}

/**
 * Evalua un cromosoma. El cromosoma tendra un valor menor cuanto mejor sea.
 * @param c Cromosoma.
 * @param *inputs, vector con las entradas
 * @param col, variable auxiliar para fijar alguna columna
 */

double evaluate_chromosome(Ptr_Chromosome chrom, double *inputs, double *outputs, int col, double *inputs_inv, double *outputs_inv)
{
    int i = 0;
    double aux = 0.0, precision=0.000001;


    for(i=0;i<CANTIDAD_N;i++){
        if(chrom->discrete_model[i]==0){
            /*if(chrom->alfas[i]==0){
                aux=10.0;
                return aux;
            }*/
            if(chrom->djk[i]!=0){
                aux=10.0;
                return aux;
            }
        }
    }
    if(check_constraint_1(chrom,outputs,col,precision,outputs_inv)!=0.0){
        aux=1.0;
        return aux;
    }
    if(check_constraint_2(chrom,inputs,col,precision)!=0.0){
        aux=2.0;
        return aux;
    }
    if(check_constraint_3(chrom,outputs,col,precision)!=0.0){
        aux=3.0;
        return aux;
    }
    if(check_constraint_4(chrom,inputs,outputs,col,precision)!=0.0){
        aux=4.0;
        return aux;
    }
    if(check_constraint_5(chrom)!=0.0){
        aux=5.0;
        return aux;
    }
    if(check_constraint_6(chrom)!=0.0){
        aux=6.0;
        return aux;
    }
    if(check_constraint_11(chrom)!=0.0){
        aux=11.0;
        return aux;
    }
    if(check_constraint_12(chrom)!=0.0){
        aux=12.0;
        return aux;
    }
    if(check_constraint_13(chrom)!=0.0){
        aux=13.0;
        return aux;
    }
    if(check_constraint_14(chrom)!=0.0){
        aux=14.0;
        return aux;
    }
    return aux;
}

/**
 * Ordena la lista de cromosomas de menor a mayor segun el valor de clasificación.
 * @param list Lista de cromosomas.
 * @param total Total de cromosomas.
 */

void classify_chromosome(Ptr_Chromosome *list, int total)
{
    int i,j;
	Ptr_Chromosome crom;
	double valor;

	for( i = 0 ; i < total ; i++ )
    {
		for( j = i+1 ; j < total ; j++ )
        {
			crom=list[i];
			valor=crom->evaluation;
			if ( ( valor > list[j]->evaluation  ) )
            {
				list[i]=list[j];
				list[j]=crom;
			}
		}
	}
	return;
}

void classify_chromosome_maximo(Ptr_Chromosome *list, double *inputs, int x, int col)
{
    int j,k;
	Ptr_Chromosome crom;
	double valor,valor2;

    for(k=0;k<x;k++){
                    for( j = k+1 ; j < x ; j++ )
                    {
                            crom=list[k];
                            valor=maximizar(crom, inputs, col);
                            valor2=maximizar(list[j],inputs, col);
                            if (valor < valor2)
                            {
                                list[k]=list[j];
                                list[j]=crom;
                            }
                        }
                    }
    return;
}
// Clasificacion por el valor de la pre-evaluacion. Cuando menor sea mejor candidato es el cromosoma.
/**
 * Ordena la lista de cromosomas de menor a mayor segun el valor de clasificación.
 * @param list Lista de cromosomas.
 * @param total Total de cromosomas.
 */
void pre_classify_chromosome(Ptr_Chromosome *list, int total)
{
	int i,j;
	Ptr_Chromosome crom;
	double valor;

	for( i = 0 ; i < total ; i++ )
    {
		for( j = i+1 ; j < total ; j++ )
        {
			crom=list[i];
			valor=crom->pre_evaluation;
			if ( ( valor > list[j]->pre_evaluation  ) )
            {
				list[i]=list[j];
				list[j]=crom;
			}
		}
	}
	return;
}


/**
 * Realiza un cruce entre dos crosmosomas. En cada cruce se generan 2 cromosomas descendientes.
 * @param a Cromosoma a. (padre)
 * @param b Cromosoma b. (madre)
 * @param c1 Cromosoma descendiente 1.
 * @param c2 Cromosoma descendiente 2.
 * @param *outputs, vector con las salidas
 * @param *inputs, vector con las entradas
 * @param col, variable auxiliar para fijar alguna columna
 */

void crossover(Ptr_Chromosome *lista, Ptr_Chromosome c1, double *inputs, double *outputs, int col,
               double *inputs_inv, double *outputs_inv, int cruzados)
{
    int i = 0, j = 0, fi, salir=0, elegido=0;
    double gamma, aux_tneg=0.0;

    elegido=rand()%cruzados;

    c1->metodo = lista[elegido]->metodo;
    c1->valido = lista[elegido]->valido;
    c1->fi= lista[elegido]->fi;
    c1->evaluation = lista[elegido]->evaluation;
    //show_chromosome(c1,inputs,outputs,col,inputs_inv,outputs_inv);
    //system("pause");

    c1->beta = lista[elegido]->beta;

    for(i=0;i<CANTIDAD_M;i++){
        c1->vik[i] = lista[elegido]->vik[i];
        c1->tneg[i] = lista[elegido]->tneg[i];
    }
    for(i=0;i<CANTIDAD_S;i++){
        c1->urk[i] = lista[elegido]->urk[i];
        c1->tpos[i] = lista[elegido]->tpos[i];
    }
    for(i=0;i<CANTIDAD_N;i++){
        c1->djk[i] = lista[elegido]->djk[i];
        c1->alfas[i] = lista[elegido]->alfas[i];
        c1->discrete_model[i] = lista[elegido]->discrete_model[i];
    }


    fi=(rand()%9)+1;
    //fi=(rand()%7)+3;

    switch(fi)
    {
        case 1:
            // Cruce BETA 1

            c1->fi=fi;
            salir=0;

            do{
            c1->beta=0.0;
            gamma = Random()/5.0;
            for(i=0;i<cruzados;i++){
            c1->beta+=lista[i]->beta;
            }
            c1->beta=c1->beta/cruzados;
            c1->beta+=gamma;
            salir++;
            }while (c1->beta>=1.0&&salir<50);

            if(c1->beta>=1.0){c1->beta=lista[rand()%cruzados]->beta;}

            recalcular(c1,inputs,outputs,col,inputs_inv,outputs_inv,fi);
            break;
        case 2:
            // Cruce BETA 2

            c1->fi=fi;
            salir=0;

            do{
            c1->beta=0.0;
            gamma = Random()/5.0;
            for(i=0;i<cruzados;i++){
            c1->beta+=lista[i]->beta;
            }
            c1->beta=c1->beta/cruzados;
            c1->beta-=gamma;
            salir++;
            }while (c1->beta>=1.0&&salir<50);

            if(c1->beta>=1.0){c1->beta=lista[rand()%cruzados]->beta;}

            recalcular(c1,inputs,outputs,col,inputs_inv,outputs_inv,fi);
            break;
        case 3:
            // Cruce TPOS-TNEG

            c1->fi=fi;
            for(j=0;j<CANTIDAD_S;j++){
                for(i=0;i<cruzados;i++){
                    aux_tneg+=lista[i]->tneg[j];
                }
                c1->tpos[j]=aux_tneg;
                aux_tneg=0.0;
            }
            recalcular(c1,inputs,outputs,col,inputs_inv,outputs_inv,fi);
            break;
        case 4:
            //CRUCE TPOS

            c1->fi=fi;
            for(i=0;i<CANTIDAD_S;i++){
                c1->tpos[i]=0.0;
            }

            for(j=0;j<CANTIDAD_S;j++){
                for(i=0;i<cruzados;i++){
                    c1->tpos[j]+=lista[i]->tpos[j];
                }
            }

            recalcular(c1,inputs,outputs,col,inputs_inv,outputs_inv,fi);
            break;
        case 5:
            // CRUCE VIK

             c1->fi=fi;
             for(i=0;i<CANTIDAD_S;i++){
                c1->vik[i]=0.0;
            }

            for(j=0;j<CANTIDAD_S;j++){
                for(i=0;i<cruzados;i++){
                    c1->vik[j]+=lista[i]->vik[j];
                }
            }

            recalcular(c1,inputs,outputs,col,inputs_inv,outputs_inv,fi);
            break;
        case 6:
            //CRUCE TPOS POR PUNTOS

            c1->fi=fi;

            for(i=0;i<CANTIDAD_S;i++){
                c1->tpos[i] = lista[rand()%cruzados]->tpos[i];
            }

            recalcular(c1,inputs,outputs,col,inputs_inv,outputs_inv,fi);
            break;
        case 7:
            //CRUCE TPOS-TNEG POR PUNTOS

            c1->fi=fi;

            for(i=0;i<CANTIDAD_S;i++){
                c1->tpos[i] = lista[rand()%cruzados]->tneg[i];
            }

            recalcular(c1,inputs,outputs,col,inputs_inv,outputs_inv,fi);
            break;
        case 8:
            //CRUCE TPOS-TNEG POR PUNTOS

            c1->fi=fi;


            for(i=0;i<CANTIDAD_S;i++){
                c1->tpos[i] = lista[rand()%cruzados]->tpos[i];
            }

            for(i=0;i<CANTIDAD_M;i++){
                c1->vik[i] = lista[rand()%cruzados]->vik[i];
            }
            recalcular(c1,inputs,outputs,col,inputs_inv,outputs_inv,fi);
            break;
        case 9:
            // Cruce TPOS-TNEG ARITMETICO

            c1->fi=fi;

            for(j=0;j<CANTIDAD_S;j++){
                aux_tneg=lista[0]->tneg[j];
                for(i=1;i<cruzados;i++){
                    if(aux_tneg<lista[i]->tneg[j]){
                        aux_tneg=lista[i]->tneg[j]-aux_tneg;
                    }else
                        {
                            aux_tneg-=lista[i]->tneg[j];
                        }
                }
                c1->tpos[j]=aux_tneg;
            }
            recalcular(c1,inputs,outputs,col,inputs_inv,outputs_inv,fi);
            break;
        default:
            break;
    }
    return;
}

/**
 * Realiza la mutación del cromosoma. Coge un término aleatorio y se modifica de forma aleatoria,
 * y se deducen los términos que dependan del término mutado
 * @param c Cromosoma a. (padre)
 * @param *outputs, vector con las salidas
 * @param *inputs, vector con las entradas
 * @param col, variable auxiliar para fijar alguna columna
 */
void mutate (Ptr_Chromosome c, double *inputs, double *outputs, int col, double *inputs_inv, double *outputs_inv)
{
	int i = 0, fi, j=0, aux=0;

    fi = (rand()%6)+1;

    //c->evaluation = BAD_CHROM;

    switch(fi){
        case 0:
            i = 0;
            j = CANTIDAD_N-1;
            // Mutación parte discreta. Parte de las b.
            while ( i < j )
            {
                c->discrete_model[i] = rand()%1;
                c->discrete_model[j] = rand()%1;
                // Al mismo tiempo intentamos que se cumplan las restricciones.
                // Hacer que alfa valga 0 y d 1 si b vale 1.

                if ( c->discrete_model[i] == 0 && aux<=(CANTIDAD_S+CANTIDAD_M-1) )
                {
                    c->alfas[i] = Random();
                    c->djk[i] = 0.0;
                    aux++;
                }else{
                    c->discrete_model[i]=1;
                    c->alfas[i] = 0.0;
                }
                if ( c->discrete_model[j] == 0 && aux<=(CANTIDAD_S+CANTIDAD_M-1)) {
                    c->alfas[j] = Random();
                    c->djk[j] = 0.0;
                    aux++;
                }else{
                    c->discrete_model[j]=1;
                    c->alfas[j] = 0.0;
                }
                i++;
                j--;
            }
            recalcular(c,inputs,outputs,col,inputs_inv,outputs_inv,fi);
            break;
        case 1:
            //Mutación de beta
            c->beta = Random();
            recalcular(c,inputs,outputs,col,inputs_inv,outputs_inv,fi);
            break;
        case 2:
            // Mutación tneg
            i = 0;
            j = CANTIDAD_M-1;
            if(j>i){
                while ( i < j )
                {
                    c->tneg[i] = (double)(rand()%1000);
                    c->tneg[j] = Random();
                    i++;
                    j--;
                }
            }else{
                c->tneg[i] = Random();
            }
            break;
        case 3:
            //Mutación tpos
            i = 0;
            j = CANTIDAD_S-1;
            if(j>i){
                while ( i < j )
                {
                    c->tpos[i] = (double)(rand()%1000);
                    c->tpos[j] = Random();
                    i++;
                    j--;
                }
            }else{
                c->tpos[i] = Random();
            }
            break;
        case 4:
            //Mutación de Vik
            i = 0;
            j = CANTIDAD_M-1;
            if(j>i){
                while ( i < j )
                {
                    c->vik[i] = (double)(rand()%1000);
                    c->vik[j] = (double)(rand()%1000);
                    i++;
                    j--;
                }
            }else{
                c->vik[i] = (double)(rand()%1000);
            }
            recalcular(c,inputs,outputs,col,inputs_inv,outputs_inv,fi);
            break;
        case 5:
            //Mutación de Urk
            i = 0;
            j = CANTIDAD_S-1;
            if(j>i){
                while ( i < j )
                {
                    c->urk[i] = (double)(rand()%1000);
                    c->urk[j] = (double)(rand()%1000);
                    i++;
                    j--;
                }
            }else{
                c->urk[i] = (double)(rand()%1000);
            }
            recalcular(c,inputs,outputs,col,inputs_inv,outputs_inv,fi);
            break;
        case 6:
            //Mutación de djk
            i = 0;
            j = CANTIDAD_N-1;
            while ( i < j )
            {
                if(c->discrete_model[i]==1){
                    c->djk[i] = (double)(rand()%1000);
                }
                if(c->discrete_model[j]==1){
                    c->djk[j] = (double)(rand()%1000);
                }
                i++;
                j--;
            }
            break;
        default:

            break;
    }



 	//printf("\n=============================== MUTACION ===================================");
    //printf("\nNum. cromosomas mutados: %d",j);
	//show_chromosome(c,inputs,outputs,col,inputs_inv,outputs_inv);


#ifdef VERBOSE
	printf("\n============================= NUEVO VALOR ==================================");
    show_chromosome(c);
#endif
	//printf("\n=============================== MUTACION ===================================");
	return;
}


////////////////////////////////////////////////////////////
//se toma aleatoriamente una variable y se cambia por otra

///////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////
double genetico(double *inputs, double *outputs, int col, FILE *fsolucion, Ptr_metaheuristica meta)
{

    int rep_mejor = 10,         // Repeticiones del mejor ?????
    fin = 0,                // ??????
    ale1 = 0,               // ??????
    iter = 0,
    validos = 0,
    frontera=0,
    no_vale=0,
    mut = 0,
    PEII = 0,
    i = 0,
    j=0,
    k = 0,
    auxiliar=0,
    auxiliarW=0,
    iguales=0,
    cruzados=0,
    NIRE = 0,
    abortar=0,
    capacidad=0,
    no_validos=0,
    *indices;
    double porcentaje, *Vx_columna, *Vy_columna, *inputs_inv, *outputs_inv, tiempo, maximo=0.0, *sumatorio_salidas, *sumatorio_entradas;
    double objetivo=0.0, objetivo2=0.0;
    //struct timeval ti, tf;

    //double *Lista_Maximos;


    Ptr_Chromosome *cromosoma_cruce=(Ptr_Chromosome*)malloc(meta->CANTIDAD_CRUCE*sizeof(Ptr_Chromosome));
    Ptr_Chromosome *cromosoma_cruce_malos=(Ptr_Chromosome*)malloc(meta->CANTIDAD_CRUCE*sizeof(Ptr_Chromosome));
    indices = (int*)malloc(meta->INEIni*sizeof(int));

    Vx_columna = (double*)malloc(CANTIDAD_N*sizeof(double));
    Vy_columna = (double*)malloc(CANTIDAD_N*sizeof(double));

    inputs_inv= (double*)malloc(CANTIDAD_M*CANTIDAD_N*sizeof(double));
    outputs_inv= (double*)malloc(CANTIDAD_S*CANTIDAD_N*sizeof(double));

    sumatorio_salidas = (double*)malloc(CANTIDAD_N*sizeof(double));
    sumatorio_entradas = (double*)malloc(CANTIDAD_N*sizeof(double));

    //Lista_Maximos = (double*)malloc(INEIni*sizeof(double));

    inversos(inputs,outputs,inputs_inv,outputs_inv);
    suma_columnas2(inputs, outputs, sumatorio_entradas, sumatorio_salidas);
/*
    //Para medir el tiempo
	struct timeval *tv;
	struct timezone *tz;
	long tf,ti,si,sf;


	tv=(struct timeval *)malloc(sizeof(struct timeval));
	tz=NULL;

    gettimeofday(tv,tz);
    si=(tv->tv_sec);
    ti=(tv->tv_usec);
    ti=si*1000000+ti;
*/
    //Reservar memoria para todos los cromosomas
    //Ptr_Chromosome *List_Chromosome=(Ptr_Chromosome*)malloc((INEIni)*sizeof(Ptr_Chromosome));
    Ptr_Chromosome *List_Chromosome = NULL;
    if(meta->INEIni>(meta->FNEIni+meta->PBBCom+meta->PWWCom))
    {
        capacidad=meta->INEIni;
        List_Chromosome=(Ptr_Chromosome*)malloc(capacidad*sizeof(Ptr_Chromosome));
    }else{
        capacidad=meta->FNEIni+meta->PBBCom+meta->PWWCom;
        List_Chromosome=(Ptr_Chromosome*)malloc(capacidad*sizeof(Ptr_Chromosome));
    }
    suma_columnas(inputs,outputs,Vx_columna,Vy_columna);
    //Creamos todos los cromosomas
    //gettimeofday(&ti, NULL);
    frontera=0;
    no_vale=0;


    printf("\n");
	for ( i = 0 ; i < meta->INEIni ; i++ )
	{
	    printf("%d",i);
		List_Chromosome[i]=create_chromosome(i,inputs,outputs,col,Vx_columna,Vy_columna, inputs_inv, outputs_inv, sumatorio_entradas, sumatorio_salidas, 1);

		fflush(stdout);
		if(List_Chromosome[i]->valido==1)
            {
                validos+=1;
            }
        printf("F = %f\n",maximizar(List_Chromosome[i],inputs,col));
        fflush(stdout);

		if(maximizar(List_Chromosome[i],inputs,col)>1.0){
            frontera++;
            no_vale++;
		}else{frontera=0;}
		if((frontera>(meta->INEIni/10)||no_vale>(meta->INEIni/3))&&validos==0){
            fin=1;
            //Lista_Maximos[0]=0.0;
		}
		/*if(validos==0)
            {
            if(List_Chromosome[i]->terminar==1)
                {
                    abortar++;
                }else {abortar=0;}

            if(abortar==INEIni/3){fin=1;}
            }*/
		if(fin==1){i=meta->INEIni;}

	}

	frontera=0;
	no_vale=0;
//    gettimeofday(&tf, NULL);   // Instante final
    //tiempo= (tf.tv_sec - ti.tv_sec)*1000 + (tf.tv_usec - ti.tv_usec)/1000.0;
    //printf("Has tardado: %g milisegundos\n", tiempo);
	porcentaje=((double)validos*100.0)/(double)meta->INEIni;
	printf("\nPorcentaje antes: %f\n",porcentaje);
	//system("pause");

	rep_mejor=0;
    //if(validos==0){fin=1;}


    fflush(stdout);

	if(fin==0){
        for(i=0;i<meta->INEIni;i++){
            if(List_Chromosome[i]->valido==1){
                if(maximo<maximizar(List_Chromosome[i],inputs,col)) {
                    maximo=maximizar(List_Chromosome[i],inputs,col);
                }
            }
        }
        printf("\nMAXIMO PRIMERA GENERACION ANTES: %f\n",maximo);
        maximo=0.0;

            for (i=0 ;i<meta->INEIni;i++)
            {
                List_Chromosome[i]->evaluation=evaluate_chromosome(List_Chromosome[i],inputs,outputs,col,inputs_inv,outputs_inv);
            }

        /*for(i=0;i<capacidad;i++){
            if(List_Chromosome[i]->valido==1){
                printf("\nList_Chromosome[%d] || ETIQ = %d || F = %f || Valido = %d || EVAL = %f ",i,List_Chromosome[i]->etiq,maximizar(List_Chromosome[i],inputs,col),List_Chromosome[i]->valido,List_Chromosome[i]->evaluation);
                show_chromosome(List_Chromosome[i],inputs,outputs,col,inputs_inv,outputs_inv);
            }
        }*/

    ///////////////////////////////////////////////////////////////////////MEJORA///////////////////////////////////////////////////////////////////////////

        printf("\n-------------------MEJORA---------------------\n");
            for(i=0;i<meta->INEIni;i++){
                if(List_Chromosome[i]->valido==0){
                    PEII=(rand()%(100))+1;
                    if ( PEII<=meta->PEIIni ){
                    printf("%d ",i);
                    mejora(List_Chromosome[i],inputs,outputs,col,inputs_inv,outputs_inv, meta->IIEIni);
                    List_Chromosome[i]->valido=good_chromosome(List_Chromosome[i],inputs,outputs,col,inputs_inv,outputs_inv);
                    fflush(stdout);
                    }
                }
            }
        printf("\n----------------------------------------------\n");

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        fflush(stdout);
        validos=0;
        for (i=0;i<meta->INEIni;i++)
        {
            if(List_Chromosome[i]->valido==1){validos+=1;}
        }

        porcentaje=((double)validos*100.0)/(double)meta->INEIni;
        printf("\nPorcentaje despues: %f\n",porcentaje);


        for ( i = 0  ; i < meta->INEIni ; i++ ){
            List_Chromosome[i]->evaluation=evaluate_chromosome(List_Chromosome[i],inputs,outputs,col,inputs_inv,outputs_inv);
		}
		classify_chromosome(List_Chromosome,meta->INEIni);
        classify_chromosome_maximo(List_Chromosome, inputs, validos,col);


        fflush(stdout);

        for(i=0;i<meta->INEIni;i++){
            if(List_Chromosome[i]->valido==1){
                if(maximo<maximizar(List_Chromosome[i],inputs,col)) {
                    maximo=maximizar(List_Chromosome[i],inputs,col);
                }
            }
        }

        printf("\nMAXIMO PRIMERA GENERACION DESPUES: %f\n",maximo);

        for(i=meta->INEIni;i<capacidad;i++){
                printf("%d",i);
                //List_Chromosome[i]=cromosoma_vacio(i);
                List_Chromosome[i]=create_chromosome(i,inputs,outputs,col,Vx_columna,Vy_columna, inputs_inv, outputs_inv, sumatorio_entradas, sumatorio_salidas, 0);
                //memcpy(List_Chromosome[i], List_Chromosome[rand()%FNEIni], sizeof(struct Chromosome));
                printf("\n");
                fflush(stdout);
        }

    }
    fflush(stdout);
    printf("!!!!!!!!!!!!! VALIDOS = %d !!!!!!!!!!!!!\n",validos);

    /*for(i=0;i<INEIni;i++){
            printf("\nCromosoma[%d] = %f || V = %d || E = %f",i,maximizar(List_Chromosome[i],inputs,col),List_Chromosome[i]->valido, List_Chromosome[i]->evaluation);
        }
        system("pause");*/

    if(validos==0){fin=1;}
    if (fin==0){/*fprintf(fsolucion, "\n\nEmpresa %d: %f |||| Validos = %f",col+1, maximo, porcentaje);*/}



	while( fin == 0 )
    {
        if(meta->PBBCom!=0 || meta->PBWCom!=0 || meta->PWWCom!=0){
        no_validos=0;
        //gettimeofday(&ti, NULL);
		//EVALUACION Y CLASIFICACION DE LOS CROMOSOMAS VALIDOS SI LOS HUBIESE.

		for ( i = 0  ; i < meta->FNEIni ; i++ ){
            List_Chromosome[i]->evaluation=evaluate_chromosome(List_Chromosome[i],inputs,outputs,col,inputs_inv,outputs_inv);
		}
		classify_chromosome(List_Chromosome,meta->FNEIni);
        validos=0;
		for(i=0;i<meta->FNEIni;i++){
            List_Chromosome[i]->valido=good_chromosome(List_Chromosome[i],inputs,outputs,col,inputs_inv,outputs_inv);
            if(List_Chromosome[i]->valido==1) validos++;
        }
        classify_chromosome_maximo(List_Chromosome, inputs, validos,col);

        /*for(i=0;i<FNEIni;i++){
            printf("\nCromosoma[%d] = %f || V = %d || E = %f",i,maximizar(List_Chromosome[i],inputs,col),List_Chromosome[i]->valido, List_Chromosome[i]->evaluation);
        }
        system("pause");*/

        /////////////////////////////////////////// SELECCION BUENOS /////////////////////////////////////////////

        if(meta->NBESel>validos){
            auxiliar=validos;
            if(validos<2){no_validos=1;}
        }else{
            auxiliar=meta->NBESel;
             }

        /////////////////////////////////////////// SELECCION MALOS /////////////////////////////////////////////

        if(meta->NWESel>meta->FNEIni-validos){
            auxiliarW=meta->FNEIni-validos;
            if(auxiliarW<=3){auxiliarW=0;}
        }else{
            auxiliarW=meta->NWESel;
             }

		////////////////////////////////////////////////////////////////// CRUCE ENTRE BUENOS //////////////////////////////////////////////////////////////////////////////
		if(no_validos==0){
		for ( i = meta->FNEIni ; i < meta->FNEIni+meta->PBBCom ; i++)
        {
            cruzados=0;
            if(meta->CANTIDAD_CRUCE>=auxiliar){
                for(j=0;j<auxiliar;j++)
                    {
                    cromosoma_cruce[j]=List_Chromosome[j];
                    cruzados++;
                    }
            }else
            {
                for(j=0;j<meta->CANTIDAD_CRUCE;j++){
                    do{
                        iguales=0;
                        indices[j]=(rand()%(auxiliar))+1;//solo se reproducen los primeros. Coge 2 aleatorios entre la primera mitad
                        for(k=0;k<j;k++){
                            if(indices[k]==indices[j]){
                               iguales=1;
                            }
                        }
                    }while(iguales==1);
                    cromosoma_cruce[j]=List_Chromosome[indices[j]];
                    cruzados++;
                }
            }

            crossover(cromosoma_cruce, List_Chromosome[i], inputs, outputs, col, inputs_inv, outputs_inv,cruzados);
            PEII=(rand()%(101));
            if ( PEII<=meta->PEIImp ){
            mejora(List_Chromosome[i],inputs,outputs,col,inputs_inv,outputs_inv, meta->IIEImp);
            fflush(stdout);
            }
		}

		}

		////////////////////////////////////////////////////////////////// CRUCE ENTRE MALOS //////////////////////////////////////////////////////////////////////////////
		/*printf("\nAUXILIAR_MALOS = %d\n",auxiliarW);
		system("pause");*/
		if(auxiliarW!=0){
            for ( i = meta->FNEIni + meta->PBBCom ; i < meta->FNEIni+meta->PBBCom+meta->PWWCom ; i++)
            {
                cruzados=0;
                if(meta->CANTIDAD_CRUCE>=auxiliarW){
                    for(j=validos;j<validos+auxiliarW;j++)
                        {
                        cromosoma_cruce_malos[j]=List_Chromosome[j];
                        cruzados++;
                        }
                }else
                {
                    for(j=0;j<meta->CANTIDAD_CRUCE;j++){
                        do{
                            iguales=0;
                            indices[j]=(rand()%(meta->FNEIni-validos))+validos;//solo se reproducen los primeros. Coge 2 aleatorios entre la primera mitad
                            for(k=0;k<j;k++){
                                if(indices[k]==indices[j]){
                                   iguales=1;
                                }
                            }
                        }while(iguales==1);
                        cromosoma_cruce_malos[j]=List_Chromosome[indices[j]];
                        cruzados++;
                    }
                }

                /*for(i=0;i<CANTIDAD_CRUCE;i++){
                    printf("\n Cromosoma cruzado ETIQUETA = %d y VALIDO = %d",cromosoma_cruce_malos[i]->etiq, cromosoma_cruce_malos[i]->valido);
                }
                system("pause");*/
                crossover(cromosoma_cruce_malos, List_Chromosome[i], inputs, outputs, col, inputs_inv, outputs_inv,cruzados);
                PEII=(rand()%(101));
                if ( PEII<=meta->PEIImp ){
                mejora(List_Chromosome[i],inputs,outputs,col,inputs_inv,outputs_inv, meta->IIEImp);
                fflush(stdout);
                }
            }
		}

		////////////////////////////////////////////////////////////////// CRUCE ENTRE BUENOS Y MALOS //////////////////////////////////////////////////////////////////////////////

		/*for ( i = FNEIni+PBBCom+PWWCom ; i < capacidad ; i++)
        {
            cruzados=0;
            if(CANTIDAD_CRUCE>=auxiliar){
                for(j=0;j<auxiliar;j++)
                    {
                    cromosoma_cruce[j]=List_Chromosome[j];
                    cruzados++;
                    }
            }else
            {
                for(j=0;j<CANTIDAD_CRUCE;j++){
                    do{
                        iguales=0;
                        indices[j]=(rand()%(auxiliar))+1;//solo se reproducen los primeros. Coge 2 aleatorios entre la primera mitad
                        for(k=0;k<j;k++){
                            if(indices[k]==indices[j]){
                               iguales=1;
                            }
                        }
                    }while(iguales==1);
                    cromosoma_cruce[j]=List_Chromosome[indices[j]];
                    cruzados++;
                }
            }

            crossover(cromosoma_cruce, List_Chromosome[i], inputs, outputs, col, inputs_inv, outputs_inv,cruzados);
            PEII=(rand()%(101));
            if ( PEII<=PEIImp ){
            mejora(List_Chromosome[i],inputs,outputs,col,inputs_inv,outputs_inv, IIEImp);
            }
		}*/


		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		//////////////////////////////////////////////////////////EVALUACION Y CLASIFICACION/////////////////////////////////////////////////////////////////////

		for ( i = 0  ; i < capacidad ; i++ )
        {
            List_Chromosome[i]->evaluation=evaluate_chromosome(List_Chromosome[i],inputs,outputs,col,inputs_inv,outputs_inv);
		}
		classify_chromosome(List_Chromosome,capacidad);
        validos=0;
		for(i=0;i<capacidad;i++){
            List_Chromosome[i]->valido=good_chromosome(List_Chromosome[i],inputs,outputs,col,inputs_inv,outputs_inv);
            if(List_Chromosome[i]->valido==1) validos++;
        }
        classify_chromosome_maximo(List_Chromosome, inputs, validos,col);

		/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        /*for(i=0;i<PBBCom;i++){
                printf("\nCromosoma[%d] = %f",i,maximizar(List_Chromosome[i],inputs,col));
        }
        system("pause");*/
		//MUTACION

		mut=(rand()%(100))+1;
		if ( mut<=meta->PEDImp ){	//el 10% de que haya mutacion
            printf("\n------------------------MUTADO----------------------------");
			ale1=(rand()%capacidad/2)+capacidad/2;
			mutate(List_Chromosome[ale1],inputs,outputs,col,inputs_inv,outputs_inv);
			mejora(List_Chromosome[ale1],inputs,outputs,col,inputs_inv,outputs_inv, meta->IIDImp);
			printf("\n");
			fflush(stdout);
			 //show_chromosome(List_Chromosome[ale1],inputs,outputs,col,inputs_inv,outputs_inv);
			 //system("pause");
		}


		//////////////////////////////////////////////////////////EVALUACION Y CLASIFICACION/////////////////////////////////////////////////////////////////////

		for ( i = 0  ; i < capacidad ; i++ )
        {
            List_Chromosome[i]->evaluation=evaluate_chromosome(List_Chromosome[i],inputs,outputs,col,inputs_inv,outputs_inv);
		}
		classify_chromosome(List_Chromosome,capacidad);
        validos=0;
		for(i=0;i<capacidad;i++){
            List_Chromosome[i]->valido=good_chromosome(List_Chromosome[i],inputs,outputs,col,inputs_inv,outputs_inv);
            if(List_Chromosome[i]->valido==1) validos++;
        }
        classify_chromosome_maximo(List_Chromosome, inputs, validos,col);

        /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

        objetivo2=maximizar(List_Chromosome[0],inputs,col);
        if(objetivo2!=objetivo){
		//fprintf(fsolucion,"\n%d ------> F = %f", iter, maximizar(List_Chromosome[0],inputs,col));
		objetivo=maximizar(List_Chromosome[0],inputs,col);
        }

        iter++;

		/*for(i=0;i<capacidad;i++){
                printf("\nCromosoma[%d] = %f || VALIDO = %d || EVAL = %f",i,maximizar(List_Chromosome[i],inputs,col), List_Chromosome[i]->valido, List_Chromosome[i]->evaluation);
            }
        system("pause");*/

		if(maximo==maximizar(List_Chromosome[0],inputs,col)){
            NIRE++;
            if(NIRE==meta->NIREnd){
                printf("\n\t\t   -----------NIREnd------------- (ITER = %d)",iter);
                iter=meta->MNIEnd;

            }
		}else{
		    maximo=maximizar(List_Chromosome[0],inputs,col);
		    NIRE=0;
		    }

        }else{iter=meta->MNIEnd;}

        if(iter==meta->MNIEnd){

            maximo=maximizar(List_Chromosome[0],inputs,col);
            printf("\nMAXIMO primera parada: %f\n",maximo);

            validos=0;

            for(i=0;i<meta->FNEIni;i++){
                List_Chromosome[i]->valido=good_chromosome(List_Chromosome[i],inputs,outputs,col,inputs_inv,outputs_inv);
                if(List_Chromosome[i]->valido==1) validos++;
            }

            porcentaje=((double)validos*100.0)/(double)meta->FNEIni;
            printf("\n\nValidos GENETICO ANTES MEJORA:%f\n",porcentaje);
            fflush(stdout);

            maximo=maximizar(List_Chromosome[0],inputs,col);
            printf("\nMAXIMO segunda parada: %f\n",maximo);

            /*for(i=0;i<capacidad;i++){
                if(List_Chromosome[i]->valido==1){
                    printf("\nList_Chromosome[%d] || ETIQ = %d || F = %f || Valido = %d || EVAL = %f ",i,List_Chromosome[i]->etiq,maximizar(List_Chromosome[i],inputs,col),List_Chromosome[i]->valido,List_Chromosome[i]->evaluation);
                }
            }*/

            /*------ MEJORO TODOS LOS CROMOSOMAS FINALES NO VALIDOS ----------*/
            printf("\n-------------------MEJORA---------------------\n");
            fflush(stdout);
            for(i=0;i<meta->FNEIni;i++){
                if(List_Chromosome[i]->valido==0){
                    PEII=(rand()%(100))+1;
                    if ( PEII<=meta->PEIImp ){
                    printf("%d ",i);
                    fflush(stdout);
                    mejora(List_Chromosome[i],inputs,outputs,col,inputs_inv,outputs_inv, meta->IIEImp);
                    List_Chromosome[i]->valido=good_chromosome(List_Chromosome[i],inputs,outputs,col,inputs_inv,outputs_inv);
                    }
                }
            }
            printf("\n----------------------------------------------\n");
            fflush(stdout);

            /* ------- VUELVO A MIRAR CUANTOS VALIDOS HAY DESPUES DE MEJORAR ---------*/
            validos=0;
            for(i=0;i<meta->FNEIni;i++){
                List_Chromosome[i]->valido=(good_chromosome(List_Chromosome[i],inputs,outputs,col,inputs_inv,outputs_inv));
                if(List_Chromosome[i]->valido==1){
                    validos++;
                }
            }

            maximo=maximizar(List_Chromosome[0],inputs,col);
            printf("\nMAXIMO tercera parada: %f\n",maximo);

            porcentaje=((double)validos*100.0)/(double)meta->FNEIni;
            printf("\n\nValidos GENETICO DESPUES MEJORA:%f\n",porcentaje);
            fflush(stdout);

            /*for(i=0;i<capacidad;i++){
                if(List_Chromosome[i]->valido==1){
                    printf("\nList_Chromosome[%d] || ETIQ = %d || F = %f || Valido = %d || EVAL = %f ",i,List_Chromosome[i]->etiq,maximizar(List_Chromosome[i],inputs,col),List_Chromosome[i]->valido,List_Chromosome[i]->evaluation);
                }
            }*/

            /*--------- CON LOS NUEVOS CROMOSOMAS, EVALUO Y CLASIFICO POR SI HAY ALGUN VALIDO NUEVO MEJOR --------*/

            for ( i = 0  ; i < meta->FNEIni ; i++ )
            {
                List_Chromosome[i]->evaluation=evaluate_chromosome(List_Chromosome[i],inputs,outputs,col,inputs_inv,outputs_inv);
            }
            classify_chromosome(List_Chromosome,meta->FNEIni);
            classify_chromosome_maximo(List_Chromosome, inputs, validos,col);

            /* ------- CREO UN VECTOR CON TODOS LOS MAXIMOS DE LOS CROMOSOMAS VALIDOS -----------*/
            k=0;
            /*for(i=0;i<validos;i++){
                    if(List_Chromosome[i]->evaluation==0.000000){
                            Lista_Maximos[k]=maximizar(List_Chromosome[i],inputs,col);
                            //printf("\nValor_Funcion[%d] = %f con BETA = %f |||VALIDO = %d|||",i,Lista_Maximos[k], List_Chromosome[i]->beta, List_Chromosome[i]->valido);
                            k++;
                    }
            }*/
            /*for(i=0;i<INEIni;i++){
                printf("\nCromosoma[%d] = %f || VALIDO = %d || EVAL = %f",i,maximizar(List_Chromosome[i],inputs,col), List_Chromosome[i]->valido, List_Chromosome[i]->evaluation);
            }*/
            //gettimeofday(&tf, NULL);   // Instante final
            //tiempo= (tf.tv_sec - ti.tv_sec)*1000 + (tf.tv_usec - ti.tv_usec)/1000.0;
            //printf("\n\nHas tardado: %g milisegundos\n", tiempo);
            //fprintf(fsolucion,"%f %d\n",Lista_Maximos[0]-maximo, List_Chromosome[0]->fi);
            //fprintf(fsolucion, "\t|| %f\t|| Metodo = %d\n", Lista_Maximos[0],List_Chromosome[0]->fi);
            //printf("\nDIFERENCIA = %f con metodo de cruce FI = %d\n\n",Lista_Maximos[0]-maximo, List_Chromosome[0]->fi);
            show_chromosome(List_Chromosome[0],inputs,outputs,col,inputs_inv,outputs_inv);
            fprintf(fsolucion,"\nEmpresa %d: %f con %f%% de validos || Fi = %d",col+1,maximizar(List_Chromosome[0],inputs,col),porcentaje, List_Chromosome[0]->fi);
            //fprintf(fsolucion,"\nEmpresa %d: %f con %f%% de validos",col+1,maximizar(List_Chromosome[0],inputs,col),porcentaje);
            maximo=maximizar(List_Chromosome[0],inputs,col);
            //system("pause");
            //free(List_Chromosome);
            /*for(i=0;i<capacidad;i++){
                free_chromosome(List_Chromosome[i]);
            }*/
            fflush(stdout);
            free(Vy_columna);
            free(Vx_columna);
            free(sumatorio_entradas);
            free(sumatorio_salidas);
            free(inputs_inv);
            free(outputs_inv);
            free(indices);
            free(List_Chromosome);
            free(cromosoma_cruce);
            free(cromosoma_cruce_malos);
            return (maximo);
            iter=0;
        }

		fflush(stdout);
	}

	free(List_Chromosome);
	free(cromosoma_cruce);
    free(cromosoma_cruce_malos);
	free(sumatorio_entradas);
    free(sumatorio_salidas);
	//free(Lista_Maximos);
    free(Vy_columna);
    free(Vx_columna);
    free(inputs_inv);
    free(outputs_inv);
    free(indices);
    //fprintf(fsolucion,"%f\n",Lista_Maximos[0]-maximo);
    //printf("\n DIFERENCIA = %f",Lista_Maximos[0]-maximo);
    return (0.0);
    printf("\n\nFIN\nrep_mejor: %d iter: %d",rep_mejor, iter);
}


/**
 * @brief Función para calcular las t- mediante la fórmula de la restricción 1.2
 * @param Chrom , es el cromosoma
 * @param *inputs, vector con las salidas
 * @param col, variable auxiliar para fijar alguna columna
 */
void modify_constraint_1 ( Ptr_Chromosome chrom, double *inputs, int col)
{
    int i,k;
    double *auxiliar;

    auxiliar = (double*)malloc(CANTIDAD_M*sizeof(double));

    for(i=0;i<CANTIDAD_M;i++){
            auxiliar[i]=0.0;
        for(k=0;k<CANTIDAD_N;k++){
            auxiliar[i]+=chrom->alfas[k]*inputs[(i*CANTIDAD_N)+k];
        }
    }

    for(i=0;i<CANTIDAD_M;i++){
        auxiliar[i]*=(-1.0);
    }

    for(i=0;i<CANTIDAD_M;i++){
        auxiliar[i]+=chrom->beta*inputs[(i*CANTIDAD_N)+col];
        chrom->tneg[i]=auxiliar[i];
    }

    free(auxiliar);
    return;
}


/**
 * @brief Función para calcular las t+ mediante la fórmula de la restricción 1.3
 * @param Chrom , es el cromosoma
 * @param *outputs, vector con las salidas
 * @param col, variable auxiliar para fijar alguna columna
 */
void modify_constraint_2 ( Ptr_Chromosome chrom, double *outputs, int col )
{
    int i,k;
    double *aux;

    aux = (double*)malloc(CANTIDAD_S*sizeof(double));

    for(i=0;i<CANTIDAD_S;i++){
            aux[i]=0.0;
        for(k=0;k<CANTIDAD_N;k++){
            aux[i]+=chrom->alfas[k]*outputs[(i*CANTIDAD_N)+k];
        }
    }
    for(i=0;i<CANTIDAD_S;i++){
        aux[i]-=chrom->beta*outputs[(i*CANTIDAD_N)+col];
        chrom->tpos[i]=aux[i];
    }

    free(aux);
    return;
}


/**
 * @brief Función para deducir los valores alfa para satisfacer las restricciones 1.10, 1.2 y 1.3
 * @param Chrom , es el cromosoma
 * @param *outputs, vector con las salidas
 * @param col, variable auxiliar para fijar alguna columna
 * @param alfamenor, indice del vector que contiene las alfas el cual apunta al alfa mas eficiente de modificar para satisfacer la restriccion 1
 */

void deducir_valores_alfa(Ptr_Chromosome chrom, double *outputs, int col, int alfamenor, double *outputs_inv, double *Vx_columna, double *Vy_columna)
{
    int valido=0,j=0,k,i;
    double sumatorio=0.0, m=99999.99;
    double *vectorP;

    vectorP = (double*)malloc(CANTIDAD_N*sizeof(double));

    if(vectorP==NULL){
        printf("Fallo en malloc vectorP");
    }

    //Hasta conseguir que todos los alfas sean mayor que 0 y se cumpla la restricción 1
    while(valido==0)
    {
        sumatorio=0.0;
        //Calculamos el vector Yrj/Yrk
        for(k=0;k<CANTIDAD_N;k++){
            vectorP[k]=0.0;
            for(j=0;j<CANTIDAD_S;j++){
                vectorP[k]+= outputs[(j*CANTIDAD_N)+k]*outputs_inv[(j*CANTIDAD_N)+col];
            }
        }

        //Calculamos la suma de todos los alfas menos el de la posición menor por los sumatorios
        for(i=0;i<CANTIDAD_N;i++){
            sumatorio+=vectorP[i]*chrom->alfas[i];
        }
        sumatorio-=vectorP[alfamenor]*chrom->alfas[alfamenor];

        //Calculamos que debe valer el alfa menor para cumplir la restricción 1
        chrom->alfas[alfamenor]=(((double)CANTIDAD_S-sumatorio)/vectorP[alfamenor]);

        //Disminuimos el valor de todos los alfas para que el alfa menor sea positivo
        if(chrom->alfas[alfamenor]<0.0)
        {
            for(i=0;i<CANTIDAD_N;i++){
                chrom->alfas[i]*=0.95;
            }
            for(i=0;i<CANTIDAD_N;i++){
                if(m>((Vy_columna[i]+Vx_columna[i]))&&(chrom->alfas[i]!=0.0)){
                    m=Vy_columna[i]+Vx_columna[i];
                    j=i;
                }
            }
            chrom->alfas[j]*=0.0;
        }else{
            valido=1;
        }
    }
    free(vectorP);
    return;
}

/**
 * @brief Función para deducir los valores u y v para satisfacer las restricciones 1.4, 1.5, 1.6 y 1.13
 * @param Chrom , es el cromosoma
 * @param *inputs, vector con las entradas
 * @param *outputs, vector con las salidas
 * @param col, variable auxiliar para fijar alguna columna
 */
void deducir_valores_uv(Ptr_Chromosome chrom, double *inputs, double *outputs, int col)
{
    double aux = 0.0,
           aux2 = 0.0,
           *matriz,
           *vector,
           val = 1.0,
           precision = 0.000001;
    int    i = 0,
           j = 0,
           k = 0,
           p = 0,
           u = 0,
           dceros = 0,
           indice = 0,
           aleatorio = 0,
           maxrepeticiones=0;
    int info, ldb, lda, nrhs, n;
    int *ipiv;

    for(i=0;i<CANTIDAD_N;i++)
    {
        if(chrom->discrete_model[i]==0)dceros++;
    }

    if(dceros==1){
        for(i=0;i<CANTIDAD_N;i++)
            {
                if(chrom->discrete_model[i]==0)
                {
                indice=i;
                }
            }
    }
    while(val!=0.0&&maxrepeticiones<50)
    {
        for ( i = 0 ; i < CANTIDAD_M ; i++ )
            chrom->vik[i]=(double)(rand()%10000)+Random();
        for ( i = 0 ; i < CANTIDAD_S ; i++ )
            chrom->urk[i]=(double)(rand()%10000)+Random();

        if(dceros==1)
        {
            aleatorio=rand()%CANTIDAD_S;
            aux=0;
            for(i=0;i<CANTIDAD_M;i++)
            {
                aux+=chrom->vik[i]*inputs[(i*CANTIDAD_N)+indice];

            }
            aux2=0;
                for(i=0;i<CANTIDAD_S;i++)
            {
                if(i!=aleatorio)aux2+=chrom->urk[i]*outputs[(i*CANTIDAD_N)+indice];
            }
            chrom->urk[aleatorio]=(aux-aux2)/outputs[(aleatorio*CANTIDAD_N)+indice];
            calculo_D(chrom,inputs,outputs);
            val=0.0;
            val+=check_constraint_4(chrom, inputs, outputs, col , precision);
            val+=check_constraint_13(chrom);
            val+=check_constraint_5(chrom);
            val+=check_constraint_6(chrom);


        //En caso de que tengamos mas ecuaciones utilizaremos la descomposicion LU mediante
        //funciones LAPACK
        }else{
            fflush(stdout);
            matriz = (double*)malloc(dceros*dceros*sizeof(double));
            vector = (double*) malloc (dceros*dceros*sizeof(double));
            ipiv = (int*) malloc (dceros*sizeof(int));
            if(matriz==NULL){
                printf("Fallo en malloc matriz\n");
            }
            if(vector==NULL){
                printf("Fallo en malloc vector\n");
            }
            fflush(stdout);
            k=0;
            for(j=0;j<CANTIDAD_N;j++)
            {
                if(chrom->discrete_model[j]==0)
                {
                    aux=0;
                    for(i=0;i<(CANTIDAD_M+CANTIDAD_S-dceros);i++)
                    {
                        if(i<CANTIDAD_M)
                        {
                            aux+=chrom->vik[i]*inputs[(i*CANTIDAD_N)+j];
                        }else{
                            aux-=chrom->urk[i-CANTIDAD_M]*outputs[((i-CANTIDAD_M)*CANTIDAD_N)+j];

                        }
                    }
                    vector[k]=aux;
                    p=0;
                    u=0;
                    for(i=0;i<dceros;i++)
                    {
                        if((dceros-i)>CANTIDAD_S)
                        {
                            matriz[(u*dceros)+k]=0-inputs[((CANTIDAD_M-(dceros-CANTIDAD_S-i))*CANTIDAD_N)+j];
                        }else{
                            matriz[(u*dceros)+k]=outputs[(p*CANTIDAD_N)+j];
                            p++;
                        }
                        u++;
                    }
                    k++;
                }
            }
            n=dceros;
            ldb=dceros;
            lda=dceros;
            nrhs=dceros;
            dgesv_( &n, &nrhs, matriz, &lda, ipiv, vector,&ldb, &info);
            k=0;
            for(i=0;i<dceros;i++)
            {
                if((dceros-i)>CANTIDAD_S)
                {
                    chrom->vik[CANTIDAD_M-(dceros-CANTIDAD_S-i)]=vector[i];
                }else{
                    chrom->urk[k]=vector[i];
                    k++;
                }
            }
           /* for(i=0;i<CANTIDAD_S;i++){
                if(chrom->urk[i]<0){
                    mejorar_u(chrom,inputs,outputs);
                }
            }*/
            //Comprobar t- si hay alguna negativa y en ese caso recalcular valores.

            calculo_D(chrom,inputs,outputs);
            for(i=0;i<CANTIDAD_N;i++){
                if(chrom->djk[i]<0.000001 && chrom->djk[i]>-0.000001) chrom->djk[i]=0.0;
            }
            val=0.0;
            val+=check_constraint_4(chrom, inputs, outputs, col , precision);
            val+=check_constraint_13(chrom);
            val+=check_constraint_5(chrom);
            val+=check_constraint_6(chrom);
            //printf("\n\nVALIDO = %f",val);

            free(matriz);
            free(vector);
            free(ipiv);
        }
        maxrepeticiones++;
        //system("pause");
    }
    //if(maxrepeticiones==50){chrom->contador+=1;}
    return;
}

void deducir_valores_uv2(Ptr_Chromosome chrom, double *inputs, double *outputs, int col)
{
    double aux = 0.0,
           aux2 = 0.0,
           *matriz,
           *vector,
           val = 1.0,
           precision = 0.000001;
    int    i = 0,
           j = 0,
           k = 0,
           p = 0,
           u = 0,
           dceros = 0,
           indice = 0,
           aleatorio = 0,
           maxrepeticiones=0;
    int info, ldb, lda, nrhs, n;
    int *ipiv;

    for(i=0;i<CANTIDAD_N;i++)
    {
        if(chrom->discrete_model[i]==0)dceros++;
    }

    if(dceros==1){
        for(i=0;i<CANTIDAD_N;i++)
            {
                if(chrom->discrete_model[i]==0)
                {
                indice=i;
                }
            }
    }
    while(val!=0.0&&maxrepeticiones<50)
    {

        if(dceros==1)
        {
            aleatorio=rand()%CANTIDAD_S;
            aux=0;
            for(i=0;i<CANTIDAD_M;i++)
            {
                aux+=chrom->vik[i]*inputs[(i*CANTIDAD_N)+indice];

            }
            aux2=0;
                for(i=0;i<CANTIDAD_S;i++)
            {
                if(i!=aleatorio)aux2+=chrom->urk[i]*outputs[(i*CANTIDAD_N)+indice];
            }
            chrom->urk[aleatorio]=(aux-aux2)/outputs[(aleatorio*CANTIDAD_N)+indice];
            calculo_D(chrom,inputs,outputs);
            val=0.0;
            val+=check_constraint_4(chrom, inputs, outputs, col , precision);
            val+=check_constraint_13(chrom);
            val+=check_constraint_5(chrom);
            val+=check_constraint_6(chrom);


        //En caso de que tengamos mas ecuaciones utilizaremos la descomposicion LU mediante
        //funciones LAPACK
        }else{
            fflush(stdout);
            matriz = (double*)malloc(dceros*dceros*sizeof(double));
            vector = (double*) malloc (dceros*dceros*sizeof(double));
            ipiv = (int*) malloc (dceros*sizeof(int));
            if(matriz==NULL){
                printf("Fallo en malloc matriz\n");
            }
            if(vector==NULL){
                printf("Fallo en malloc vector\n");
            }
            fflush(stdout);
            k=0;
            for(j=0;j<CANTIDAD_N;j++)
            {
                if(chrom->discrete_model[j]==0)
                {
                    aux=0;
                    for(i=0;i<(CANTIDAD_M+CANTIDAD_S-dceros);i++)
                    {
                        if(i<CANTIDAD_M)
                        {
                            aux+=chrom->vik[i]*inputs[(i*CANTIDAD_N)+j];
                           // printf("\nInputs[%d] = %f",(i*CANTIDAD_N)+j,inputs[(i*CANTIDAD_N)+j]);
                        }else{
                            aux-=chrom->urk[i-CANTIDAD_M]*outputs[((i-CANTIDAD_M)*CANTIDAD_N)+j];

                        }
                    }
                    vector[k]=aux;
                    //printf("\n\n Vector[%d] = %f \n",k, vector[k]);
                    p=0;
                    u=0;
                    for(i=0;i<dceros;i++)
                    {
                        if((dceros-i)>CANTIDAD_S)
                        {
                            matriz[(u*dceros)+k]=0-inputs[((CANTIDAD_M-(dceros-CANTIDAD_S-i))*CANTIDAD_N)+j];
                        }else{
                            matriz[(u*dceros)+k]=outputs[(p*CANTIDAD_N)+j];
                            //printf("\nOutputs[%d] = %f",(p*CANTIDAD_N)+j,outputs[(p*CANTIDAD_N)+j]);
                            p++;
                        }
                        u++;
                    }
                    k++;
                }
            }
            n=dceros;
            ldb=dceros;
            lda=dceros;
            nrhs=dceros;
            dgesv_( &n, &nrhs, matriz, &lda, ipiv, vector,
                        &ldb, &info);
            k=0;
            for(i=0;i<dceros;i++)
            {
                if((dceros-i)>CANTIDAD_S)
                {
                    chrom->vik[CANTIDAD_M-(dceros-CANTIDAD_S-i)]=vector[i];
                }else{
                    chrom->urk[k]=vector[i];
                    //printf("\n URK[%d] = %f",k,chrom->urk[k]);
                    k++;
                }
            }
           /* for(i=0;i<CANTIDAD_S;i++){
                if(chrom->urk[i]<0){
                    mejorar_u(chrom,inputs,outputs);
                }
            }*/
            //Comprobar t- si hay alguna negativa y en ese caso recalcular valores.

            calculo_D(chrom,inputs,outputs);
            for(i=0;i<CANTIDAD_N;i++){
                if(chrom->djk[i]<0.000001 && chrom->djk[i]>-0.000001) chrom->djk[i]=0.0;
            }
            val=0.0;
            val+=check_constraint_4(chrom, inputs, outputs, col , precision);
            val+=check_constraint_13(chrom);
            val+=check_constraint_5(chrom);
            val+=check_constraint_6(chrom);
            //printf("\n\nVALIDO = %f",val);

            free(matriz);
            free(vector);
            free(ipiv);
        }
        maxrepeticiones++;
        //system("pause");
    }
    return;
}

void mejorar_u(Ptr_Chromosome chrom, double *inputs, double *outputs){

    int i;
    double maximo=0.0;

    for(i=0;i<CANTIDAD_S;i++){
        if(maximo<chrom->urk[i]) {
                maximo=chrom->urk[i];
        }
    }
    maximo/=10.0;
    maximo+=1.0;
    //printf("\nMaximo URK = %f ////Indice = %d", maximo,indice);
    //system("pause");
    return;
}

/**
 * @brief Función para la modificación de betas y alfas para cumplir las restricciones 1.2,1.3,1.11 y 1.12,
 * si hay alguna t+ o t- negativa aux1 o aux2 valdrá uno y por lo tanto seguirán modificandose alfas y betas.
 * Si tenemos tanto alguna t+ como alguna t- negativa pasará a modificar directamente las alfas
 * @param Chrom , es el cromosoma
 * @param *inputs, vector con las entradas
 * @param *outputs, vector con las salidas
 * @param col, variable auxiliar para fijar alguna columna
 * @param alfamenor, indice del vector que contiene las alfas el cual apunta al alfa mas eficiente de modificar para satisfacer la restriccion 1
 * @param *Vx_columna, vector con los sumatorios en columnas de las entradas
 * @param *Vy_columna, vector con los sumatorios en columnas de las salidas
 */
void modify_beta(Ptr_Chromosome chrom, double *inputs, double *outputs, int col,
                  int alfamenor, double *Vx_columna, double *Vy_columna, double *outputs_inv)
{
    int aux1 = 0,
        aux2 = 0,
        d = 1;
    double factor = Random()/5.0;
    //printf("\nfactor = %f", factor);
    //Comprobamos que las t+ y t- son todas >=0 y vemos el número de alfas distintas de 0
    aux1=check_constraint_11(chrom);
    aux2=check_constraint_12(chrom);
    //printf("\n -----------------aux1 = %d-----------------------------",aux1);
    //printf("\n -----------------aux2 = %d-----------------------------",aux2);
    //system("pause");
    while((aux1!=0||aux2!=0)&&d<100){
        if(aux1!=0&&aux2==0){
               // printf("\nMODIFICAR BETA SUMANDO (BETA = %f)", chrom->beta);
                //printf("\nBETA1 = %f", chrom->beta);
            chrom->beta += factor;
            /*m2=0;
            if(m1==1)d=10;*/
            if(chrom->beta>=1.00){
                    //printf("\nEEEERROOOOOOR!!!!");
                chrom->beta=0.5;
                //printf("\nDENTRO");
               // m1 = 1;
            }

        }
        if(aux1==0&&aux2!=0){
                //printf("\nMODIFICAR BETA RESTANDO (BETA = %f)", chrom->beta);
                //printf("\nBETA2 = %f", chrom->beta);
            chrom->beta -= factor;
           /* m1=0;
            if(m2==1)d=10;*/
            if(chrom->beta<=0.00){
                    //printf("\nEEEERROOOOOOR!!!!");
                chrom->beta=0.5;
                //m2=1;
            }
        }
        if(aux1!=0&&aux2!=0){
            //printf("\nMODIFICAR BETA DOBLE (BETA = %f)", chrom->beta);
            /*modify_alfas_2(chrom, inputs, outputs, col, Vx_columna, Vy_columna);
            deducir_valores_alfa(chrom, outputs, col, alfamenor,outputs_inv,Vx_columna,Vy_columna);
            d++;*/
            chrom->beta+=factor;
        }

        //Volvemos a cálcular las t+ y t-
        modify_constraint_1(chrom, inputs, col);
        modify_constraint_2(chrom, outputs, col);

        //Volvemos a comprobar las condiciones
        d++;
        aux1=0;
        aux2=0;
        aux1=check_constraint_11(chrom);
        aux2=check_constraint_12(chrom);
        //printf("\n |||||||||||||||||||| aux1 = %d-----------------------------",aux1);
        //printf("\n |||||||||||||||||||| aux2 = %d-----------------------------",aux2);
        //system("pause");
    }
    return;
}


/**
 * @brief Función para la modificación de betas,t+ y t-. Es llamada desde la función de cruce por lo que
 * no se recalculan las alfas.
 * @param Chrom , es el cromosoma
 * @param *inputs, vector con las entradas
 * @param *outputs, vector con las salidas
 * @param col, variable auxiliar para fijar alguna columna
 */
void modify_beta_cruce(Ptr_Chromosome chrom, double *inputs, double *outputs, int col)
{
    int aux1 = 0,
        aux2 = 0,
        d = 1,
        m1 = 0,
        m2 = 0;
    double factor = 0.55;
    //Comprobamos que las t+ y t- son todas >=0 y vemos el número de alfas distintas de 0
    aux1=check_constraint_11(chrom);
    aux2=check_constraint_12(chrom);

    while((aux1!=0||aux2!=0)&&d<100){
        if(aux1!=0&&aux2==0){
            chrom->beta += factor/d;
            m2=0;
            if(m1==1)d=100;
            if(chrom->beta>=1.00){
                chrom->beta=0.999;
                m1 = 1;
            }

        }
        if(aux1==0&&aux2!=0){
            chrom->beta -= factor/d;
            m1=0;
            if(m2==1)d=100;
            if(chrom->beta<=0.00){
                chrom->beta=0.0;
                m2=1;
            }
        }
        if(aux1!=0&&aux2!=0){
            d=99;
        }

        //Volvemos a cálcular las t+ y t-
        modify_constraint_1(chrom, inputs, col);
        modify_constraint_2(chrom, outputs, col);

        //Volvemos a comprobar las condiciones
        d++;
        aux1=0;
        aux2=0;
        aux1=check_constraint_11(chrom);
        aux2=check_constraint_12(chrom);
    }
    return;
}
/**
 * @brief Función donde sacaremos los peores casos (entre las columnas) de las entradas y las salidas respecto a los sumatorios
 * de las restricciones 2 y 3. Luego sacaremos el mejor caso de entre los peores que sacamos anteriormente que será el que modificaremos para
 * cumplir con las restricciones.
 * @param Chrom , es el cromosoma
 * @param *inputs, vector con las entradas
 * @param *outputs, vector con las salidas
 */
int mejor_caso(Ptr_Chromosome chrom, double *inputs, double *outputs)
{
    int i,alfamenor,j;
    double *Peor_columnasx,
           *Peor_columnasy,
           aux = 999999;

    Peor_columnasx = (double*)malloc(CANTIDAD_N*sizeof(double));
    Peor_columnasy = (double*)malloc(CANTIDAD_N*sizeof(double));

    //Cogemos los peores casos para cada columna
    for(i=0;i<CANTIDAD_N;i++)
    {
        Peor_columnasx[i]=inputs[i];
        for(j=0;j<CANTIDAD_M;j++)
        {
            if(Peor_columnasx[i]<inputs[(j*CANTIDAD_N)+i])Peor_columnasx[i]=inputs[(j*CANTIDAD_N)+i];
        }
    }

    for(i=0;i<CANTIDAD_N;i++)
    {
        Peor_columnasy[i]=outputs[i];
        for(j=0;j<CANTIDAD_S;j++)
        {
            if(Peor_columnasy[i]<outputs[(j*CANTIDAD_N)+i])Peor_columnasy[i]=outputs[(j*CANTIDAD_N)+i];
        }
    }

    //Cogemos el mejor caso para modificar el alfa para la restricción 1 y lo guardamos en alfamenor
    for(i=0;i<CANTIDAD_N;i++)
    {
        if(((Peor_columnasx[i]+Peor_columnasy[i])<aux)&&(chrom->alfas[i]!=0.0))
        {
            aux=Peor_columnasx[i]+Peor_columnasy[i];
            alfamenor=i;
        }
    }
    free(Peor_columnasx);
    free(Peor_columnasy);

    return alfamenor;
}


/**
 * @brief Función donde queremos conseguir tener en dos vectores todos los sumatorios para todas las entradas y salidas de las
 * restrucciones 1.2 y 1.3
 * @param Chrom , es el cromosoma
 * @param *inputs, vector con las entradas
 * @param *outputs, vector con las salidas
 * @param *Vx2, vector donde tendremos los sumatorios de la restricción 1.2
 * @param *Vy2, vector donde tendremos los sumatorios de la restriccion 1.3
 */
void sumatorios(Ptr_Chromosome chrom, double *inputs, double *outputs, double *Vx, double *Vy)
{
    int i,k;

    for(i=0;i<CANTIDAD_M;i++){
            Vx[i]=0.0;
        for(k=0;k<CANTIDAD_N;k++){
            Vx[i]+=chrom->alfas[k]*inputs[(i*CANTIDAD_N)+k];
        }
    }
    for(i=0;i<CANTIDAD_S;i++){
            Vy[i]=0.0;
        for(k=0;k<CANTIDAD_N;k++){
            Vy[i]+=chrom->alfas[k]*inputs[(i*CANTIDAD_N)+k];
        }
    }

    return;
}

/**
 * Calculo de los sumatorios en columnas de las entradas y salidas
 * @param Chrom , es el cromosoma
 * @param *inputs, vector con las entradas
 * @param *outputs, vector con las salidas
 * @param *Vx_fila, vector con los sumatorios en filas de las entradas
 */
void suma_filas(Ptr_Chromosome chrom, double *inputs, double *Vx_fila)
{
    int i,j;

    //Calculo de los sumatorios de los valores de entrada (de cada fila)
    for(i=0;i<CANTIDAD_M;i++)
    {
        Vx_fila[i]=0.0;
        for(j=0;j<CANTIDAD_N;j++)
        {
            if(chrom->discrete_model[j]==0)Vx_fila[i] = Vx_fila[i] + inputs [(i*CANTIDAD_N)+j];
        }
    }
    return;
}

/**
 * Calculo de los sumatorios en columnas de las entradas y salidas
 * @param *inputs, vector con las entradas
 * @param *outputs, vector con las salidas
 * @param *Vx_columna, vector con los sumatorios en columnas de las entradas
 * @param *Vy_columna, vector con los sumatorios en columnas de las salidas
 */
void suma_columnas(double *inputs, double *outputs, double *Vx_columna, double *Vy_columna)
{
    int i,j;
    //Calculo de los sumatorios de los valores de entrada y salida (de cada columna), normalizados
    for(i=0;i<CANTIDAD_N;i++)
    {
        Vx_columna[i]=0.0;
        for(j=0;j<CANTIDAD_M;j++)
        {
            Vx_columna[i]+=inputs[(j*CANTIDAD_N)+i];
        }
        Vx_columna[i]=Vx_columna[i]/CANTIDAD_M;
    }

    for(i=0;i<CANTIDAD_N;i++)
    {
        Vy_columna[i]=0.0;
        for(j=0;j<CANTIDAD_S;j++)
        {
            Vy_columna[i]+=outputs[(j*CANTIDAD_N)+i];
        }
        Vy_columna[i]=Vy_columna[i]/CANTIDAD_S;
    }
    return;
}


/**
 * Calculo del término mas eficiente para la modificación de los alfas. Se modifica el alfa que
 * mas influya en la restricción a cumplir y que menos influya en la restricción opuesta (1.2 y 1.3).
 * @param Chrom , es el cromosoma
 * @param *inout, vector con las entradas o salidas.
 * @param *Vxy_columna, vector con los sumatorios en columnas de las entradas o salidas.
 * @param i, variable utilizada en los bucles for donde se llama a esta función.
 */
int termino_inout(Ptr_Chromosome chrom, double *inout, double *Vxy_columna, int i)
{
    int j,
        columna,
        mayor=-999999;;

    for(j=0;j<CANTIDAD_N;j++)
    {
        if(((inout[(i*CANTIDAD_N)+j]-Vxy_columna[j])>mayor)&&(chrom->alfas[j]>0.000002))
        {
            mayor=inout[(i*CANTIDAD_N)+j]-Vxy_columna[j];
            columna=j;
        }
    }
    return columna;
}


/**
 * Calculo de los terminos Djk
 * @param Chrom , es el cromosoma
 * @param *inputs, vector con las entradas
 * @param *outputs, vector con las salidas
 */
void calculo_D(Ptr_Chromosome chrom, double *inputs, double *outputs)
{
    int i,k;
    double *vecaux;

    vecaux = (double*)malloc(CANTIDAD_N*sizeof(double));

    //Utilizamos esta función para la multiplicación de entradas por vik y salidas por urk, para el calculo de las D en la 1.4
    for(i=0;i<CANTIDAD_N;i++){
            chrom->djk[i]=0.0;
        for(k=0;k<CANTIDAD_M;k++){
            chrom->djk[i]+=chrom->vik[k]*inputs[(k*CANTIDAD_N)+i];
        }
    }

    for(i=0;i<CANTIDAD_N;i++){
            vecaux[i]=0.0;
        for(k=0;k<CANTIDAD_S;k++){
            vecaux[i]+=chrom->urk[k]*outputs[(k*CANTIDAD_N)+i];
        }
    }

    //Sumamos los resultados para el calculo de todas las djk
    for(i=0;i<CANTIDAD_N;i++){
        chrom->djk[i]-=vecaux[i];
    }

    free(vecaux);
    return;
}

void inversos(double *inputs, double *outputs, double *inputs_inv, double *outputs_inv)
{
    int i;

    for(i=0;i<CANTIDAD_M*CANTIDAD_N;i++){
        inputs_inv[i]=1/inputs[i];
    }
    for(i=0;i<CANTIDAD_S*CANTIDAD_N;i++){
        outputs_inv[i]=1/outputs[i];
    }
    return;
}

void mejora(Ptr_Chromosome chrom,double *inputs,double *outputs,int col,double *inputs_inv,double *outputs_inv, int intensificacion)
{
    int i,j=0,aux=0,ind_aux=0,evaluacion_max=1;

    evaluacion_max=(int)(evaluate_chromosome(chrom,inputs,outputs,col,inputs_inv,outputs_inv));

    while(evaluacion_max!=0 && j<10){

        //printf("\n-----------1----------------");
        //show_chromosome(chrom,inputs,outputs,col,inputs_inv,outputs_inv);
        //system("pause");

            switch(evaluacion_max){

                case 1:
                    while(ind_aux==0 && aux<intensificacion)
                    {
                    chrom->beta+=0.01;
                    recalcular(chrom,inputs,outputs,col,inputs_inv,outputs_inv,1);
                    aux++;
                    ind_aux=good_chromosome(chrom,inputs,outputs,col,inputs_inv,outputs_inv);
                    }
                    aux=0;
                break;
                case 2:
                    while(ind_aux==0 && aux<intensificacion)
                    {
                    chrom->beta+=0.01;
                    recalcular(chrom,inputs,outputs,col,inputs_inv,outputs_inv,1);
                    aux++;
                    ind_aux=good_chromosome(chrom,inputs,outputs,col,inputs_inv,outputs_inv);
                    }
                    aux=0;
                break;
                case 3:
                    while(ind_aux==0 && aux<intensificacion)
                    {
                    chrom->beta+=0.01;
                    recalcular(chrom,inputs,outputs,col,inputs_inv,outputs_inv,1);
                    aux++;
                    ind_aux=good_chromosome(chrom,inputs,outputs,col,inputs_inv,outputs_inv);
                    }
                    aux=0;
                break;
                case 4:
                    while(ind_aux==0 && aux<intensificacion)
                    {
                    deducir_valores_uv(chrom,inputs,outputs,col);
                    aux++;
                    ind_aux=good_chromosome(chrom,inputs,outputs,col,inputs_inv,outputs_inv);
                    }
                    aux=0;
                break;
                case 5:
                    while(ind_aux==0 && aux<intensificacion)
                    {
                    deducir_valores_uv(chrom,inputs,outputs,col);
                    aux++;
                    ind_aux=good_chromosome(chrom,inputs,outputs,col,inputs_inv,outputs_inv);
                    }
                    aux=0;
                break;
                case 6:
                    while(ind_aux==0 && aux<intensificacion)
                    {
                    deducir_valores_uv(chrom,inputs,outputs,col);
                    aux++;
                    ind_aux=good_chromosome(chrom,inputs,outputs,col,inputs_inv,outputs_inv);
                    }
                    aux=0;
                break;
                case 11:
                    while(ind_aux==0 && aux<intensificacion)
                    {
                    for(i=0;i<CANTIDAD_S;i++){
                        chrom->tpos[i]-=0.01;
                    }
                    recalcular(chrom,inputs,outputs,col,inputs_inv,outputs_inv,3);
                    aux++;
                    ind_aux=good_chromosome(chrom,inputs,outputs,col,inputs_inv,outputs_inv);
                    }
                    aux=0;
                break;
                case 12:
                    //chrom->tneg[ind_aux]-=0.001;
                break;
                case 13:
                    while(ind_aux==0 && aux<intensificacion)
                    {
                    deducir_valores_uv(chrom,inputs,outputs,col);
                    aux++;
                    ind_aux=good_chromosome(chrom,inputs,outputs,col,inputs_inv,outputs_inv);
                    }
                    aux=0;
                break;
                case 14:
                    chrom->tpos[ind_aux]-=0.001;
                    recalcular(chrom,inputs,outputs,col,inputs_inv,outputs_inv,3);
                break;
                default:
                break;
            }
        //printf("\nValor de mejora = %d",evaluacion_max);
        if(good_chromosome(chrom,inputs,outputs,col,inputs_inv,outputs_inv)==1){
            chrom->fi += 1000;
            printf("(%d) ", evaluacion_max);
        }
        j++;
        evaluacion_max=(int)(evaluate_chromosome(chrom,inputs,outputs,col,inputs_inv,outputs_inv));
    }
    return;
}


void recalcular(Ptr_Chromosome chrom, double *inputs, double *outputs, int col,double *inputs_inv,double *outputs_inv, int nivel)
{

double precision = 0.000001;


    switch(nivel){
        case 1:
            chrom->tpos[0]=chrom->tpos[0]+((check_constraint_1(chrom, outputs, col, precision,outputs_inv)*outputs[col]))*CANTIDAD_S;
            calculo_alfas(chrom,outputs,col);
            modify_constraint_1(chrom, inputs, col);
            deducir_valores_uv(chrom,inputs,outputs,col);
            break;
        case 2:
            chrom->tpos[0]=chrom->tpos[0]+((check_constraint_1(chrom, outputs, col, precision,outputs_inv)*outputs[col]))*CANTIDAD_S;
            calculo_alfas(chrom,outputs,col);
            modify_constraint_1(chrom, inputs, col);
            deducir_valores_uv(chrom,inputs,outputs,col);
            break;
        case 3:
            generate_beta2(chrom, outputs,col);
            calculo_alfas(chrom,outputs,col);
            modify_constraint_1(chrom, inputs, col);
            break;
        case 4:
            generate_beta2(chrom, outputs,col);
            calculo_alfas(chrom,outputs,col);
            modify_constraint_1(chrom, inputs, col);
            break;
        case 5:
            deducir_valores_uv2(chrom, inputs, outputs, col);
            break;
        case 6:
            generate_beta2(chrom, outputs,col);
            calculo_alfas(chrom,outputs,col);
            modify_constraint_1(chrom, inputs, col);
            break;
        case 7:
            generate_beta2(chrom, outputs,col);
            calculo_alfas(chrom,outputs,col);
            modify_constraint_1(chrom, inputs, col);
            break;
        case 8:
            generate_beta2(chrom, outputs,col);
            calculo_alfas(chrom,outputs,col);
            modify_constraint_1(chrom, inputs, col);
            deducir_valores_uv2(chrom, inputs, outputs, col);
            break;
        case 9:
            generate_beta2(chrom, outputs,col);
            calculo_alfas(chrom,outputs,col);
            modify_constraint_1(chrom, inputs, col);
            break;
        default:
            break;
    }
    return;
}

void generate_beta (Ptr_Chromosome chrom, double *outputs, int col)
{
    int i = 0,
        j = 0,
        u = 0,
        a=0,
        decision=0,
        valido = 0;
    double aux = 10.0,
           aux2 = 0.0,
           sum = 0.0;

 decision=rand()%10;

 if(decision<5){
     for(i=0;i<CANTIDAD_S;i++){

                aux=Random()*((double)((rand()%1000)+1));
                aux2=Random();
                chrom->tpos[i]=(aux2/aux);
                }
        if(decision<2){
            chrom->tpos[rand()%CANTIDAD_S]=0.0;
        }

        while(valido==0&&u<100){
            for(i=0;i<CANTIDAD_S;i++){
                sum+=(chrom->tpos[i]/outputs[(CANTIDAD_N*i)+col]);
            }
            sum/=CANTIDAD_S;
            chrom->beta=1.0-sum;
            if(chrom->beta>0.0&&chrom->beta<=1.0){
                valido=1;

            }

            a=rand()%CANTIDAD_S;
            if(chrom->beta<=0.0)chrom->tpos[a]-=Random();
            if(chrom->beta>1.0)chrom->tpos[a]+=Random();

            j++;
            if(j==100){
                j=0;
                aux=Random()*((double)((rand()%1000)+1));
                for(i=0;i<CANTIDAD_S;i++){
                    chrom->tpos[i]=(aux2/aux);
                }
                u++;
            }
        }
 }else{

    for(i=0;i<CANTIDAD_S;i++){

            chrom->tpos[i]=(double)(rand()%1000);
            }

    if(decision<7){
            chrom->tpos[rand()%CANTIDAD_S]=0.0;
        }

    while(valido==0&&u<100){
        for(i=0;i<CANTIDAD_S;i++){
            sum+=(chrom->tpos[i]/outputs[(CANTIDAD_N*i)+col]);
        }
        sum/=CANTIDAD_S;
        chrom->beta=1.0-sum;
        if(chrom->beta>0.0&&chrom->beta<=1.0){
            valido=1;

        }

        a=rand()%CANTIDAD_S;
        if(chrom->beta<=0.0)chrom->tpos[a]-=Random();
        if(chrom->beta>1.0)chrom->tpos[a]+=Random();

        j++;
        if(j==100){
            j=0;
            for(i=0;i<CANTIDAD_S;i++){
            chrom->tpos[i]=(double)(rand()%1000);
            }
            u++;
        }
    }
 }
    return;
}

void generate_beta2 (Ptr_Chromosome chrom, double *outputs, int col)
{
    int i = 0,
        j = 0,
        u = 0,
        a=0,
        valido = 0;
    double sum = 0.0;


    while(valido==0&&u<100){
        for(i=0;i<CANTIDAD_S;i++){
            sum+=(chrom->tpos[i]/outputs[(CANTIDAD_N*i)+col]);
        }
        sum/=CANTIDAD_S;
        chrom->beta=1.0-sum;
        if(chrom->beta>=0.0&&chrom->beta<=1.0){
            valido=1;
            //system("pause");
        }

        a=rand()%CANTIDAD_S;
        if(chrom->beta<=0.0)chrom->tpos[a]/=10.0;
        if(chrom->beta>1.0)chrom->tpos[a]*=10.0;
        //printf("\TPOS[0] = %f", chrom->tpos[0]);
        j++;
        if(j==100){
            j=0;
            u++;
        }
    }
    return;
}

void calculo_alfas (Ptr_Chromosome chrom, double *outputs, int col) {
    int i = 0,
        j = 0,
        alfnocero=0,
        k=0;
    int *indices;
    double *matriz,*vector;
    int info, ldb, lda, nrhs, n;
    int *ipiv;

    //Sumatorios alfa*entrada y alfa*salida
    indices = (int*)calloc((CANTIDAD_M),sizeof(int));

    k=0;
    for(i=0;i<CANTIDAD_N;i++)
    {
        chrom->alfas[i]=0.0;
        if(chrom->discrete_model[i]==0){
            alfnocero++;
            indices[k]=i;
            //printf("\nIndices[%d] = %d",k,i);
            k++;
        }
    }

    matriz = (double*)malloc(CANTIDAD_S*CANTIDAD_S*sizeof(double));
    vector = (double*) calloc (CANTIDAD_S*CANTIDAD_S,sizeof(double));
    ipiv = (int*) malloc (CANTIDAD_S*sizeof(int));

    if(alfnocero==CANTIDAD_S)
    {
        for(i=0;i<CANTIDAD_S;i++){
            vector[i]=(chrom->beta*outputs[i*CANTIDAD_N+col])+chrom->tpos[i];
        }

        for(i=0;i<CANTIDAD_S;i++)
        {
            for(j=0;j<CANTIDAD_S;j++){
                matriz[(j*alfnocero)+i]=outputs[(i*CANTIDAD_N)+indices[j]];
                //printf("\nMatriz[%d] = output [%d]",(j*alfnocero)+i,(i*CANTIDAD_N)+indices[j]);
            }
        }
    //system("pause");
    }else{
        k=CANTIDAD_S;
        for(i=0;i<(alfnocero-CANTIDAD_S);i++){
            chrom->alfas[indices[i]]=Random();
        }
        for(i=0;i<CANTIDAD_S;i++){
            vector[i]=(chrom->beta*outputs[i*CANTIDAD_N+col])+chrom->tpos[i];
            for(j=0;j<(alfnocero-CANTIDAD_S);j++){
                vector[i]-=(chrom->alfas[indices[j]]*outputs[(i*CANTIDAD_N)+indices[j]]);
            }
        }
        for(i=0;i<CANTIDAD_S;i++)
        {
            k=alfnocero-CANTIDAD_S;
            for(j=0;j<CANTIDAD_S;j++){
                matriz[(j*CANTIDAD_S)+i]=outputs[(i*CANTIDAD_N)+indices[k]];
                k++;
            }
        }
    }

    n=CANTIDAD_S;
    ldb=CANTIDAD_S;
    lda=CANTIDAD_S;
    nrhs=CANTIDAD_S;
    dgesv_( &n, &nrhs, matriz, &lda, ipiv, vector,
                &ldb, &info);

    k=alfnocero-CANTIDAD_S;
    for(i=0;i<CANTIDAD_S;i++)
    {
        chrom->alfas[indices[k]]=vector[i];
        k++;
        //printf("\nAlfa [%d] = %f",indices[i],vector[i]);
    }
    k=0;
    //system("pause");
    free(matriz);
    free(vector);
    free(ipiv);
    free(indices);

    return;
}

double maximizar(Ptr_Chromosome chrom, double *inputs, int col)
{
    double aux = 0.0;
    int i = 0;
    for ( i = 0 ; i < CANTIDAD_M ; i++ )
    {
        aux += chrom->tneg[i]/inputs[(i*CANTIDAD_N)+col];
    }
    aux /= CANTIDAD_M;
    aux = chrom->beta - aux;

return aux;
}

void suma_columnas2(double *inputs, double *outputs, double *sumatorio_entradas, double *sumatorio_salidas)
{

    int i,j;
    double sumatorio=0.0;

    for(i=0;i<CANTIDAD_N;i++){
        sumatorio=0.0;
            for(j=0;j<CANTIDAD_S;j++){
                sumatorio+=outputs[CANTIDAD_N*j+i];
                //printf("\nOutputs[%d]",CANTIDAD_N*j+i);
            }
            sumatorio_salidas[i]=sumatorio;
            //printf("\nSumatorio = %f > Maximo = %f ?????",sumatorio,maximo);
         }
        for(i=0;i<CANTIDAD_N;i++){
        sumatorio=0.0;
            for(j=0;j<CANTIDAD_M;j++){
                sumatorio+=inputs[CANTIDAD_N*j+i];
                //printf("\nOutputs[%d]",CANTIDAD_N*j+i);
            }
            sumatorio_entradas[i]=sumatorio;
            //printf("\nSumatorio = %f > Maximo = %f ?????",sumatorio,maximo);
         }
    return;
}

Ptr_Chromosome cromosoma_vacio(etiq){

    Ptr_Chromosome chrom;

    chrom = (Ptr_Chromosome)malloc(sizeof(struct Chromosome));

    //gettimeofday(&ti, NULL);
    chrom->etiq = etiq;
    chrom->metodo = 0;
    chrom->valido = 0;
    chrom->fi=1;
    chrom->evaluation = BAD_CHROM;
    chrom->vik = (double*)malloc(CANTIDAD_M*sizeof(double));
    chrom->urk = (double*)malloc(CANTIDAD_S*sizeof(double));
    chrom->tneg = (double*)malloc(CANTIDAD_M*sizeof(double));
    chrom->tpos = (double*)malloc(CANTIDAD_S*sizeof(double));
    chrom->djk = (double*)malloc(CANTIDAD_N*sizeof(double));
    chrom->alfas = (double*)malloc(CANTIDAD_N*sizeof(double));
    chrom->discrete_model = (int*)malloc(CANTIDAD_N*sizeof(int));

    chrom->contador=0;
    chrom->terminar=0;

    return chrom;
}
