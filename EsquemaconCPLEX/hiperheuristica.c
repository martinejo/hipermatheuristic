
#include <stdio.h>
#include <stdlib.h>
#include "genetico.h"

void generar_meta(Ptr_metaheuristica *list, Ptr_hiperheuristica hiper){

    int j,k,i,p,aux=0;
	Ptr_metaheuristica meta2;
	double valor,valor2;

	/*----------- ORDENADOR DE MAYOR A MENOS FITNESS -----------*/

    for(k=0;k<hiper->NIM;k++){
                    for( j = k+1 ; j < hiper->NIM ; j++ )
                    {
                            meta2=list[k];
                            valor= meta2->fitness;
                            valor2=list[j]->fitness;
                            if (valor < valor2)
                            {
                                list[k]=list[j];
                                list[j]=meta2   ;
                            }
                        }
                    }

     /*-------------------- COMBINAR ---------------------------*/

    p=0;

    for(k=hiper->NIM;k<(hiper->NIM+hiper->NFM);k++){

    Ptr_metaheuristica meta = NULL;
    meta = (Ptr_metaheuristica)malloc(sizeof(struct metaheuristica));

    meta->INEIni = list[rand()%(hiper->NMS)]->INEIni;
    meta->FNEIni = list[rand()%(hiper->NMS)]->FNEIni;
    meta->NBESel = list[rand()%(hiper->NMS)]->NBESel;
    meta->NWESel = list[rand()%(hiper->NMS)]->NWESel;
    meta->PBBCom = list[rand()%(hiper->NMS)]->PBBCom;
    meta->PBWCom = list[rand()%(hiper->NMS)]->PBWCom;
    meta->PWWCom = list[rand()%(hiper->NMS)]->PWWCom;
    meta->MNIEnd = list[rand()%(hiper->NMS)]->MNIEnd;
    meta->NIREnd = list[rand()%(hiper->NMS)]->NIREnd;
    meta->PEIIni = list[rand()%(hiper->NMS)]->PEIIni;
    meta->IIEIni = list[rand()%(hiper->NMS)]->IIEIni;
    meta->PEIImp = list[rand()%(hiper->NMS)]->PEIImp;
    meta->IIEImp = list[rand()%(hiper->NMS)]->IIEImp;
    meta->PEDImp = list[rand()%(hiper->NMS)]->PEDImp;
    meta->IIDImp = list[rand()%(hiper->NMS)]->IIDImp;
    meta->CANTIDAD_CRUCE=2;

    list[k]=meta;
    p++;
}


return;
}
