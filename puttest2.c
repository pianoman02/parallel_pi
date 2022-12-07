#include<stdio.h>
#include<stdlib.h>
#include <bsp.h>

long P;

void test(){
    bsp_begin(P);
    long s = bsp_pid();
    long *x = NULL;
    if (s==0){
        x = malloc(P*sizeof(long));
        bsp_push_reg(x,P*sizeof(long));
    }
    else
        bsp_push_reg(&x,0); // Ook als ik bsp_push_reg(x,P*sizeof(long)) uitvoer geeft het een foutmelding.
    bsp_sync();

    if (s==0){
        x[0] = s;
    }
    else{
        bsp_put(0,&s,&x,s*sizeof(long),sizeof(long));
    }

    bsp_sync();

    if (s==0)
    {
        for (long i=0; i<P; i++)
            printf("number %ld is %ld",i,x[i]);
    }
    bsp_sync();

    bsp_pop_reg(x);
    if (s==0)
        free(x);
    bsp_end();
}

int main(int argc, char **argv){
    bsp_init(test, argc, argv);
 
    /* Sequential part */
    printf("How many processors do you want to use?\n"); fflush(stdout);
    scanf("%ld",&P);
    if (P > bsp_nprocs()){
        printf("Sorry, only %u processors available.\n", bsp_nprocs());
        fflush(stdout);
        exit(EXIT_FAILURE);
    }
 
    /* SPMD part */
    test();
 
    /* Sequential part */
    exit(EXIT_SUCCESS);
 
}