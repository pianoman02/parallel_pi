#include<stdio.h>
#include<string.h>
#include "bspedupack.h"

long P; // number of processors requested

// Functions from bspfft
void bspfft(double complex *x, long n, bool forward, double complex *w,
                long *rho_np, long *rho_p);
void bspfft_init(long n, double complex *w, long *rho_np, long *rho_p);

void printvariable(double complex *x, char* st,long np){
    /* This function is used for debugging reasons */
    long s = bsp_pid();
    for (long i=0; i<np; i++){
        printf("In %ld index %ld the value of %s is %lf + %lf i\n",s,i,st,creal(x[i]),cimag(x[i]));
    }
}  /*end printvariable*/

void multiply(double complex *x, double complex *y,long n, double complex *w){
    /* Mulitplies two numbers x and y, which are stored in cyclic radix fasion.
       It is assumed that there are enough zeros in both x and y,
       (for now, at least n/2 zeros in both x and y)
       such that the multiplication doesn't cause any cyclic overflow
       
       x and y must have been registered before calling this function.
       Moreover, the weights w have to be initialised using bspfft_init*/
    


} /*end multiply*/


void carry_add(double complex *x, long *x_carry, long n){
    /* Performs one carry-add operation on the number x
       Moreover, it rounds the number to get rid of any numerical errors.
       It is assumed that x is registered before calling this function
       x_carry is a container which should be free, and registered*/
    long p = bsp_nprocs();
    long s = bsp_pid();
    long np = n/p;
    for (long i=0; i<np; i++){
        x_carry[i] = lround(creal(x[i]));
        x[i] = x_carry[i]%100;
        x_carry[i] = x_carry[i]/100;
    }
    bsp_put((s+1)%p,x_carry, x_carry,0,np*sizeof(long));
    bsp_sync();
    for (long i=0; i<np; i++){
        x[i] += x_carry[i];
    }
} /*end carry_add*/


void run(){
    bsp_begin(P);
    long p= bsp_nprocs();
    long s= bsp_pid();
    long n = 8;
    long decimalsPerRadix = 2;
    long np = n/p;
    // NOTE: We assume that p devides n

    /* Determine the number of computation supersteps, by computing
       the smallest integer t such that (n/p)^t >= p */
    long t= 0;
    for (long c=1; c<p; c *= np)
        t++;

    /* Allocate, register,  and initialize vectors */
    double complex *w= vecallocc((t+1)*np);
    double complex *x= vecallocc(np);
    long *x_carry = vecalloci(np);
    double complex *y= vecallocc(np);
    long *y_carry = vecalloci(np);
    bsp_push_reg(x,np*sizeof(double complex));
    bsp_push_reg(y,np*sizeof(double complex));
    bsp_push_reg(x_carry, np*sizeof(long));
    bsp_push_reg(y_carry, np*sizeof(long));
    long *rho_np= vecalloci(np);
    long *rho_p=  vecalloci(p);

    for (long j=0; j<np; j++){
        x[j]= 0;
        y[j]= 0;
    }

    double complex *xstart = NULL;

    bsp_sync();

    if (s==0){
        // Reading the numbers
        char xstr[n*decimalsPerRadix/2+1], ystr[n*decimalsPerRadix/2+1];
        int xlen, ylen;
        printf("For the numbers, please use at most %ld digits\n",n*decimalsPerRadix/2);
        printf("Please enter the number x: \n");
        scanf("%s",xstr);
        printf("Please enter the number y: \n");
        scanf("%s",ystr);
        // TODO: Print an error message if number is to large
        xlen = strlen(xstr);
        ylen = strlen(ystr);
//        printf("x= %s\n",xstr);
//        printf("y= %s\n",ystr);

        // Converting the numbers to complex doubles.
        xstart = vecallocc(n); // because the x is also used for the return 
        bsp_push_reg(xstart,n*sizeof(double complex));
        double complex *ystart = vecallocc(n/2);
        for (long i=0; i<n; i++){
            xstart[i] = 0;
        }
        for (long i=0; i<n/2; i++){
            ystart[i] = 0;
        }
        long j=xlen-1;
        for (long i=0; i<n/2; i++){
            // TODO: Make more general using decimalsPerRadix
            if (j == -1)
                break;
            xstart[i] = xstr[j]-'0';
            j--;
            if (j == -1)
                break;
            xstart[i] += 10*(xstr[j]-'0');
            j--;          
        }
//        for (long i=0; i<n/2; i++)
//            printf("the %ld number is %lf\n",i,creal(xstart[i]));
        j=ylen-1;
        for (long i=0; i<n/2; i++){
            // TODO: Make more general using decimalsPerRadix
            if (j == -1)
                break;
            ystart[i] = ystr[j]-'0';
            j--;
            if (j == -1)
                break;
            ystart[i] += 10*(ystr[j]-'0');   
            j--;         
        }
//        for (long i=0; i<n/2; i++)
//            printf("the %ld number is %lf\n",i,creal(ystart[i]));

        // Now broadcasting the numbers.
        // TODO: Make this more efficient, by bundeling the numbers in bulk
        for (long i=0; i<n/2; i++){
            bsp_put(i%p,&(xstart[i]),x,(i/p)*sizeof(double complex),sizeof(double complex));
            bsp_put(i%p,&(ystart[i]),y,(i/p)*sizeof(double complex),sizeof(double complex));
        }
    }
    bsp_push_reg(xstart,n*sizeof(double complex));
    bsp_sync();
//    for (long i=0; i<np; i++){
//        printf("In %ld index %ld the value of x is %lf\n",s,i,creal(x[i]));
//    }

    /* Initialize the weight and bit reversal tables */
    bspfft_init(n,w,rho_np,rho_p);
    bsp_sync();
  
    /* Perform the FFTs */
    bspfft(x,n,true,w,rho_np,rho_p);
    bspfft(y,n,true,w,rho_np,rho_p);

    for (long i=0; i<np; i++){
        printf("In %ld index %ld the value of x is %lf + %lf i\n",s,i,creal(x[i]),cimag(x[i]));
    }
    for (long i=0; i<np; i++){
        printf("In %ld index %ld the value of y is %lf + %lf i\n",s,i,creal(y[i]),cimag(y[i]));
    }

    for (int i=0; i<np; i++){
        x[i] *= y[i];
    }
    for (long i=0; i<np; i++){
        printf("In %ld index %ld the value of x2 is %lf + %lf i\n",s,i,creal(x[i]),cimag(x[i]));
    }
    bspfft(x,n,false,w,rho_np,rho_p);
    carry_add(x,x_carry,n);
    // Now, perform rounding and carry-adds.
    
    for (long i=0; i<np; i++){
        printf("In %ld index %ld the value of x is %lf\n",s,i,creal(x[i]));
    }
    
    /////////////// Now, bring all the numbers back to processor 0 to print them.
    // TODO: Make communication more efficient
    for (long i=0; i<np; i++)
        bsp_put(0,&(x[i]),xstart,(s+i*p)*sizeof(double complex), sizeof(double complex));
    bsp_sync();
    if (s==0){
        // Handeling the rollovers
        char xystr[n*decimalsPerRadix+1];
        xystr[n*decimalsPerRadix] = '\0';
        printf("The product is\n");
        printf("%s",xystr);
        bsp_pop_reg(xstart);
        vecfreec(xstart);
        xstart=NULL;
    }
    
    fflush(stdout);
    bsp_sync();

    bsp_pop_reg(x);
    bsp_pop_reg(y);
    bsp_pop_reg(x_carry);
    bsp_pop_reg(y_carry);
    bsp_sync();

    vecfreei(rho_p);
    vecfreei(rho_np);
    vecfreec(x);
    vecfreec(y);
    vecfreei(x_carry);
    vecfreei(y_carry);
    vecfreec(w);

    bsp_end();
 }

int main(int argc, char **argv){
    bsp_init(run, argc, argv);
 
    /* Sequential part */
    printf("How many processors do you want to use?\n"); fflush(stdout);
    scanf("%ld",&P);
    if (P > bsp_nprocs()){
        printf("Sorry, only %u processors available.\n", bsp_nprocs());
        fflush(stdout);
        exit(EXIT_FAILURE);
    }
 
    /* SPMD part */
    run();
 
    /* Sequential part */
    exit(EXIT_SUCCESS);
 
}
