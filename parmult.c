#include<stdio.h>
#include<string.h>
#include "bspedupack.h"

long P; // number of processors requested

// Functions from bspfft
void bspfft(double complex *x, long n, bool forward, double complex *w,
                long *rho_np, long *rho_p);
void bspfft_init(long n, double complex *w, long *rho_np, long *rho_p);

void readvariable(double complex *x, char* st, long n, long decimalsPerRadix){
    /* We assume x has already been registered and has size n/p
       A bsp_sync is required to actually obtain the variable*/
    long p = bsp_nprocs();
    long s = bsp_pid();
    if (s==0){
        char xstr[n*decimalsPerRadix+1];
        int xlen;
        printf("Enter the number %s, please use at most %ld digits\n",st, n*decimalsPerRadix);
        scanf("%s",xstr);
        // TODO: Print an error message if number is to large
        xlen = strlen(xstr);

        // Converting the numbers to complex doubles.
        double complex *xstart = vecallocc(n);
        for (long i=0; i<n; i++){
            xstart[i] = 0;
        }
        long j=xlen-1;
        for (long i=0; i<n; i++){
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
        // Now broadcasting the numbers.
        // TODO: Make this more efficient, by bundeling the numbers in bulk
        for (long i=0; i<n; i++){
            bsp_put(i%p,&(xstart[i]),x,(i/p)*sizeof(double complex),sizeof(double complex));
        }
        vecfreec(xstart);
    }
} /*end readvariable*/

void carry_add_seq(double complex *x, long radix, long n){
    /*
       x = the variable to carry_add, stored on only one processor
       n = size of x
       radix = the radix with which the number is stored

       This function can deal with the fact that some coefficients are negative,
       but it gives undefined behaviour on vectors which represent negative numbers.
    */
    long carry = 0;
    long temp = 0;
    bool firstroundofCarryAddsDone = false;
    long i=0;
    while(!(firstroundofCarryAddsDone && carry==0)){
        temp = lround(creal(x[i]));
        temp += carry;
        x[i] = temp % radix;
        carry = temp / radix;
        if (creal(x[i])<0){
            x[i] += radix;
            carry --;
        }
        if (i==n-1)
        {
            i=0;
            firstroundofCarryAddsDone = true;
        }
        else
            i++;
    }
} /*end carry_add_seq*/

void prettyprinting(double complex *x, char* st, long n, long decimalsPerRadix){
    /* For prettyprinting, all numbers are brought back to processor 0
       Therefore, it is slow, but good for final output.

       x = the variable to print, distributed in cyclic distribution over the processors
       n = size of x
       st = name of x to be printed
       decimalsPerRadix = which power of 10 the radix is. Ex: if the radix is 100, we have decimalsPerRadix = 2*/
    long p = bsp_nprocs();
    long np = n/p;
    long s = bsp_pid();

    double complex *x0 = NULL;
    if (s==0){
        x0 = vecallocc(n);
        bsp_push_reg(x0,n*sizeof(long double));
    }
    else
        bsp_push_reg(&x0,0);

    ///////////////////////
    bsp_sync();
    ///////////////////////

    if (s==0){
        for (long i=0; i<np;i++)
            x0[i*p] = x[i];
    }
    else{
        for (long i=0; i<np; i++)
            bsp_put(0,&(x[i]),&x0,(s+i*p)*sizeof(double complex), sizeof(double complex));
    }
    ///////////////////////
    bsp_sync();
    ///////////////////////

    if (s==0){
        carry_add_seq(x0,100,n);
        // Handeling the rollovers
        char xystr[n*decimalsPerRadix+1];
        xystr[n*decimalsPerRadix] = '\0';
        long val = 0;
        for (long i=0; i<n; i++){
            val = lround(creal(x0[i]));
            xystr[decimalsPerRadix*n-1-decimalsPerRadix*i] = ((char)val%10)+'0';
            xystr[decimalsPerRadix*n-2-decimalsPerRadix*i] = ((char)val/10)+'0';
        }
        printf("%s = %s\n",st,xystr);

        bsp_pop_reg(x0);
        vecfreec(x0);
    }
    else
        bsp_pop_reg(&x0);
    fflush(stdout);
} /*end prettyprinting*/

void printvariable(double complex *x, char* st,long n){
    /* This function is used for debugging reasons */
    long np = n/bsp_nprocs();
    long s = bsp_pid();
    for (long i=0; i<np; i++){
        printf("In %ld index %ld the value of %s is %lf + %lf i\n",s,i,st,creal(x[i]),cimag(x[i]));
    }
}  /*end printvariable*/

void copy(double complex *x, double complex *y, long n){
    /* Copies x to y*/
    long np = n/bsp_nprocs();
    for (int i=0; i<np; i++){
        y[i] = x[i];
    }
} /*end copy*/

void multiply(double complex *x, double complex *y,long n, double complex *w, long*rho_np, long*rho_p, bool fast){
    /* Mulitplies two numbers x and y, which are stored in cyclic radix fasion.
       The result is stored in x.
       It is assumed that there are enough zeros in both x and y,
       (for now, at least n/2 zeros in both x and y)
       such that the multiplication doesn't cause any cyclic overflow
       
       x and y must have been registered before calling this function.
       Moreover, the weights w have to be initialised using bspfft_init*/
    long np = n/bsp_nprocs();
    double complex *y2;
    if (fast == false){
        y2 = vecallocc(np);
        bsp_push_reg(y2,np*sizeof(double complex));
        copy(y,y2,n);
    }
    else
        y2 = y;

    bspfft(x,n,true,w,rho_np,rho_p);
    bspfft(y2,n,true,w,rho_np,rho_p);

    for (int i=0; i<np; i++){
        x[i] *= y2[i];
    }
    bspfft(x,n,false,w,rho_np,rho_p);
    if (fast == false){
        bsp_pop_reg(y2);
        vecfreec(y2);
    }
} /*end multiply*/

void add(double complex *x, double complex *y, long n){
    /* Adds the two numbers by performing a pointwise addition
       The result is stored in x*/
    long np = n/bsp_nprocs();
    for (int i=0;i<np; i++){
        x[i] += y[i];
    }
} /*end add*/

void minus(double complex *x, double complex *y, long n){
    /* Calculates x-y
       Stores the result in y (NOTE!!)  */
    long np = n/bsp_nprocs();
    for (int i=0; i<np; i++){
        y[i] = x[i] - y[i];
    }
} /*end minus*/

void set_half_zero(double complex *x, long n){
    // TODO: make a range in which x is allowed to be nonzero.
    // for now, I assume the numbers will never get larger then 100 (a single digit in radix notation).
    // Therefore, i will make all the indices 1,\ldots n/2 zero.
    long s = bsp_pid();
    long np = n/bsp_nprocs();
    if (s==0){
        for (long i=1; i<= np/2; i++)
            x[i] = 0;
    }
    else
        for (long i=0; i<np/2; i++)
            x[i] = 0;
}

void carry_add(double complex *x, long *x_carry, long n){
    /* Performs one carry-add operation on the number x
       Moreover, it rounds the number to get rid of any numerical errors.
       It is assumed that x is registered before calling this function
       x_carry is a container which should be free, and registered*/
    
    /// TODO: Make the carry_add possible for negative numbers, and add the zero zone
    long p = bsp_nprocs();
    long s = bsp_pid();
    long np = n/p;
    if (s!= p-1){
        for (long i=0; i<np; i++){
            x_carry[i] = lround(creal(x[i]));
            x[i] = x_carry[i]%100;
            x_carry[i] = x_carry[i]/100;
        }
    }
    else{
        for (long i=1; i<np; i++){
            x_carry[i] = lround(creal(x[i-1]));
            x[i-1] = x_carry[i]%100;
            x_carry[i] = x_carry[i]/100;
        }
        x_carry[0] = lround(creal(x[np-1]));
        x[np-1] = x_carry[0]%100;
        x_carry[0] = x_carry[0]/100;
    }
    bsp_put((s+1)%p,x_carry, x_carry,0,np*sizeof(long));
    //////////////////////////////
    bsp_sync();
    //////////////////////////////
    for (long i=0; i<np; i++){
        x[i] += x_carry[i];
    }
} /*end carry_add*/

void setToHalf(double complex *x, long n, long radix){
    /* Sets the number x to be a half*/
    long p = bsp_nprocs();
    long s = bsp_pid();
    long np = n/p;
    for (int i=0; i<np; i++){
        x[i] = 0;
    }
    if (s==p-1)
        x[np-1] = radix/2;
} /*end setToHalf*/

void setToTwo(double complex *x, long n){
    /* Sets the number x to be 2*/
    long np = n/bsp_nprocs();
    long s = bsp_pid();
    for (int i=0; i<np; i++){
        x[i] = 0;
    }
    if (s==0)
        x[0] = 2;
}

void setToZero(double complex *x, long n){
    long np = n/bsp_nprocs();
    for (int i=0; i<np; i++){
        x[i] = 0;
    }
} /*end setToZero*/

void one_over_square_root(double complex *a, double complex *x, long *carry, long n,double complex *w, long*rho_np, long*rho_p){
    /* Returns at x the square root of x.
    */
    long M = 10;
    long np = n/bsp_nprocs();
    setToHalf(x,n, 100); // set x to initial value 1/2
    double complex *half = vecallocc(np);
    setToHalf(half,n,100);
    double complex *xprev = vecallocc(np);
    bsp_push_reg(xprev,np*sizeof(double complex));
    bsp_push_reg(half,np*sizeof(double complex));

    //////////////////////////////////
    bsp_sync();
    //////////////////////////////////

    for (int i=0; i<M; i++){
        copy(x,xprev,n);
        multiply(x,xprev,n,w,rho_np,rho_p, false);
        set_half_zero(x,n);
        carry_add(x,carry,n);
        carry_add(x,carry,n);
        multiply(x,xprev,n,w,rho_np,rho_p, false);
        set_half_zero(x,n);
        carry_add(x,carry,n);
        carry_add(x,carry,n);
        multiply(x,a,n,w,rho_np,rho_p, false);
        set_half_zero(x,n);
        carry_add(x,carry,n);
        carry_add(x,carry,n);
        minus(xprev,x,n);
        multiply(x,half,n,w,rho_np,rho_p, false);
        set_half_zero(x,n);
        carry_add(x,carry,n);
        carry_add(x,carry,n);
        add(x,xprev,n);
        carry_add(x,carry,n);
    }
    bsp_pop_reg(xprev);
    bsp_pop_reg(half);
    vecfreec(xprev);
    vecfreec(half);
} /*end one_over_sqare_root*/

void square_root(double complex *a, double complex *x, long *carry, long n,double complex *w, long*rho_np, long*rho_p){
    one_over_square_root(a,x,carry,n,w,rho_np,rho_p);
    multiply(x,a,n,w,rho_np,rho_p,false);
    carry_add(x,carry,n);
    carry_add(x,carry,n);
} /*end sqrt_root*/

void runsqrt(){
    bsp_begin(P);
    long p = bsp_nprocs();
    long n= 32; // NOTE: n should be a power of two, as well as p
    long decimalsPerRadix = 2;
    long np = n/p;

    double complex *a = vecallocc(np);
    double complex *x = vecallocc(np);
    long *carry = vecalloci(np);
    setToTwo(a,n);
    setToZero(x,n);
    bsp_push_reg(a,np*sizeof(double complex));
    bsp_push_reg(x,np*sizeof(double complex));
    bsp_push_reg(carry,np*sizeof(long));

    // Initialize the weight and bit reversal tables
    /* First, determine the number of computation supersteps, by computing
       the smallest integer t such that (n/p)^t >= p */
    long t= 0;
    for (long c=1; c<p; c *= np)
        t++;
    double complex *w= vecallocc((t+1)*np);
    long *rho_np= vecalloci(np);
    long *rho_p=  vecalloci(p);
    bspfft_init(n,w,rho_np,rho_p);

    ///////////////////
    bsp_sync();
    ///////////////////

    square_root(a,x,carry,n,w,rho_np,rho_p);
    prettyprinting(x,"sqrt(2)",n,2);

    bsp_pop_reg(carry);
    bsp_pop_reg(a);
    bsp_pop_reg(x);

    vecfreei(carry);
    vecfreec(a);
    vecfreec(x);
    vecfreei(rho_p);
    vecfreei(rho_np);
    vecfreec(w);


}

void runmult(){
    bsp_begin(P);
    long p= bsp_nprocs();
    long n = 256; // NOTE: n should be a power of two, as well as p
    long decimalsPerRadix = 2;
    long np = n/p;

    // Allocate, register, and initialize vectors 
    double complex *x= vecallocc(np);
    double complex *y= vecallocc(np);
    for (long j=0; j<np; j++){
        x[j]= 0;
        y[j]= 0;
    }
    bsp_push_reg(x,np*sizeof(double complex));
    bsp_push_reg(y,np*sizeof(double complex));
    
    // Allocate space for the carry-adds
    long *x_carry = vecalloci(np);
    long *y_carry = vecalloci(np);
    bsp_push_reg(x_carry, np*sizeof(long));
    bsp_push_reg(y_carry, np*sizeof(long));

    // Initialize the weight and bit reversal tables
    /* First, determine the number of computation supersteps, by computing
       the smallest integer t such that (n/p)^t >= p */
    long t= 0;
    for (long c=1; c<p; c *= np)
        t++;
    double complex *w= vecallocc((t+1)*np);
    long *rho_np= vecalloci(np);
    long *rho_p=  vecalloci(p);
    bspfft_init(n,w,rho_np,rho_p);

    /////////////////////////
    bsp_sync();
    /////////////////////////

    readvariable(x,"x",n,decimalsPerRadix);
    readvariable(y,"y",n,decimalsPerRadix);

    /////////////////////////
    bsp_sync();
    /////////////////////////

    multiply(x,y,n,w,rho_np,rho_p, true);
    carry_add(x,x_carry,n);
    prettyprinting(x,"x*y",n,decimalsPerRadix);

    bsp_pop_reg(x);
    bsp_pop_reg(y);
    bsp_pop_reg(x_carry);
    bsp_pop_reg(y_carry);

    vecfreec(x);
    vecfreec(y);
    vecfreei(x_carry);
    vecfreei(y_carry);
    vecfreei(rho_p);
    vecfreei(rho_np);
    vecfreec(w);

    bsp_end();
 }

int main(int argc, char **argv){
    bsp_init(runsqrt, argc, argv);
 
    /* Sequential part */
    printf("How many processors do you want to use?\n"); fflush(stdout);
    scanf("%ld",&P);
    if (P > bsp_nprocs()){
        printf("Sorry, only %u processors available.\n", bsp_nprocs());
        fflush(stdout);
        exit(EXIT_FAILURE);
    }
 
    /* SPMD part */
    runsqrt();
 
    /* Sequential part */
    exit(EXIT_SUCCESS);
 
}
