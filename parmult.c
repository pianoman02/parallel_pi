#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "bspedupack.h"

long P; // number of processors requested
long N; // size of the number
long DECIMALS; // amount of decimals of radix;
long M; 

// Functions from bspfft
void bspfft(double complex *x, long n, bool forward, double complex *w,
                long *rho_np, long *rho_p);
void bspfft_init(long n, double complex *w, long *rho_np, long *rho_p);

// NB: this function works not yet for radix != 100
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

void printvariable(double complex *x, char* st,long n){
    /* This function is used for debugging reasons */
    long np = n/bsp_nprocs();
    long s = bsp_pid();
    for (long i=0; i<np; i++){
        printf("In %ld index %ld the value of %s is %lf + %lf i\n",s,i,st,creal(x[i]),cimag(x[i]));
    }
}  /*end printvariable*/

void printvariable_seq(double complex *x, char* st, long n){
    for (long i=0; i<n; i++){
        printf("In index %ld the value of %s is %lf + %lf i\n",i,st,creal(x[i]),cimag(x[i]));
    }
}

void prettyprinting(double complex *x, char* st, long n, long decimalsPerRadix,long radix, bool checkpi){
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
        carry_add_seq(x0,radix,n);
        char xystr[n*decimalsPerRadix+1];
        xystr[n*decimalsPerRadix] = '\0';
        long val = 0;
        long devideby = 1;
        for (long i=0; i<n; i++){
            val = lround(creal(x0[i]));
            devideby=1;
            for (long j=1; j<= decimalsPerRadix; j++){
                xystr[decimalsPerRadix*n-j-decimalsPerRadix*i] = (char)((val/devideby)%10)+'0';
                devideby *=10;
            }
        }
        // Now print the first 1000 decimals
        char prstr[1001];
        prstr[1000] = '\0';
        for (long i=0; i<1000 && i<n*decimalsPerRadix; i++){
            prstr[i] = xystr[i];
        }
        printf("%s = %s\n",st,prstr);

        if (checkpi){
            // how to open a file is from https://www.geeksforgeeks.org/c-program-to-read-contents-of-whole-file/
            FILE* ptr;
            char ch;
            ptr = fopen("pi.txt","r");

            if (ptr == NULL){
                printf("file can't be opened \n");
            }
            long i=0;
            do {
                ch = fgetc(ptr);
                if (ch != xystr[i])
                    break;
                else
                    i++;
                // Checking if character is not EOF.
                // If it is EOF stop reading.
            } while (ch != EOF);
        
            // Closing the file
            fclose(ptr);
            printf("Congratulations, we've got %ld good decimals\n",i);
        }
        
        bsp_pop_reg(x0);
        vecfreec(x0);
    }
    else
        bsp_pop_reg(&x0);
    fflush(stdout);
} /*end prettyprinting*/

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
        // Note, we don't need to do a sync since that is already done in the first fft.
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

void square(double complex *x, long n, double complex *w, long *rho_np, long*rho_p){
    long np = n/bsp_nprocs();
    bspfft(x,n,true,w,rho_np,rho_p);
    for (int i=0; i<np; i++){
        x[i] *= x[i];
    }
    bspfft(x,n,false,w,rho_np,rho_p);
} /*end square*/

void multiply_by_2(double complex *x, long n){
    long np = n/bsp_nprocs();
    for (int i=0;i<np; i++){
        x[i] *= 2;
    }
} /*end multiply_by_2*/

void divide_by_2(double complex *x, long*carry, long n, long radix){
    /*Divides the number by two, takes one synchronisation.
      Fun fact: the carry-adds are here used in the reverse direction!*/
    long p = bsp_nprocs();
    long s = bsp_pid();
    long np = n/bsp_nprocs();
    if (s!=0){
        for (long i=0; i<np; i++){
            carry[i] = lround(creal(x[i]));
            x[i] = carry[i]/2;
            carry[i] = (carry[i]%2)*radix/2;
        }
        bsp_put(s-1,carry, carry,0,np*sizeof(long));
    }
    else{
        for (long i=0; i<np-1; i++){
            carry[i] = lround(creal(x[i+1]));
            x[i+1] = carry[i]/2;
            carry[i] = (carry[i]%2)*radix/2;
        }
        carry[np-1] = lround(creal(x[0]));
        x[0] = carry[np-1]/2;
        carry[np-1] = (carry[np-1]%2)*radix/2;
        bsp_put(p-1,carry, carry,0,np*sizeof(long));
    }
    //////////////////////////////
    bsp_sync();
    //////////////////////////////
    for (long i=0; i<np; i++){
        x[i] += carry[i];
    }
} /*end divide_by_2*/

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
       Stores the result in x  */
    long np = n/bsp_nprocs();
    for (int i=0; i<np; i++){
        x[i] -= y[i];
    }
} /*end minus*/

void minus_2(double complex *x, double complex *y, long n){
    /* Calculates x-y
       Stores the result in y (NOTE!!)  */
    long np = n/bsp_nprocs();
    for (int i=0; i<np; i++){
        y[i] = x[i] - y[i];
    }
} /*end minus_2*/

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

void carry_add(double complex *x, long *carry,long radix, long n){
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
            carry[i] = lround(creal(x[i]));
            x[i] = carry[i]%radix;
            carry[i] = carry[i]/radix;
        }
    }
    else{
        for (long i=1; i<np; i++){
            carry[i] = lround(creal(x[i-1]));
            x[i-1] = carry[i]%radix;
            carry[i] = carry[i]/radix;
        }
        carry[0] = lround(creal(x[np-1]));
        x[np-1] = carry[0]%radix;
        carry[0] = carry[0]/radix;
    }
    bsp_put((s+1)%p,carry, carry,0,np*sizeof(long));
    //////////////////////////////
    bsp_sync();
    //////////////////////////////
    for (long i=0; i<np; i++){
        x[i] += carry[i];
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

void setToInt(double complex *x, long n, int v){
    /* Sets the number x to be v. We must have v<radix*/
    long np = n/bsp_nprocs();
    long s = bsp_pid();
    for (int i=0; i<np; i++){
        x[i] = 0;
    }
    if (s==0)
        x[0] = v;
}

void one_over_square_root(double complex *a, double complex *x, long *carry,long radix, long n,double complex *w, long*rho_np, long*rho_p){
    /* Returns at x the 1/sqrt(a).
    */
    long np = n/bsp_nprocs();
    setToHalf(x,n, radix); // set x to initial value 1/2
    double complex *xprev = vecallocc(np);
    double complex *afft = vecallocc(np);
    double complex *xfft = vecallocc(np);
    bsp_push_reg(xprev,np*sizeof(double complex));
    bsp_push_reg(afft, np*sizeof(double complex));
    bsp_push_reg(xfft, np*sizeof(double complex));

    //////////////////////////////////
    bsp_sync();
    //////////////////////////////////
    copy(a,afft,n);
    bspfft(afft,n,true,w,rho_np,rho_p);

    for (int i=0; i<M; i++){
        copy(x,xprev,n);
        // Calculating ax^3
        copy(x,xfft,n);
        bspfft(xfft,n,true,w,rho_np,rho_p);
        for (int i=0; i<np; i++){
            x[i] =xfft[i]*xfft[i];
        }
        bspfft(x,n,false,w,rho_np,rho_p);
        set_half_zero(x,n);
        carry_add(x,carry,radix,n);
        carry_add(x,carry,radix,n);
        bspfft(x,n,true,w,rho_np,rho_p);
        for (int i=0; i<np; i++){
            x[i] =x[i]*xfft[i];
        }
        bspfft(x,n,false,w,rho_np,rho_p);
        set_half_zero(x,n);
        carry_add(x,carry,radix,n);
        carry_add(x,carry,radix,n);
        bspfft(x,n,true,w,rho_np,rho_p);
        for (int i=0; i<np; i++){
            x[i] =x[i]*afft[i];
        }
        bspfft(x,n,false,w,rho_np,rho_p);
        set_half_zero(x,n);
        carry_add(x,carry,radix,n);
        carry_add(x,carry,radix,n);

        // The rest
        minus_2(xprev,x,n);
        divide_by_2(x,carry,n,radix);
        set_half_zero(x,n);
        carry_add(x,carry,radix,n);
        carry_add(x,carry,radix,n);
        add(x,xprev,n);
        carry_add(x,carry,radix,n);
    }
    bsp_pop_reg(xprev);
    bsp_pop_reg(afft);
    bsp_pop_reg(xfft);
    vecfreec(xprev);
    vecfreec(afft);
    vecfreec(xfft);
} /*end one_over_sqare_root*/

void square_root(double complex *a, double complex *x, long *carry,long radix, long n,double complex *w, long*rho_np, long*rho_p){
    one_over_square_root(a,x,carry,radix,n,w,rho_np,rho_p);
    multiply(x,a,n,w,rho_np,rho_p,false);
    set_half_zero(x,n);
    carry_add(x,carry,radix,n);
    carry_add(x,carry,radix,n);
} /*end sqrt_root*/

void one_over(double complex *a, double complex *x, long *carry,long radix, long n, double complex *w, long*rho_np, long*rho_p){
    /* Returns at x the number 1/a.
    */
    long np = n/bsp_nprocs();
    setToHalf(x,n, radix); // set x to initial value 1/2
    double complex *xprev = vecallocc(np);
    double complex *afft = vecallocc(np);
    bsp_push_reg(xprev,np*sizeof(double complex));
    bsp_push_reg(afft,np*sizeof(double complex));

    //////////////////////////////////
    bsp_sync();
    //////////////////////////////////

    copy(a,afft,n);
    bspfft(afft,n,true,w,rho_np,rho_p);

    for (int i=0; i<M; i++){
        copy(x,xprev,n);
        // Calculate ax^2
        square(x,n,w,rho_np,rho_p);
        set_half_zero(x,n);
        carry_add(x,carry,radix,n);
        carry_add(x,carry,radix,n);
        bspfft(x,n,true,w,rho_np,rho_p);
        for (int i=0; i<np; i++){
            x[i] =x[i]*afft[i];
        }
        bspfft(x,n,false,w,rho_np,rho_p);
        set_half_zero(x,n);
        carry_add(x,carry,radix,n);
        carry_add(x,carry,radix,n);
        // The rest
        multiply_by_2(xprev,n);
        minus_2(xprev,x,n);
        set_half_zero(x,n);
        carry_add(x,carry,radix,n);
    }
    bsp_pop_reg(xprev);
    bsp_pop_reg(afft);
    vecfreec(xprev);
    vecfreec(afft);
} /*end one_over*/

void calculatepi(){
    bsp_begin(P);
    long p = bsp_nprocs();
    long n= N; // NOTE: n should be a power of two, as well as p
    // NOTE: if you choose n too small, then the number of steps in the loops is too large in comparison to the numbers, so you get nonsense perhaps
    long decimalsPerRadix = DECIMALS; // max 5
    long radix =1;
    for (long i=0; i<decimalsPerRadix; i++){
        radix *= 10;
    }
    long np = n/p;

    // Initialise and register
    double complex *a = vecallocc(np);
    double complex *aprev = vecallocc(np);
    double complex *b = vecallocc(np);
    double complex *c = vecallocc(np);
    double complex *d = vecallocc(np);
    double complex *power2 = vecallocc(np);

    bsp_push_reg(a,np*sizeof(double complex));
    bsp_push_reg(aprev, np*sizeof(double complex));
    bsp_push_reg(b,np*sizeof(double complex));
    bsp_push_reg(c,np*sizeof(double complex));
    bsp_push_reg(d,np*sizeof(double complex));
    bsp_push_reg(power2, np*sizeof(double complex));

    long *carry = vecalloci(np);
    bsp_push_reg(carry,np*sizeof(long));

    // Initialize the weight and bit reversal tables
    // First, determine the number of computation supersteps, by computing
    //   the smallest integer t such that (n/p)^t >= p 
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

    // Initialising values
    setToInt(aprev,n,2);
    square_root(aprev,a,carry,radix,n,w,rho_np,rho_p);
    setToInt(power2,n,1);
    setToInt(b,n,1);
    setToInt(d,n,1);
    long M = 25;
    for (long i=0; i<M; i++){
        copy(a,aprev,n);
        add(a,b,n);

        divide_by_2(a,carry,n,radix);
        set_half_zero(a,n);
        carry_add(a,carry,radix,n);

        multiply(aprev,b,n,w,rho_np,rho_p, false); 
        set_half_zero(aprev,n);
        carry_add(aprev,carry,radix,n);
        carry_add(aprev,carry,radix,n);

        square_root(aprev,b,carry,radix,n,w,rho_np,rho_p);

        multiply_by_2(power2,n);
        carry_add(power2,carry,radix,n);

        copy(a,aprev,n);
        multiply(aprev,a,n,w,rho_np,rho_p, false);
        set_half_zero(aprev,n);
        carry_add(aprev,carry,radix,n);
        carry_add(aprev,carry,radix,n);

        copy(b,c,n);
        multiply(c,b,n,w,rho_np,rho_p,false);
        set_half_zero(c,n);
        carry_add(c,carry,radix,n);
        carry_add(c,carry,radix,n);

        minus_2(aprev,c,n);

        multiply(c,power2,n,w,rho_np,rho_p, false);
        set_half_zero(c,n);
        carry_add(c,carry,radix,n);
        carry_add(c,carry,radix,n);

        minus(d,c,n);
        
    }
    square(a,n,w,rho_np,rho_p);
    set_half_zero(a,n);
    carry_add(a,carry,radix,n);
    carry_add(a,carry,radix,n);

    multiply_by_2(a,n);
    set_half_zero(a,n);
    carry_add(a,carry,radix,n);
    carry_add(a,carry,radix,n);

    one_over(d,aprev,carry,radix,n,w,rho_np,rho_p);

    multiply(a,aprev,n,w,rho_np,rho_p,true);
    set_half_zero(a,n);
    carry_add(a,carry,radix,n);
    carry_add(a,carry,radix,n);
    prettyprinting(a,"pi",n,decimalsPerRadix,radix,true);
    // Now check how many decimals I've got right

    // Deregister and free
    bsp_pop_reg(a);
    bsp_pop_reg(aprev);
    bsp_pop_reg(b);
    bsp_pop_reg(c);
    bsp_pop_reg(d);
    bsp_pop_reg(power2);

    vecfreec(a);
    vecfreec(aprev);
    vecfreec(b);
    vecfreec(c);
    vecfreec(d);
    vecfreec(power2);

    bsp_pop_reg(carry);
    vecfreei(carry);

    vecfreei(rho_p);
    vecfreei(rho_np);
    vecfreec(w);
}

void runone_over(){
    bsp_begin(P);
    long p = bsp_nprocs();
    long n= 16; // NOTE: n should be a power of two, as well as p
    long decimalsPerRadix = 3;
    long radix = 1000;
    long np = n/p;

    double complex *a = vecallocc(np);
    double complex *x = vecallocc(np);
    long *carry = vecalloci(np);
    setToInt(a,n,3);
    setToInt(x,n,0);
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

    one_over(a,x,carry,radix,n,w,rho_np,rho_p);
    printvariable(x,"1/3",n);
    prettyprinting(x,"1/3",n,decimalsPerRadix,radix,false);

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
    long radix = 100;

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
    carry_add(x,x_carry,radix,n);
    prettyprinting(x,"x*y",n,decimalsPerRadix,radix,false);

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
    bsp_init(calculatepi, argc, argv);
 
    /* Sequential part */
    printf("How many processors do you want to use?\n"); fflush(stdout);
    scanf("%ld",&P);
    if (P > bsp_nprocs()){
        printf("Sorry, only %u processors available.\n", bsp_nprocs());
        fflush(stdout);
        exit(EXIT_FAILURE);
    }
    printf("What precision n do you want to use?\n"); fflush(stdout);
    scanf("%ld",&N);
    printf("How digits should the radix have?\n"); fflush(stdout);
    scanf("%ld",&DECIMALS);
    printf("How many iterations (M)?\n");
    scanf("%ld",&M);

    /* SPMD part */
    calculatepi();
 
    /* Sequential part */
    exit(EXIT_SUCCESS);
 
}
