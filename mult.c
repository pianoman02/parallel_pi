 #include<stdio.h>
 #include"bspedupack.h"
// This program calculates the inproduct sequentially.
// First, we represent the number using chars.
#define MAXDIGITS 1000

int getline (char s[]){
    int i,c;
    for (i=0; i<MAXDIGITS-1 && (c = getchar())!=EOF && c!='\n'; i++){
        s[i] = c;
    }
    s[i] = '\0';
    return i;
}

int main(){
    char s[MAXDIGITS];
    unsigned long long a=0,b=0,c = 0; // At least 64 bits= 8 bytes
    printf("What should the first number be? (in decimal)\n");
    int sizea, sizeb;
    
    scanf("%lld", &a);
    // First, read 17 digits. Then, convert that to a long long. Take 7 bytes out of that. Then continue.


    printf("What should the second number be? (in decimal)\n");
    scanf("%lld", &b);
    // now storing it in an array of size 16 
    char digitsa[16];
    char digitsb[16];
    char digitsc[16];
    // initialising
    for (int i=0; i<16; i++){
        digitsa[i] = 0;
        digitsb[i] = 0;
        digitsc[i] = 0;
    }
    for (int i=0; i<8; i++){
        digitsa[i] = a>>(8*i) & 255;
        digitsb[i] = b>>(8*i) & 255;
    }
    // Multiplying
    int aid=0;
    int i=0;
    for (int k=0; k<=14; k++){
        for (int j=0; j<=k; j++){
            i = k-j;
            aid = digitsc[k];
            aid += digitsa[i]*digitsb[j];
            digitsc[k] = aid & 255;
            digitsc[k+1] += aid>>8;
        }
    }
    for (int k=0; k<= 14; k++){
        c += ((unsigned long long)digitsc[k])>>(k*8);
    }
    printf("The product is\n%lld\n",c);
    return 0;
}
///////////////////////////////////
// From here on, parallel FFT
///////////////////////////////////

