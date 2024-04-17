#include <stdlib.h>
#include <stdio.h>

/*
    chread.c
    ==============================

    -Extracts correlator values for a particular value.

*/

void main(){
    
    char path[255];
    char buff[255];
    FILE *fp;
    char *strptr;
            
    int Nt;
    int ch_no;
    
    scanf("%d",&Nt);
    scanf("%d",&ch_no);

    int nmin = ch_no*(Nt+1);    //Start of correlator
    int nmax = nmin+Nt;         //End of correlator

    scanf("%s",&path);

    if ((fp = fopen(path,"r")) == NULL){
        printf("Can't find file!");
        exit(1);
    }

    double G[Nt];

    //Iterate through lines and extract correct channel
    int counter = 0;
    while (fgets(buff,255, (FILE*)fp)) {
        if (counter>nmin && counter<=nmax){
            int i = counter-(nmin+1);
            G[i] = strtold(buff, &strptr);      //Interpret string as float
        }
        counter++;
    }
    fclose(fp);
   
    for (int i=0; i<Nt; i++){
        printf("%.17g\n",G[i]);
    } 
}
