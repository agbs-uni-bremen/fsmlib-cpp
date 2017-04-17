#include <stdlib.h>
#include <string.h>
#include <stdio.h>


extern void sut_init();

extern void sut_reset();

extern const char* sut(const char* input);



void getNextIO(char** p, char** x, char** y) {
    
    *x = NULL;
    *y = NULL;
    
    char* aux;
    
    while ( **p != 0 && **p != '(' ) (*p)++;
    if ( **p == 0 ) return;
    
    (*p)++;
    aux = strchr(*p,'/');
    if ( aux == NULL ) return;
    
    *x = *p;
    *aux = 0;
    
    *p = aux+1;
    aux = strchr(*p,')');
    if ( aux == NULL ) return;
    
    *y = *p;
    *aux = 0;
    
    *p = aux + 1;
    
}



void executeTestCase(const char* tcId, char* line) {
    
    char* p = line;
    char* x = 0;
    char* y = 0;
    printf("%s",tcId);
    
    while ( *p ) {
        
        if ( p > line ) printf(".");
        
        getNextIO(&p,&x,&y);
        
        if ( x != NULL && y != NULL ) {
            const char* r = sut(x);
            if ( strcmp(r,y) != 0 ) {
                printf(" after input %s: expected %s - observed %s: FAIL\n",
                       x,y,r);
                return;
            }
            else {
                printf("(%s,%s)",x,r);
            }
        }
        
    }
    
    printf(" PASS\n");
    
}


void executeTestCases(const char* fname) {
    
    const int lineSize = 100000;
    char* line = (char*)calloc(lineSize,1);
    FILE* f = fopen(fname,"r");
    if ( f == NULL ) {
        fprintf(stderr,"Could not open file %s - exit.\n",fname);
        exit(1);
    }
    
    int tcNum = 0;
    while ( fgets(line,lineSize,f) ) {
        
        size_t len = strlen(line);
        
        // Replace newline by null character
        if ( len > 1 ) {
            line[len-1] = 0;
            char tcId[100];
            *tcId = 0;
            sprintf(tcId,"TC-%d: ",++tcNum);
            executeTestCase(tcId,line);
            sut_reset();
        }
        
    }
    
}




int main(int argc, char** argv) {
    
    if ( argc < 2 ) {
        fprintf(stderr,"Missing file name of test suite file - exit.\n");
        exit(1);
    }
    
    sut_init();
    
    executeTestCases(argv[1]);
    
    exit(0);
    
}
