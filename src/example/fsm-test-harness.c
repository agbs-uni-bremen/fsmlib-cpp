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



/**
 *  Execute one test case
 *  @param tcId Test case idntifier
 *  @param line Null-terminated string containing one
 *         test case formatted as
 *         (x1/y1).(x2/y2)...(xn/yn),
 *         where xi is the ith input to the SUT and
 *         yi is the ith expected output.
 *
 *
 */
void executeTestCase(const char* tcId, char* line) {
    
    char* p = line;
    char* x = 0;
    char* y = 0;
    printf("%s",tcId);
    
    while ( *p ) {
        
        if ( p > line ) printf(".");
        
        // After successful completion,
        // x points to the next input string, and
        // y points to the next expected output, also
        // represented as string.
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
    
    // Process test suite file line by line.
    // Every line contains exactly one test case,
    // specified in format
    //  (x1/y1).(x2/y2)...(xn/yn)
    // where xi is the ith input to the SUT and
    // yi is the ith expected output.
    while ( fgets(line,lineSize,f) ) {
        
        size_t len = strlen(line);
        
        if ( len > 1 ) {
            // Replace newline by null character
            line[len-1] = 0;
            char tcId[100];
            *tcId = 0;
            sprintf(tcId,"TC-%d: ",++tcNum);
            // Execute one test case, specified in the current line
            executeTestCase(tcId,line);
            // Reset SUT via wrapper function, since this
            // test case has been completely executed
            sut_reset();
        }
        
    }
    
}



/**
 *  The test harness is invoked with
 *  @param Name of test suite text file.
 *  The test suite file is expected to be in the format
 *  of test suites created by program fsm-generator.
 *
 *  The test harness is linked to the SUT wrapper and the
 *  SUT code, as described in
 * http://www.informatik.uni-bremen.de/agbs/jp/papers/test-automation-huang-peleska.pdf
 *  Appendix B.5.
 */
int main(int argc, char** argv) {
    
    if ( argc < 2 ) {
        fprintf(stderr,"Missing file name of test suite file - exit.\n");
        exit(1);
    }
    
    // Initialise the SUT via wrapper function sut_init()
    sut_init();
    
    // Execute all test cases contained in the
    // test suite file.
    executeTestCases(argv[1]);
    
    exit(0);
    
}
