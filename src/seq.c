#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

//matrix: input 
//sel consens
//n number of consens
void cons(int *matrix, int *len1, int *len2, int *sel, int *n, int *katMax) {
    printf("matrix: %dx%d\n", len1[0], len2[0]);

    int i, j, key;    
    for (i=0; i < len2[0]; i++) {    
	int count[katMax[0]+1];    
	for (j=0; j<katMax[0]+1; j++) {    
	    count[j] = 0;    
	}    

	//Test per column    
	for (j=0; j < len1[0]; j++) {     
	    key = matrix[i*len1[0]+j];    
	    count[key] = count[key]+1;    
	}    

	//Get most freq       
	sel[i] = count[0];    
	n[i] = 0;    
	for (j=1; j<katMax[0]+1; j++) {    
	    if (sel[i] < count[j]) {    
		sel[i] = count[j];    
		n[i] =j;    
	    }    
	}    
    }
}

void mtchSub(int *len, char **pattern, char **string, int *seqN, int *posM) {    
    char *pch = strstr(*string, *pattern);    
    //printf("Pattern: %s\n", pattern[0]);    
                              
    if(pch) {    
        // get first occurence         
        int pos = -1*(*string - pch);    
	posM[0] = pos;
        //printf("found @ %d\n", pos);    
        int i;                             
        //printf("len %d\n", len[0]);      
        for (i=pos; i<len[0]; i++) {                    
            char *c = string[0];        
            c = c+i;    
            if (c[0] == 'A') {    
                seqN[i] = 1;      
            } else if (c[0] =='C') {    
                seqN[i] = 2;    
            } else if (c[0] =='T') {    
                seqN[i] = 4;    
            } else if (c[0] =='G') {    
                seqN[i] = 3;    
            }    
        }        
    }        
}    

void mtchSubRev(int *len, char **pattern, char **string, int *seqN, int *posM) {    
    char *pch = strstr(*string, *pattern);    
    //printf("Pattern: %s\n", pattern[0]);    
                              
    if(pch) {    
        // get first occurence         
        int pos = -1*(*string - pch);    
	posM[0] = pos;
        //printf("found @ %d\n", pos);    
        int i;                             
        //printf("len %d\n", len[0]);      
        for (i=0; i<pos+len[0]; i++) {                    
            char *c = string[0];        
            c = c+i;    
            if (c[0] == 'A') {    
                seqN[i] = 1;      
            } else if (c[0] =='C') {    
                seqN[i] = 2;    
            } else if (c[0] =='T') {    
                seqN[i] = 4;    
            } else if (c[0] =='G') {    
                seqN[i] = 3;    
            }    
        }        
    }        
}    

void mtch(char **pattern, char **string, int *res) {
    char *pch = strstr(*string, *pattern);
    if(pch) { 
	res[0]=1 ;
    }
}

