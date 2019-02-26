#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

int findeNachbar(int *m, int *visited, int start, int invStart, int len,  int size, int found, int clID, int star);
int findCluster(int *m, int *visited, int len, int size, int id, int star);

void startAdj(int *matrix, int *len, int *n, int *size, int *erg, int *maxTry, int *star) {
    printf("Matrix: %d x %d \n", len[0], len[0]);
    printf("Pos [0,0]: %d\n", matrix[0]);

    //adj matrix m 
    int m[len[0]][len[0]];

    //numbers of Connections *2
    int nCon = 0;

    //fill adjacency matrix
    for (int i=0; i < len[0]; i++) {
        for (int j = 0; j < len[0]; j++) {
            m[i][j] = matrix[i*len[0]+j];
            if (m[i][j] == 0) {
                //printf(".");
            } else {
                //printf("|");
                nCon++;
            }
        }
        //printf("\n");
    }
    printf("No of connections: %d\n", nCon);

    //Visited-Matrix
    //int visited[len[0]*len[0]];
    int *visited = malloc(len[0]*len[0]*sizeof(visited));
    memset(visited, 0, len[0]*len[0]*sizeof(visited));

    // Finde n zusammenhängende Cluster
    // 0 is used to mark not-visited cells
    for (int i=1; i <= n[0]; i++) {
        printf("   Cluster: %d\n", i);
        int versuche = maxTry[0]; 
	int success;
        do {
            versuche--;
	    success = findCluster(m, visited, len[0], size[0], i, star[0]);
        } while (!success && versuche > 0);
	//printf("SUCCESS CLUSTER: %d\n", success);
    }

    
    memcpy(erg, visited, len[0]*len[0]*sizeof(int));
    free(visited);
}


//Findet zusammenhängede Cluster der Größe size
//Suche in adj, teste, ob in visited schon besucht,
//markiere als besucht in visited, return bei Erfolg
int findCluster(int *m, int *visited, int len, int size, int id, int star) {
    //Finde zufälligen startpunkt
    int maxIndex = len*len;
    int start = rand() % (maxIndex + 1);
    int invStart; 
    printf("   -> Search Startindex: %d ", start);
    for (;start < maxIndex; start++) {
        if (m[start] == 1 && visited[start] == 0) {
            break;
        }
    }
    invStart = len*(start%len)+start/len;
    printf("Found at: %d / %d. ", start, invStart);
    //Convert to row/index
    int rowIndex=start % len;
    int colIndex=start / len;
    printf("Converted to [%d,%d]\n.", rowIndex, colIndex);

    //Erzeuge Arbeitskopie von visited
    int *tmpVisited = malloc(maxIndex*sizeof(tmpVisited));
    memcpy(tmpVisited, visited, maxIndex*sizeof(tmpVisited));

    //Markiere startposition
    //printf("START: %d, INVST: %d\n", start, invStart);
    tmpVisited[start] = id;
    tmpVisited[invStart] = id;

    // Finde nachbar
    int success = findeNachbar(m, tmpVisited, start, invStart, len, size, 1, id, star);
    if (success == 1) {
	printf("Copy back\n");
	memcpy(visited, tmpVisited, maxIndex*sizeof(tmpVisited));
    }
    free(tmpVisited);

    int count = 0;
    int x,y;
    for (int i=0;i < maxIndex; i++) {
	if (visited[i] == id) {
	    count++;
	    ///
	    x = i/len;
	    y = i%len;
	    printf(" CL %d: %d [%d,%d]\n", id, i,x,y);
	    //printf("v [%d,%d]=>%d\n", i,j,visited[i*j]);
	}
    }
    printf("Counted in Cluster: %d\n", count);

    return(success);
}


int findeNachbar(int *m, int *visited, int start, int invStart, int len, int size, int found, int clID, int star) {
    //printf("SIZE %d, FOUND %d\n", size, found);
    if (size == found) { 
	return(1);
    } 

    // Find neighbor
    int pos = start;
    
    int rowIndex=start % len;
    int colIndex=start / len;
    printf("Starting from: [%d,%d]\n.", rowIndex, colIndex);

    int c, max;

    // Search per row
    c = (int)floor(start / len);
    if (c < 0) { 
	c = 0;
    }
    max = (c+1)*len;
    if (max > len*len) { 
	max = len*len;
    }
    //printf("FROM :%d %d, TO %d\n", c*len, c, max);
    for (int i = c*len; i < max; i++) {
	if (m[i] == 1 && visited[i] == 0) {
	    //Found 
	    pos = i;
	    //printf("FOUND I: i=%d ADJ %d VISITED %d\n", i, m[i], visited[i]);
	    rowIndex=i % len;
	    colIndex=i / len;
	    printf("->  [%d,%d]\n.", rowIndex, colIndex);
	    break;
	}
    }
    
    if ( pos == start && star == 0) {
	c = (int)floor(invStart / len);
	if (c < 0) { 
	    c = 0;
	}
	max = (c+1)*len;
	if (max > len*len) { 
	    max = len*len;
	}

	// Search per col
	for (int i = c*len; i < max; i++) {
	    if (m[i] == 1 && visited[i] == 0) {
		//Found 
		pos = i;
		//printf("FOUND I: i=%d ADJ %d VISITED %d\n", i, m[i], visited[i]);
		break;
	    }
	}
    } 

    //Did we find a new neighbor?
    if (pos != start) {
	int nachbar, invNachbar;
	nachbar = pos;
	invNachbar = len*(nachbar%len)+nachbar/len;
	visited[nachbar] = clID;
	visited[invNachbar] = clID;
	printf("\tREGISTER: [%d, %d], [%d, %d]\n", nachbar%len, nachbar / len, invNachbar%len, invNachbar / len);

	found++;
	//printf("Currently found: %d\n", found);
	return(findeNachbar(m, visited, nachbar, invNachbar, len, size, found, clID, star));
    } 
    return(0);
}
