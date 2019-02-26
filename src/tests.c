#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

double mean(double *val, int len);
double var(double *val, double mean, int len);
double ws(double var1, double var2, int n1, int n2);
double welch(double mean1, double mean2, int n1, int n2, double var1, double var2);
void swap(int *a, int *b);
void randomize(int arr[], int n);
int cmpGrp(const void *x, const void *y);
int cmp(const void *x, const void *y);


//Global vars
double *array; //Needed for sorting
int *arrayGrp;

//without WY p correction
void tTest_unpaired_unequalVar_simple(double *data, int *lenRows, int *grp, int *lenGrp, double *df, double *statistic) {
    //store teststatistic per row
    double U[lenRows[0]];
    int i;
    int j;

    //Get Indices for groups
    int pos = 0; // First index of second group
    for (i=1; i < lenGrp[0]; ++i) {
	//printf("VAL %f\n", data[i]);
	if (grp[i-1] != grp[i]) {
	    pos = i;
	}
    }

    //Calculate actual testtstiatic
    //double actTest[lenRows[0]]; 
    //double actTestAbs[lenRows[0]]; //absolute value
    //int index[lenRows[0]];
    for (i=0; i < lenRows[0]; i++) {
	//Pointer to first element of each row
	double *start = &data[i*lenGrp[0]]; 

	//calc means and vars
	double m1 = mean(start, pos); 
	double m2 = mean(start+pos, lenGrp[0]-pos); 
	//printf("MEAN 1 %f, MEAN 2 %f \n", m1, m2);
	double v1 = var(start, m1, pos);
	double v2 = var(start+pos, m2, lenGrp[0]-pos);
	//printf("VAR 1 %f, VAR 2 %f \n", v1, v2);

	//Welch's test
	double stat = welch(m1, m2, pos, lenGrp[0]-pos, v1, v2);
	//printf("Teststat: %f\n", stat);

	//Degrees of freedom
	double dg = ws(v1, v2, pos, lenGrp[0]-pos);
	//printf("Degrees of Freedom: %f\n", dg);

	//actTest[i] = stat; 
	//actTestAbs[i] = fabsl(stat); 
	//index[i] = i;
	
	statistic[i] = stat;
	df[i] = dg;
    }
}


void tTest_unpaired_unequalVar(double *data, int *lenRows, int *grp, int *lenGrp, int *B, double *pval) {
   //Grouping
   printf("LenGrp: %d\n", lenGrp[0]);
   printf("LenRows: %d\n", lenRows[0]);
   printf("Perm: %d\n", B[0]);
   
   //Init random gen
   srand(time(NULL));

   //Store all testatistics
   double U[B[0]*lenRows[0]];

   int i;
   int j;
   int b;

   //Get Indices for groups
   int pos = 0; // First index of second group
   for (i=1; i < lenGrp[0]; ++i) {
	//printf("VAL %f\n", data[i]);
	if (grp[i-1] != grp[i]) {
	  pos = i;
	}
   }

   //Calculate actual testtstiatic
   double actTest[lenRows[0]]; 
   double actTestAbs[lenRows[0]]; //absolute value
   int index[lenRows[0]];
   for (i=0; i < lenRows[0]; i++) {
     //Pointer to first element of each row
     double *start = &data[i*lenGrp[0]]; 

     //calc means and vars
     double m1 = mean(start, pos); 
     double m2 = mean(start+pos, lenGrp[0]-pos); 
     //printf("MEAN 1 %f, MEAN 2 %f \n", m1, m2);
     double v1 = var(start, m1, pos);
     double v2 = var(start+pos, m2, lenGrp[0]-pos);
     //printf("VAR 1 %f, VAR 2 %f \n", v1, v2);

     //Welch's test
     double stat = welch(m1, m2, pos, lenGrp[0]-pos, v1, v2);
     //printf("Teststat: %f\n", stat);

     //Degrees of freedom
     //double df = ws(v1, v2, pos, lenGrp[0]-pos);
     //printf("Degrees of Freedom: %f\n", df);
     
     actTest[i] = stat; 
     actTestAbs[i] = fabsl(stat); 
     index[i] = i;
   }

   // Reorder data: decreasing absolute actua test statistic vals
   array = actTestAbs;
   qsort(index, lenRows[0], sizeof(*index), cmp);
   double actTestAbsOrd[lenRows[0]]; //ordered teststat
   double dataOrd[lenRows[0]*lenGrp[0]]; //ordered data
   for  (i=0; i < lenRows[0]; i++) {
     actTestAbsOrd[i] = actTestAbs[index[i]];
     //Pointer to first element of each row, original data
     double *start = &data[i*lenGrp[0]];
     memcpy(&dataOrd[index[i]*lenGrp[0]], start, sizeof(double)*lenGrp[0]);  
   }

   //--------------------------------------------------------
   for (b=0; b < B[0]; b++) {
       // Permute group assignment
       int grpNew[lenGrp[0]];
       memcpy(&grpNew, grp, sizeof(int)*lenGrp[0]);
       randomize(grpNew, lenGrp[0]);
       int grpIndex[lenGrp[0]];
       for (i=0; i < lenGrp[0]; i++) {
         //printf("%d %d\n", grp[i], grpNew[i]);
         grpIndex[i] = i;
       }
       arrayGrp = grpNew;
       qsort(grpIndex, lenGrp[0], sizeof(*grpIndex), cmpGrp);
       
       // Sort data
       double dataNew[lenRows[0]*lenGrp[0]];
       for (i=0; i < lenGrp[0]; i++) {
         for (j=0; j < lenRows[0]; j++) {
           //Pointer to first element of each column
           double *from = &data[grpIndex[i]]+j*lenGrp[0];
           double *to = &dataNew[i] + j*lenGrp[0];
           memcpy(to, from, sizeof(double));
         }
       }
       for (i=0; i < lenGrp[0]; i++) {
         double *start = &dataNew[i];
         printf("%d->%d  ", grpIndex[i], arrayGrp[grpIndex[i]]);
       }
       printf("\n");
    
       //Calculate test statics 
       double newTestStat[lenRows[0]];
       for (i=0; i < lenRows[0]; i++) {
         double *start = &dataNew[i*lenGrp[0]];
         double m1 = mean(start, pos);
         double m2 = mean(start+pos, lenGrp[0]-pos);
         double v1 = var(start, m1, pos);
         double v2 = var(start+pos, m2, lenGrp[0]-pos);
         double stat = welch(m1, m2, pos, lenGrp[0]-pos, v1, v2);
         newTestStat[i] = stat;
       }
    
       //Order 
       double u[lenRows[0]];
       for (i=0; i < lenRows[0]; i++) {
         u[i] = fabsl(newTestStat[lenRows[0]]);
       }
       
       for (i=lenRows[0]-1; i>=0; i--) {
         if (u[i+1] > fabsl(newTestStat[i])) {
           u[i] = u[i+1];
         } else {
           u[i] = fabsl(newTestStat[i]);
         }
         //printf("%f\n", u[i]);
       }

       printf("%d %d \n ", (int)sizeof(*U), B[0]*lenRows[0]);

       double *start = &U[b*lenRows[0]];
       memcpy(start, u, sizeof(double)*lenRows[0]); 
   }

   // ---------------------------------------------
   // Get results from different permutations
   double pAdj[lenRows[0]];
   for (i=0; i < lenRows[0]; i++) {
     int count = 0;
     for (b=0; b < B[0]; ++b) {
       //double *start = U[b];
       double *start = &U[b*lenRows[0]];
       if ((start+i)[0] >= actTestAbsOrd[i]) {
         count++;
       }
       //printf("i %d, b %d, val: %f\n", i, b, (start+i)[0]);
     }
     pAdj[i] = (double)count/B[0];
     //printf("pAdj %f\n", pAdj[i]);
   }

   //Assure monotony
   for (i = 2; i < lenRows[0]; i++) {
     if (pAdj[i-1] > pAdj[i]) {
       pAdj[i] = pAdj[i-1]; 
     }
   }

   //Reorder 
   for (i=0; i < lenRows[0]; i++) {
     pval[index[i]] = pAdj[i];
   }

}

void perm_tTest_unpaired_unequalVar(double *data, int *len, int *grp, int *B, double *pval) {
   //Grouping
   //printf("Len: %d\n", len[0]);
   //printf("Perm: %d\n", B[0]);

   //Init random gen
   srand(time(NULL));

   //Store all teststistics
   double U[B[0]];

   //Count vars
   int i;
   int j;
   int b;

   //Get Indices for groups
   int pos = 0; // First index of second group
   for (i=1; i < len[0]; ++i) {
        //printf("VAL %f\n", data[i]);
        if (grp[i-1] != grp[i]) {
          pos = i;
        }
   }

   //Calculate actual test statistics
   double *start = data;
   //calc means and vars
   double m1 = mean(start, pos);
   double m2 = mean(start+pos, len[0]-pos);
   //printf("MEAN 1 %f, MEAN 2 %f \n", m1, m2);
   double v1 = var(start, m1, pos);
   double v2 = var(start+pos, m2, len[0]-pos);
   //printf("VAR 1 %f, VAR 2 %f \n", v1, v2);
   //Welch's test
   double actStat = welch(m1, m2, pos, len[0]-pos, v1, v2);
   //printf("Teststat: %f\n", actStat);
   //-------------------------------
   for (b=0; b < B[0]; b++) {
       // Permute group assignment
       int grpNew[len[0]];
       memcpy(&grpNew, grp, sizeof(int)*len[0]);
       randomize(grpNew, len[0]);
       int grpIndex[len[0]];
       for (i=0; i < len[0]; i++) {
         //printf("%d %d\n", grp[i], grpNew[i]);
         grpIndex[i] = i;
       }
       arrayGrp = grpNew;
       qsort(grpIndex, len[0], sizeof(*grpIndex), cmpGrp);

       // Sort data
       double dataNew[len[0]];
       for (i=0; i < len[0]; i++) {
           double *from = &data[grpIndex[i]];
           double *to = &dataNew[i];
           memcpy(to, from, sizeof(double));
       }
       for (i=0; i < len[0]; i++) {
         double *start = &dataNew[i];
         //printf("%f ", dataNew[i]);
         //printf("%d->%d  ", grpIndex[i], arrayGrp[grpIndex[i]]);
       }
       //printf("\n");

       //Calculate test statics
       double *start = dataNew;
       double m1 = mean(start, pos);
       double m2 = mean(start+pos, len[0]-pos);
       double v1 = var(start, m1, pos);
       double v2 = var(start+pos, m2, len[0]-pos);
       double stat = welch(m1, m2, pos, len[0]-pos, v1, v2);
       //printf("m1 %f, m2 %f v1 %f v2 %f\n", m1, m2, v1, v2);

       double *st = &U[b];
       memcpy(st, &stat, sizeof(double));
   }

   //----------------------------
   int count = 0;
   for (b=0; b<B[0]; b++) {
     //two sided 
     //printf("B[b] %f, actStat %f\n", U[b], actStat);
     if (fabsl(U[b]) >= fabsl(actStat)) {
       count++;
     }
   }

   //calc p value
   pval[0] = (double)count/B[0];
   //printf("p-value: %f; stat: %f \n",  pval[0], actStat);
}




/*
 * Welch-Satterthwaite approximation
 * to get the degrees of freedom
 */
double ws(double var1, double var2,
        int n1, int n2) {
  double nom = pow(var1 + var2, 2);
  double denom = pow(var1, 2)/((n1-1)) + pow(var2, 2)/((n2-1));

  return (nom/denom);
}

/*
 * Calculates the Welch's t teststatistic
 * for unequal variances
 */
double welch(double mean1, double mean2,
                int n1, int n2,
                double var1, double var2) {
  return ((mean1-mean2)/sqrt(var1/n1 + var2/n2));
}

/*
 * Calculates the variance
 */
double var(double *val, double mean, int len) {
  double res = 0.0;
  int i;
  for (i=0; i<len; i++) {
    res += pow(val[i]-mean, 2);
  }
  return (res/(len-1));
}

/*
 *  Calculates mean for a given number
 *  of doubles
 */
double mean(double *val, int len) {
  double ret = 0.0;
  int i;
  for (i=0; i<len; i++) {
    ret += val[i];
  }

  return (ret/len);
}

int cmpGrp(const void *x, const void *y) {
  int xx = *(int*)x;
  int yy = *(int*)y;
  return arrayGrp[xx] < arrayGrp[yy] ? -1 : arrayGrp[xx] > arrayGrp[yy];
}
int cmp(const void *x, const void *y) {
  int xx = *(int*)x;
  int yy = *(int*)y;
  return array[xx] > array[yy] ? -1 : array[xx] < array[yy];
}
void swap(int *a, int *b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}

void randomize(int arr[], int n) {
    int i;
    for(i = n-1; i > 0; i--) {
        int j = rand() % (i+1);
        swap(&arr[i], &arr[j]);
    }
}

