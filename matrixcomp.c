/* 
   CS288 HOMEWORK 8
   Your program will take in two command-line parameters: n and error
   command: jacobi 5 0.0001
   command: jacobi 10 0.000001
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>

#define N 100
#define MAX_ITER 10000

int jacobi();
void init();
int convergence();
void srand();
void print_vector();
void print_equation();
void final_check();

float a[N][N], b[N];
float x[N], buf[N];
int n;
float error;

int main(int argc, char **argv){
  int n_iter;			/* number of iterations */
  n = atoi(argv[1]);
  error = atof(argv[2]);
  init();		   /* initalize a, x0 and b - DO not change */
  print_equation();
  n_iter = jacobi();
  printf("...Solved after %d iterations\n", n_iter);
  print_equation();
  final_check();
  return 0;
}

int jacobi(){
  int i=0;int j = 0; int k=0;
  float sum;
  while(!convergence(k)){
  //loop through row
  for (i=0; i<n; i++){
  float diag = a[i][i];
  sum = b[i];
  //loop through column
  for(j=0; j<n; j++){
    if(i!=j){//if not a diagnal element in matrix
    sum-= a[i][j]*x[j];//subtract every remainder times X from b (b-Rx=sum)
    } 
	}
  x[i]=(sum/diag);//get new x (Xnew=sum/diagnal)
  printf("X%d=%f ", i, x[i] );
    }//end for
  printf("pass %d\n", k);
  printf("\n");
  k++;
  }//end while
  return k;
}

// returns 1 if converged else 0
int convergence(int iter){
  int i,j,flag=0;// init to false
  float sum;
   for(i=0; i<n; i++){//loop through row
  float k=0;
    for(j=0; j<n; j++){//loop through column
       k+=(a[i][j] * x[j]);//multiply every A in a row by every X in the vector and add them all up [a1  a2 a3 a4] *[x1 x2 x3 x4 ] = b
	}
  float temp=(k - b[i]);
  if(temp<0)//if temp is negative make positive (correctional step)
  temp = temp*(-1);
  sum+=(temp);//That way sum is always positive and error cannot cancel other error
  }//end for loop
  sum = sum/n;
  if(!(sum<error)){
   flag=0;
   //break;
    }
  else if(sum<error){//converged!
   flag=1;
    }
  return flag;
}

// Try not to change this. Use it as is.
void init(char **argv){
  int i,j,k,flag=0;
  float sum;
  int seed = time(0) % 100;	/* seconds since 1/1/1970 */

  srand(seed);
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      a[i][j] = rand() & 0x7;
      if ((rand() % 2)==1) a[i][j] = -a[i][j]; //if odd make negative
    }//for each row
    sum = 0;
    for (j=0;j<n;j++) if(i!=j) sum = sum + abs(a[i][j]);
    if (a[i][i] < sum) a[i][i] = sum + a[i][i];
  }

  for (i=0;i<n;i++) x[i]=1;

  srand(seed);
  for (i=0;i<n;i++){
    b[i]=rand() & 0x7;
    if ((rand() & 0x1)) b[i] = -b[i];
  }

  //print_equation();

}

void print_equation(){
  int i,j;

  printf("A*x=b\n");
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) printf("%2d ",(int)(a[i][j]));
    printf(" * x%d = %d\n",i,(int)(b[i]));
  }
  printf("\n");
  print_vector(x);
}

void print_vector(float *l){
  int i;
  for (i=0; i<n; i++) printf("%.6f ",l[i]);
  printf("\n");
}

void final_check(){
printf("Final check...\n");
int i, j;
float val;
  for (i=0;i<n;i++){
  val=0;
    for (j=0; j<n; j++){
	val+= (a[i][j] * x[j]);
    }
  printf("val:%f vs b: %f\n", val, b[i]);
  }
}

// end of file

