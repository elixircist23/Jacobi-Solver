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

float a[N][N], b[N];
float x[N], buf[N];
int n;
float error;

int main(int argc, char **argv){
  int n_iter;
  n = atoi(argv[1]);
  error = atof(argv[2]);

  init();
  n_iter = jacobi();
  print_vector(buf);
  printf("iterations: %i\n", n_iter);

  return 0;
}

int jacobi(){
  int i,j,k;
  int iter = 0;
  float sum;

  while(!convergence(iter) && iter < MAX_ITER){
    for(int i = 0; i < n; i++){
      sum = 0;
      for(int j = 0; j < n; j++){
        if(j != i){
          sum = sum + (a[i][j] * x[j]);
        }
      }
      buf[i] = (b[i] - sum) / (a[i][i]);
    }
    for(k = 0; k < n; k++)
      x[k] = buf[k];
    iter++;
  }
  return k;
}

int convergence(int iter){
  int i,j,flag=1;
  float k = 0;
  for(i = 0; i < n; i++){
    k = 0;
    for(j = 0; j < n; j++){
      k = k + (a[i][j] * x[j]);
    }
    if((k - b[i]) > error)
      return 0;
  }

  return flag;
}

void init(char **argv){
  int i,j,k,flag=0;
  float sum;
  int seed = time(0) % 100;	/* seconds since 1/1/1970 */

  srand(seed);
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
      a[i][j] = rand() & 0x7;
      if (rand() & 0x1) a[i][j] = -a[i][j];
    }
    sum = 0;
    for (j=0;j<n;j++) if(i!=j) sum = sum + abs(a[i][j]);
    if (a[i][i] < sum) a[i][i] = sum + a[i][i];
  }

  for (i=0;i<n;i++) x[i]=1;

  srand(seed);
  for (i=0;i<n;i++){
    b[i]=rand() & 0x7;
    if (rand() & 0x1) b[i] = -b[i];
  }

  print_equation();

}

void print_equation(){
  int i,j;

  printf("A*x=b\n");
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) printf("%2d ",(int)(a[i][j]));
    printf(" * x%d = %d\n",i,(int)(b[i]));
  }
  printf("\n");
}

void print_vector(float *l){
  int i;
  for (i=0; i<n; i++) printf("%.6f ",l[i]);
  printf("\n");
}
