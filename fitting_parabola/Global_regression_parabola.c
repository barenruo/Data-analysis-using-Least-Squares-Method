#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <time.h>
#include <errno.h>
#include <float.h>


#define MAX_S   3
#define N 101
#define E 2.7182818284590452354
#define PI 3.1415926

int GaussJordanInv( int ln, double AA[MAX_S][MAX_S] );

/**--------------------------------------------------------------+
                                 MAIN()
    File LIST: "section4_2.txt.txt" 
    GLOBAL   : MAX_S (max. size of the normal matrix)
                    N     (the number of observation)
    CALLS: GausJordanInv()
    Created: 2020-12-03; Renruo Ba
 +---------------------------------------------------------------*/
int main()
{
double a[N]={0};            //time variable
double b[N]={0};            //observation vector
double A[N][MAX_S]={0};     //design matrix
double AT[MAX_S][N]={0};    //transposed design matrix
double AA[MAX_S][MAX_S]={0};//normal matrix 
double Al[MAX_S]={0};       //right hand side of normal system
double X[MAX_S]={0};        //unknown parameter contains 'k'&'b'
double trend[N]={0};
double detrend[N]={0};
double Fnreal[N]={0}; 
double Fnimagine[N]={0}; // real and imaginary parts of F(n)
double Am[N];           // amplitude of sin wave at kth frequency
double Ph[N];          //Phase of the sin wave with A to be amplitude 
double largestamplitude= 0;  //largest amplitude from FFT
double secondlargeamplitude= 0;  //second largest amplitude from FFT
int largest = 0;
int secondlarge = 0;
int ret = 0;


FILE *fpRead = fopen("section4_2.txt","r");//open data file
int num = 1;
double m =0;
 while (fgetc(fpRead) != EOF)
 {     num ++;
       fscanf(fpRead,"%lf",&m); 
 }
 fclose(fpRead);
int n = num/2;
printf("there are %d observations\n",n);


FILE *fp = fopen("section4_2.txt","r");//prepare for the time vector & observation vector
printf("Here shows our data\n");
for(int i =0;i < num;i++)
 {     if(i%2==0)
        {
        fscanf(fp,"%lf",&a[i/2]); //input time variable to a[]
        printf("%lf,",a[i/2]);
        }
      
    else
         {
         fscanf(fp,"%lf",&b[(i-1)/2]);//input measurements to b[] 
         printf("%lf\n",b[(i-1)/2]);       
         }
      
 }
fclose(fp);

/*--------------------------------------------------------------
     function model:  u = a + b * (t-t0) + c * (t-t0)*(t-t0), 
     Unknown parameter: X = [a;b;c]
     observation equation: [u]=[A][X]
     design matrix: A=[1 (t-t0) (t-t0)*(t-t0)]
---------------------------------------------------------------*/
//printf("Our Design Matrix\n");
for(int i =0;i < N;i++)//calculate design matrix A
 {     
       A[i][0]=1;//input the first column of design matrix  
       A[i][1]=a[i]-a[0];//input the second column of design matrix
       A[i][2]=(a[i]-a[0])*(a[i]-a[0]);//input the third column of design matrix
       printf("A[%d][0]=%lf,A[%d][1]=%lf,A[%d][2]=%lf\n",i,A[i][0],i,A[i][1],i,A[i][2]);
 }


for(int j =0;j < MAX_S;j++)//transpose the design matrix to get AT
{
      for(int i =0;i < N;i++)
      {
         AT[j][i]=A[i][j];
         
      }
      
}
printf("Our normal Matrix\n");
for(int i =0;i < MAX_S;i++)//calculate the normal matrix
{
      for(int j =0;j < MAX_S;j++)
       {
          for(int k =0;k < N;k++)
          {
             AA[i][j]=AA[i][j]+AT[i][k]*A[k][j];
             
          }
          printf("%lf ",AA[i][j]);

       }
          printf("\n");
}

printf("Our inversed normal Matrix\n");
ret = GaussJordanInv( MAX_S, AA );//inverse of the normal matrix

  for(int i=0; i<MAX_S; i++) {
    for (int j=0; j<MAX_S; j++ )
      printf("%lf   ", AA[i][j] );
    printf("\n");
  }   

//printf("right hand side of normal equation\n");
for(int i =0;i < MAX_S;i++)//calculate the right hand side of normal equation
{
 
       for(int k =0;k < N;k++)
          {
             Al[i]=Al[i]+AT[i][k]*b[k];
             
          }
          printf("%lf\n",Al[i]);


}

printf("Our unknown paramaters\n");
for(int i =0;i < MAX_S;i++)//calculate the unknown parameter
{
 
       for(int j =0;j < MAX_S;j++)
          {
             X[i]=X[i]+AA[i][j]*Al[j];
             
          }
          
}
printf("a = %lf\nb = %lf\nc = %lf\n",X[0],X[1],X[2]);


printf("the parabola u = %lf + %lf(t-t0) + %lf(t-t0)*(t-t0) is the approximation of the data aquaired\n",X[0],X[1],X[2]);


return 1;

}







int GaussJordanInv( int ln, double AA[MAX_S][MAX_S] ) 
{
int i=0, j=0, k=0, m=0, nr=0, ret=0;
double Piv=0.0e0, tmpV=0.0e0;
int PCol[ln]; // vector for storing of the pivotal columns
  for(i=0; i<ln; i++ ) PCol[i]=0;
  for( nr=0; nr<ln; nr++ ) {
    Piv=0.0E0;  
    for ( j=nr; j<ln; j++ ) {
      if ( Piv >= fabs( AA[nr][j]) ) continue;
      Piv = fabs( AA[nr][j] );
      PCol[nr] = j; // the pivotal column
    }//for(j=rn;j<ln;j++) 
    j = PCol[nr];
    for ( k=0; k<ln; k++ ) {
      tmpV = AA[k][j];
      AA[k][j] = AA[k][nr];
      AA[k][nr] = tmpV;
    }//for(k=0;k<ln;k++)
    tmpV = 1.0E0/AA[nr][nr];
    AA[nr][nr] = 1.0E0;
    for ( j=0; j<ln; j++ )
      AA[nr][j] = tmpV * AA[nr][j];
    for ( k=0; k<ln; k++ ) {
      if ( nr == k ) continue;
      tmpV = AA[k][nr];
      AA[k][nr] = 0.0E0;
      for ( j=0; j<ln; j++ )
        AA[k][j] = AA[k][j] - tmpV * AA[nr][j];
    }// for ( k=0; k<ln; k++ )
  }//for(l=0;l<N;l++)
  for ( i=0; i<ln; i++ ) {
    m = ln - i - 1;
    k = PCol[m];
    for ( j=0; j<ln; j++ ) {
      tmpV = AA[m][j];
      AA[m][j] = AA[k][j];
      AA[k][j] = tmpV;
    }//for(j=0;j<ln;j++)
  } //for(i=0;i<=ln;i++)
  return(ret);
}// GaussJordanInv()
