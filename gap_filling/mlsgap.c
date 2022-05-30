#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <time.h>
#include <errno.h>
#include <float.h>


#define N 184
#define MAX_S 3
#define tau 33

int GaussJordanInv( int ln, double AA[MAX_S][MAX_S] ) ;
double fillgapls(double bfill[N*2],double afill[N*2],double intertime,int index,int miss);
int findgap(double a[N],int i);
double power(double t,int i);

int main()
{
double a[N]={0};            //time variable
double b[N]={0};            //observation vector
double intertime[N]={0};
double intervalue[N]={0};
FILE *fpRead = fopen("datagap4_4.txt","r");//open data file
int num = 1;
double m =0;
double bfill[2*N]={0};
double afill[2*N]={0};

 while (fgetc(fpRead) != EOF)
 {     num ++;
       fscanf(fpRead,"%lf",&m); 
 }
 fclose(fpRead);
int n = num/2;
//printf("there are %d observations\n",n);


FILE *fp = fopen("datagap4_4.txt","r");//prepare for the time vector & observation vector
//printf("Here shows our data\n");
for(int i =0;i < num;i++)
 {     if(i%2==0)
        {
        fscanf(fp,"%lf",&a[i/2]); //input time variable to a[]
        //printf("%lf,",a[i/2]);
        }
      
    else
         {
         fscanf(fp,"%lf",&b[(i-1)/2]);//input measurements to b[] 
        // printf("%lf\n",b[(i-1)/2]);       
         }
      
 }
fclose(fp);
int miss = a[N-1]-N;
for(int i=0;i<miss;i++)
{
intertime[i]= findgap(a,i);
//printf("insert to %lf day\n",intertime[i]);
}
//intervalue[i]= fillgap(a,b,intertime[i]);

for(int i=0;i<N;i++)
{
  if (i<intertime[0]-1)
{bfill[i]=b[i];
 afill[i]=a[i];
}
  else
{bfill[i+miss]=b[i];
 afill[i+miss]=a[i];
}
}
for(int i=0;i<miss;i++)
{ 
bfill[(int)intertime[i]-1]=0;
afill[(int)intertime[i]-1]=0;
}

for(int i=0;i<N+miss;i++)
{
//printf("%d,%lf\n",i,bfill[i]);

}

for(int i=0;i<miss;i++)
{int index=0;
int k = 0;
  if (i%2 == 0)
        { k = i/2;}
    else
      { k = miss - (i+1)/2;}
    index = (int)intertime[k]-1;

  intervalue[k]= fillgapls(bfill,afill,intertime[k],index,miss);
  
  bfill[index]=intervalue[k];
  afill[index]=intertime[k];
}

for(int i=0;i<miss;i++)
{printf("insert to %lf day, %lf\n",intertime[i],intervalue[i]);
}
return 1;
}


double fillgapls(double bfill[N*2],double afill[N*2],double intertime,int index,int miss)
{
double a[tau] = {0};
double b[tau] = {0};
double A[tau][MAX_S]={0};     //design matrix
double AT[MAX_S][tau]={0};    //transposed design matrix
double AA[MAX_S][MAX_S]={0};//normal matrix 
double Al[MAX_S]={0};       //right hand side of normal system
double X[MAX_S]={0};        //unknown parameter contains 'k'&'b'
int size = 0;
double L_hat=0;
int ret=0;
for(int j=0;j<2*miss+1;j++)
{ 
if(bfill[index-miss+j]!=0)
{
b[size] =bfill[index-miss+j];
a[size] =afill[index-miss+j];
size=size+1;
}

//printf("%d,%lf,%lf\n",size,a[size-1],b[size-1]);
}
//printf("Our Design Matrix\n");
for(int i =0;i < size;i++)//calculate design matrix A
 {     for(int j=0;j<MAX_S;j++)
       {A[i][j]=power(a[i],j);//input the first column of design matrix  
       }
       //printf("A[%d][0]=%lf,A[%d][1]=%lf,A[%d][2]=%lf\n",i,A[i][0],i,A[i][1],i,A[i][2]);
 }
for(int j =0;j < MAX_S;j++)//transpose the design matrix to get AT
{ for(int i =0;i < size+1;i++)
      { AT[j][i]=A[i][j];         
      }     
}
//printf("Our normal Matrix\n");
for(int i =0;i < MAX_S;i++)//calculate the normal matrix
{for(int j =0;j < MAX_S;j++)
       { for(int k =0;k < size;k++)
          { AA[i][j]=AA[i][j]+AT[i][k]*A[k][j];
          }
          //printf("%lf ",AA[i][j]);
       }
          //printf("\n");
}
//printf("Our inversed normal Matrix\n");
ret = GaussJordanInv( MAX_S, AA );//inverse of the normal matrix
for(int i=0; i<MAX_S; i++) 
{ for (int j=0; j<MAX_S; j++ )
      {
//printf("%lf   ", AA[i][j] );
    //printf("\n");
      }
  }   
//printf("right hand side of normal equation\n");
for(int i =0;i < MAX_S;i++)//calculate the right hand side of normal equation
{ for(int k =0;k < size+1;k++)
          { Al[i]=Al[i]+AT[i][k]*b[k];
          }
          //printf("%lf\n",Al[i]);
}
//printf("Our unknown paramaters\n");
for(int i =0;i < MAX_S;i++)//calculate the unknown parameter
{for(int j =0;j < MAX_S;j++)
          {X[i]=X[i]+AA[i][j]*Al[j];
          }
}
//printf("a = %lf\nb = %lf\nc = %lf\n",X[0],X[1],X[2]);
//printf("the parabola u = %lf + %lf(t-t0) + %lf(t-t0)*(t-t0) is the approximation of the data aquaired\n",X[0],X[1],X[2]);
 for(int k =0;k < MAX_S;k++)
          { L_hat=L_hat+power(intertime,k)*X[k];
          }          
return(L_hat);
//printf("wma(%d)%lf\n",i,wma[i]);

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

int findgap(double a[N],int i)
{double mt[N] ={0};
int td = 0;
int index =0;
for(int j=1;j<N;j++) 
{td=0;
   if(a[j]-a[j-1]>1)
     {td=  a[j]-a[j-1]-1; 
       for(int k=0;k<td;k++)
       { mt[index]=a[j-1]+k+1;
         index = index+1;
       }
      }
}
return mt[i];
}

double power(double t,int i)
{
double tt=1;
if(i>0)
{
for(int j=0;j<i;j++)
{
tt=tt*t;
}
}

return tt;
}
