#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <time.h>
#include <errno.h>
#include <float.h>


#define N 200
#define tau 21

int main()
{
double a[N]={0};            //time variable
double b[N]={0};            //observation vector
double wma[N-tau+1]={0};
double window[tau]={0};
double weight[tau]={0};
double sum=0;
FILE *fpRead = fopen("section4_3.txt","r");//open data file
int num = 1;
double m =0;
 while (fgetc(fpRead) != EOF)
 {     num ++;
       fscanf(fpRead,"%lf",&m); 
 }
 fclose(fpRead);
int n = num/2;
//printf("there are %d observations\n",n);


FILE *fp = fopen("section4_3.txt","r");//prepare for the time vector & observation vector
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
for(int j=0;j<tau;j++)
{
weight[j] = j+1;
 if(j>(tau-1)/2)
{weight[j]=weight[j]+2*((tau-1)/2-j);
}
sum+=weight[j];
printf("%lf\n",weight[j]);
}


for(int i=0;i<N-tau+1;i++)
{
 
for(int j=0;j<tau;j++)
{
 window[j] = b[j+i];
 wma[i] = wma[i]+window[j]*weight[j]/sum;

}

printf("wma(%d)%lf\n",i,wma[i]);
}

 


return 1;

}






























