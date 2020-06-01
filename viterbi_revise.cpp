#include <stdio.h>
#include <stdlib.h>
#include <math.h>

# define IA 16807
# define IM 2147483647
# define AM (1.0/IM)
# define IQ (127773)
# define IR 2836
# define NTAB 32
# define NDIV (1+(IM-1)/NTAB)
# define EPS 1.2e-7
# define RNMX (1.0-EPS)

double ran1(long *);
long *idum;

int main(void)
{
    
    int u[63]={0}; u[0]=1;                             /* u=generator information */
    unsigned int N;                                             /* N=total bits*/
    int i,j,k,m,p,q,t;                                 /*loop variable*/                      
    int inf;                                           /*inf=information*/  
    int s[6]={0};                                      /*s=state*/
    int temp_out1,temp_out2;                           /*temporary out*/
    int temps[6]={0};                 
    double *eno1,*eno2;                                /*enoput of encoder*/
    double x1,x2,cri,n1,n2,sigma;                      /*variables for normal function*/
    float SNR;
    long SEED;        
    int ps[64][32]={0};                                /*ps[previous state][time t-1]=information*/
    int cs[64][32];                                    /*cs[current state][time t]=information*/
    double metric[64]={0},pmetric[64]={0};             /*metric[state]*/  
    double *enc1,*enc2;
    int temp0_out1,temp0_out2,temp1_out1,temp1_out2;
    double out0,out1;                                  /*record the metric for each state*/ 
    int *finalout;                                           /*export the final output*/
    int csb[6]={0};                                    /*change the current state to binary*/  
    int rem;                                           /*reminder*/
    int small;                              /*for sorting*/
    double temp_small;
    int decision;
    double BER;
    int er=0;
    double num;
    
    printf("Please enter the randon seed(seed must be a NEGATIVE number): ");
    scanf("%d",&SEED);
    printf("Please enter the total number of decoded bits: ");
    scanf("%d",&N);
    printf("Please enter the SNR(in dB):");
    scanf("%f",&SNR);
    printf("Hard(1) or soft(2) decision: ");
    scanf("%d",&decision);

    idum=(long*)malloc(sizeof(long));
    *idum=SEED;                               
    eno1=(double*)malloc(sizeof(double)*N);
    eno2=(double*)malloc(sizeof(double)*N);  
    finalout=(int*)malloc(sizeof(int)*N);    
    sigma=sqrt(1/pow(10,(SNR/10)));

    for(j=0;j<63;j++)        /*generate sequence */
       u[j+6]=u[j+1]^u[j];                
                   
          
    for(i=0;i<N;i++)                          /*enceode information*/
    {
     do{                                       /*noise*/
        x1=ran1(idum);
        x2=ran1(idum);
        x1=2*x1-1;
        x2=2*x2-1;
        cri=x1*x1+x2*x2;
     }while (cri>=1.0);
     n1=sigma*x1*sqrt((-2)*log(cri)/cri);
     n2=sigma*x2*sqrt((-2)*log(cri)/cri);                                         
     
     inf=u[(i%63)];                                  /*information oput*/
     temp_out1=inf^s[1]^s[2]^s[4]^s[5];
     temp_out2=inf^s[0]^s[1]^s[2]^s[5];
     
/*map 0 to 1 ; 1 to -1 */             
     *(eno1+i)=(-2)*temp_out1+1+n1;
     *(eno2+i)=(-2)*temp_out2+1+n2;            
                                                                                                                                       
     temps[0]=s[0];temps[1]=s[1];temps[2]=s[2];temps[3]=s[3];temps[4]=s[4];temps[5]=s[5];   /*state transistion*/ 
     s[0]=inf;s[1]=temps[0];s[2]=temps[1];s[3]=temps[2];s[4]=temps[3];s[5]=temps[4];    
     
     if(decision==1)                        /*hard decision*/
       {
        if(*(eno1+i)>0)
           *(eno1+i)=1;
        else
           *(eno1+i)=-1;
        if(*(eno2+i)>0)
           *(eno2+i)=1;
        else
           *(eno2+i)=-1;
       }
    }    
    
/******************************************decoding***************************************/    
	
	
	
	
	for(t=0;t<6;t++)                              /*trallies diagram is not full*/
    {
       for(i=0;i<pow(2,t+1);i++)                  
       {
           csb[5]=i>>5;rem=i%32;
           csb[4]=rem>>4;rem=rem%16;
           csb[3]=rem>>3;rem=rem%8;
           csb[2]=rem>>2;rem=rem%4;
           csb[1]=rem>>1;
           csb[0]=i%2;
           temp0_out1=csb[0]^csb[2]^csb[3]^csb[5]^0;
           temp0_out2=csb[0]^csb[1]^csb[2]^csb[3]^0;
           temp0_out1=(-2)*temp0_out1+1;
           temp0_out2=(-2)*temp0_out2+1;
           out0=(temp0_out1-*(eno1+t))*(temp0_out1-*(eno1+t))+(temp0_out2-*(eno2+t))*(temp0_out2-*(eno2+t));
           out1=(temp1_out1-*(eno1+t))*(temp1_out1-*(eno1+t))+(temp1_out2-*(eno2+t))*(temp1_out2-*(eno2+t));
           cs[i][31]=csb[0];
           metric[i]=pmetric[(i>>1)]+out0;
           
           for(j=0;j<31;j++)
           {cs[i][j]=ps[(i>>1)][j+1];}
           
       }         
       for(m=0;m<64;m++)          
       {
           pmetric[m]=metric[m];
           for(p=0;p<32;p++)    
              ps[m][p]=cs[m][p];
       }
    }


   for(t=6;t<N;t++)                                   /*trallies is full path*/
   {   
       if(t>31)                                       /*start output*/
       {
          temp_small=metric[0];
          small=0;
          for(p=0;p<64;p++)                           /*find the smallest metric*/
          {
              if(metric[p]<temp_small)
              {
                 temp_small=metric[p];
                 small=p;
              }  
          }
           *(finalout+(t-32))=ps[small][0];
       }
       
       for(i=0;i<64;i++)
       {
           csb[5]=i>>5;rem=i%32;
           csb[4]=rem>>4;rem=rem%16;
           csb[3]=rem>>3;rem=rem%8;
           csb[2]=rem>>2;rem=rem%4;
           csb[1]=rem>>1;
           csb[0]=i%2;
           temp0_out1=csb[0]^csb[2]^csb[3]^csb[5]^0;
           temp0_out2=csb[0]^csb[1]^csb[2]^csb[3]^0;
           temp1_out1=csb[0]^csb[2]^csb[3]^csb[5]^1;
           temp1_out2=csb[0]^csb[1]^csb[2]^csb[3]^1;
           temp0_out1=(-2)*temp0_out1+1;
           temp0_out2=(-2)*temp0_out2+1;
           temp1_out1=(-2)*temp1_out1+1;
           temp1_out2=(-2)*temp1_out2+1;
           out0=(temp0_out1-*(eno1+t))*(temp0_out1-*(eno1+t))+(temp0_out2-*(eno2+t))*(temp0_out2-*(eno2+t));
           out1=(temp1_out1-*(eno1+t))*(temp1_out1-*(eno1+t))+(temp1_out2-*(eno2+t))*(temp1_out2-*(eno2+t));
           if((pmetric[(i>>1)]+out0)<=(pmetric[((i>>1)+32)]+out1))                   /*choose the upper path*/
             {
                 cs[i][31]=csb[0];                                            
                 metric[i]=pmetric[(i>>1)]+out0;
                 for(j=0;j<31;j++)
                     cs[i][j]=ps[(i>>1)][j+1];
             } 
            else
             {
                 cs[i][31]=csb[0];
                 metric[i]=pmetric[((i>>1)+32)]+out1;
                 for(j=0;j<31;j++)
                     cs[i][j]=ps[((i>>1)+32)][j+1];
            
             }
                      
       }
            
       for(m=0;m<64;m++)                                   
       {
           pmetric[m]=metric[m];
           for(p=0;p<32;p++)    
              ps[m][p]=cs[m][p];
       }
       
   
   }
   
/************************finish decoding*****************************************************/	
	
	
	for(i=0;i<(N-32);i++) 
	    er=er+(u[i%63]^*(finalout+i));
 
     num=(N-31);
     BER=er/num;
     
     printf("error=%d  BER=%e \n",er,BER);
	
	
	
	                                   
    free(eno1);
    free(eno2);
    free(finalout);

    system("pause");
    return 0;
    
}    



double ran1(long *idum)  
{
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  double temp;
 
  
  if(*idum <=0 || !iy){
     if (-(*idum) < 1) *idum=1;
     else *idum = -(*idum);
     for (j=NTAB+7;j>=0;j--){
         k=(*idum)/IQ;
         *idum=IA*(*idum-k*IQ)-IR*k;
         if (*idum<0) *idum += IM;
         if (j<NTAB)  iv[j] = *idum;
     }
     iy=iv[0];
  }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}  
