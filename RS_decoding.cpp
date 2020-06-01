#include <stdio.h>
#include <stdlib.h>


int logtable[31]={1,2,4,8,16,5,10,20,13,26,17,7,14,28,29,31,
                  27,19,3,6,12,24,21,15,30,25,23,11,22,9,18};
/*int gx[17]={14,3,21,15,29,15,16,20,25,24,2,8,13,1,28,15,1}; 
*/
int gflog(int);
int invlog(int);
int gfmul(int,int);
int gfdiv(int,int);
int polydeg(int*);
int *polydiv(int*,int*,int*);
int *polyadd(int*,int*);
int *polymul(int*,int*); 
void Euclid(int*,int*,int,int,int*,int*);


int main()
{
    /*input data*/
	char input[31]={14 ,20 ,13 ,25 ,4 ,8 ,26 ,11 ,31 ,6 ,31 ,18 ,1 ,5 ,17 ,11 ,5 ,27 ,15 ,0 ,9 ,0 ,1 ,11 ,17, 0, 0 ,28 ,27 ,18, 22};
	int decode=0,count=0;  /*decode==1 =>true  decode==-1=>false*/      
	int i,j,u;/*loop variable*/
	int Rprime[31]={0};
	int *sigma=(int*)malloc(sizeof(int)*32);   /*sigma=sigma0(x)erasure locator polynomial*/
	int sigma1[32]={0};/*sigma1(x)*/
	int *sigmax=(int*)malloc(sizeof(int)*32); /*sigma(x)=sigma0(x)*sigma1(x)*/
	int temp[32]={0};
	int temp2[32]={0};temp2[0]=1;
	int syndrome[32]={0};/*syndrome(x)=s(x)*/
	int eranum=0;/*erasure number*/   
    int *sx=(int*)malloc(sizeof(int)*32); /*sx=s0(x)=sigma0(x)*syndrome(x)*/
	int mu,nu;
	int xr[32]={0};xr[16]=1; /*xr=x^r*/
	int deg_sx;
	int omega[32]={0};/*oemga(x)*/
	int sum_sig=0,sum_sigp=0,sum_ome=0;
	int sigmaprime[32]={0};
	int Er[32]={0};/*Error*/
	int output[31]={0};
	int tempq[32]={0};
	int degsyn=0;
	
	/*compute the erasure locator polynomial*/   	    
	for(i=0;i<31;i++)
	   {
		if(input[i]=='*')
		  {
		   eranum=eranum+1;	
		   temp[0]=1;
		   temp[1]=logtable[i];
		   sigma=polymul(temp2,temp);
		   Rprime[i]=0;
		   for(j=0;j<32;j++)
		      {
			   temp2[j]=sigma[j];	
		      }
		  }
		else
		   Rprime[i]=input[i];
	   }
	if(eranum==0)  /*no errasure sigma0=1*/
      {
       sigma[0]=1;          
       for(i=1;i<32;i++)          
       sigma[i]=0;
       }             
	/*compute syndrome s[0]=s1 s[1]=s2 etc...*/
    if(eranum<17)/*  check whwther the erasure is out of the bound*/
	  {
	   for(j=0;j<16;j++)
	      {
		   for(i=0;i<31;i++)
		       syndrome[j]=syndrome[j]^gfmul(Rprime[i],logtable[((i*(j+1))%31)]);	 		   
	      }
	   for(i=15;i>=0;i--)
	      {
		   if(syndrome[i]!=0)
		     {
			  degsyn=i;
			  break;		
			 }   
		  }
	   if(degsyn==0 && syndrome[0]==0)  
	     {
		  		
		  decode=1;
		  for(j=0;j<31;j++)
		      Er[j]=0;	
	     }
	   else
	     {	 
	      sx=polydiv(polymul(syndrome,sigma),xr,tempq);
          mu=(16-eranum)/2;
	      nu=15-mu;    
          Euclid(xr,sx,mu,nu,sigma1,omega);
    
          u=sigma1[0]; 
          for(i=0;i<32;i++)
             {
		      sigma1[i]=gfdiv(sigma1[i],u);
		      omega[i]=gfdiv(omega[i],u);
	         }	

          sigmax=polymul(sigma,sigma1);   
    
	      for(i=0;i<31;i=i+2)
		      sigmaprime[i]=sigmax[i+1];		
        
          /*time-domain completion*/
          if(sigmax[0]==0 || polydeg(omega)>=eranum+polydeg(sigma1))
	         decode=-1;
       
          else
	        {
	         count=0;
	         for(i=0;i<31;i++)
	            {
		         sum_sig=0;
		         sum_sigp=0;
		         sum_ome=0;	
		         for(j=0;j<32;j++)
		            {			   	
			         sum_sig=sum_sig^gfmul(sigmax[j],logtable[(j*(31-i))%31]);		
		             sum_sigp=sum_sigp^gfmul(sigmaprime[j],logtable[(j*(31-i))%31]);
		             sum_ome=sum_ome^gfmul(omega[j],logtable[(j*(31-i))%31]);
			        }
			    
		         if(sum_sig==0 && sum_sigp!=0)
		           {
			        count=count+1;
			        Er[i]=gfdiv(sum_ome,sum_sigp);
			       }
		         else
		            Er[i]=0; 		
		        }
	         if(count==polydeg(sigmax))
	            decode=1;
	         else
	            decode=-1;  	     
	        }    
    
          }
      }
	else/*out of boubd*/
	  decode=-1;

    /*check is finish*/	    
    if(decode==1)
      {
	   printf("TRUE\n");
       for(i=0;i<31;i++)
          {
		   output[i]=Rprime[i]^Er[i];
           printf("%d ",output[i]);
		  }
	  }
    else
	  printf("FALSE");
	    
	   
    system("pause"); 
    return 0;
}
/*comoute the degree of the polynomail*/
int polydeg(int *poly)
{
	int i;
	int degree;
	
	for(i=31;i>=0;i--)
	   {
		if(poly[i]!=0)
		  {
		   degree=i;
		   break;
		  }		
	   }
	return degree;
}
/*perform the Euclid's algorithm*/	
void Euclid(int *ax,int *bx,int mu,int nu,int *sigmax,int *omega)
{
	int *rix=(int*)malloc(sizeof(int)*32);
	int *r1x=(int*)malloc(sizeof(int)*32);
	int *r2x=(int*)malloc(sizeof(int)*32);
	int *vix=(int*)malloc(sizeof(int)*32);
	int *v1x=(int*)malloc(sizeof(int)*32);
	int *v2x=(int*)malloc(sizeof(int)*32);
	int qix[32]={0};
	int i,j,p,q;/*loop variable*/
	int degr=30,degv=0;
	
	for(i=0;i<32;i++)
	   {
		rix[i]=0;
	    r2x[i]=ax[i];		
	    r1x[i]=bx[i];
	    vix[i]=0;
	    v2x[i]=0;		
	    v1x[i]=0;
	   }
	v1x[0]=1;
	
	degv=polydeg(v1x);
	degr=polydeg(r1x);
	/*Chenck mu and nu first*/   	
	if(degr<=nu && degv<=mu)
	  {		
	   for(i=0;i<32;i++)
	      {	
		   sigmax[i]=v1x[i];
		   omega[i]=r1x[i];
		  }
	  }
	else
	  {			
		for(j=0;j<32;j++)
		   {
			rix=polydiv(r2x,r1x,qix);/*r_i=r_(i-2)-r_(i-1)*qi*/
			vix=polyadd(v2x,polymul(qix,v1x));/*v_i=v_(i-2)-v_(i-1)*qi*/
			
            degr=polydeg(rix);
            degv=polydeg(vix);

			if(degr<=nu && degv<=mu)
	          {
	           for(i=0;i<32;i++)
	              {
		           sigmax[i]=vix[i];
		           omega[i]=rix[i];
		          }
			   break;
			  }
			else
			  {
			   for(i=0;i<32;i++)
	              {
		           r2x[i]=r1x[i];
			       r1x[i]=rix[i];
			       v2x[i]=v1x[i];
			       v1x[i]=vix[i];
			       qix[i]=0;
		          }		
			  }			     
		   }
	  }		   	   	  		  		
      free(rix);
      free(r1x);
	  free(r2x);
	  free(vix);
	  free(r1x);
	  free(r2x);
}	

/*polynomail division*/
int *polydiv(int *pdividend,int *pdivisor,int *pquotient)
{
	int *temp1=(int*)malloc(sizeof(int)*32);
	int *remainder=(int*)malloc(sizeof(int)*32);
	int i,j;/*loop variabele*/
	int degpdividend,degpdivisor,degremainder;
	
	for(i=0;i<32;i++)
	   {
		temp1[i]=0;
	    remainder[i]=0;
	   }
	
	remainder=pdividend;
	
    degpdividend=polydeg(pdividend);
    degpdivisor=polydeg(pdivisor);
	degremainder=degpdividend;   

    for(j=(degpdividend-degpdivisor);j>=0;j--)
	   {
		pquotient[j]=gfdiv(remainder[degremainder],pdivisor[degpdivisor]);
		temp1[j]=pquotient[j];
		remainder=polyadd(polymul(pdivisor,temp1),remainder);
		degremainder=degremainder-1;
		/*compare the degree of remainder and divisor*/
		if(degremainder<degpdivisor)
		   break;
		else
		   temp1[j]=0;
	   }
	free(temp1);
	return remainder;
}				
		 			
/*find corresponding log value of a number*/
int gflog(int inval)
{
    int outpow; /*for input value and output log*/
    int p=0; /*loop variable*/
        
    for(p=0;p<31;p++)
    {  
        if(inval==logtable[p])
        {
           outpow=p;
           p=32;
        }
    }                    
            
    return outpow;
}

/*inverse the log computation*/
int invlog(int inpow)
{
    int outval;
     
    outval=logtable[(inpow%31)];
    return outval;
}

/*perform the addition in GF(32)*/
int gfmul(int mult1,int mult2)
{
    int result;
    
    if(mult1==0 || mult2==0)
       result=0;
    else
       result=invlog(gflog(mult1)+gflog(mult2));
 
    return result; 
}
/*division in GF(32)*/
int gfdiv(int dividend,int divisor)
{
    int quotient;
    int temp;
    
    if(dividend==0)
       quotient=0;
    else if(divisor==0)
       quotient=-1;
    else
       {
        temp=gflog(dividend)-gflog(divisor);
        if(temp<0)
          {
           temp=temp+31;
           quotient=invlog(temp);
          }
          else
          quotient=invlog(temp);
       }
    return quotient;  
}
/*add the two polynomials*/
int *polyadd(int *addpoly1,int *addpoly2)
{
	int *sumpoly=(int*)malloc(sizeof(int)*32);
	int i;
	
	for(i=0;i<32;i++)
		sumpoly[i]=addpoly1[i]^addpoly2[i];
		
	return sumpoly;
}	
/*multiply two polynomails*/
int *polymul(int *mulpoly1,int *mulpoly2)
{
    int *outpoly=(int*)malloc(32*sizeof(int));
    int temp; /*the coefficient of the specific x power*/
    int i,j;/*loop variable*/
         
    
    for(i=0;i<32;i++)
       {
        temp=0;
        for(j=0;j<i+1;j++)
            temp=temp^gfmul(mulpoly1[j],mulpoly2[i-j]);
                             
        outpoly[i]=temp;
       }                
    return outpoly;
}
