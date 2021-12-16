#include <stdio.h>
#include <math.h>

double q12,q22,gama,w,alpha,N,wL,w2,Omega;
double d,dd,h,t,Tr,fr,phi,VarreduraDoppler;
double Tp,a10,a20,area,x,y,LarguraDoppler,PassoDoppler;
int Pulsos,i,j,k,m,n,psi,PassoRKfento,g;
double a[4],b[4],c[4],k1[4],k2[4],k3[4],k4[4];
double Pi=3.141592653589793;

double f(double a1,double a2,double a12,double b12,int i)
{
    if (i==1) return            +2*Omega*(b12*cos(alpha)-a12*sin(alpha)) +q22*a2-(a1-1)*gama;
    if (i==2) return            -2*Omega*(b12*cos(alpha)-a12*sin(alpha)) -q22*a2    -a2*gama;
    if (i==3) return  -(d-dd)*b12-Omega*(a2-a1)*sin(alpha)               -q12*a12  -a12*gama;
    if (i==4) return   (d-dd)*a12+Omega*(a2-a1)*cos(alpha)               -q12*b12  -b12*gama;
    
}      

main() 
{
    FILE *arquivo;
    arquivo=fopen("dados.dat","w");
    
    //16-10-10
    //em função da frequência
    //Introdução da fase phi
    //Introdução da fase do caminho na cavidade wL*Tr
    //calculos numericos na presença do campo
    //calculos analiticos no decaimento
    
    gama=0.0025*1e9*0;
    PassoDoppler=0.2;
    VarreduraDoppler=3000;
    q22=(2*Pi)*5e6;
    q12=0.5*q22;
    fr=100e6;
    Tp=100e-15; 
    Omega=1*q22/(fr*Tp);
    PassoRKfento=100;
    LarguraDoppler=200;
    Pulsos=500;             
    a10=1;                  
    a20=0;                  
    w2= (2*Pi)*400e12;
    wL= (2*Pi)*400e12;
    phi=(2*Pi)*0.4*0;

    d=w2-wL;
    Tr=1/fr;                
    
    for (n=-(int)VarreduraDoppler;n<=(int)VarreduraDoppler;n++)
    {
        dd=(2*Pi)*PassoDoppler*1e6*n;
        N=-1;
        t=0;
        
        a[1]=a10; a[2]=a20; a[3]=0;  a[4]=0;
        b[1]=0;   b[2]=0;   b[3]=0;  b[4]=0;
    
    for (i=1;i<=2*Pulsos+1;i++)
       {
        if (i % 2 == 0)
             {
              N=N+1;    
              g=PassoRKfento;
              h=(1/double(g))*Tp;
                            
              alpha=-N*wL*Tr+N*phi;                       
              
                      for (k=1;k<=g;k++)
                       {
                        t=t+h;                
                                        
                        for (j=1;j<=4;j++){
                        k1[j]=f(a[1],a[2],a[3],a[4],j); }
                                                           
                        for (j=1;j<=4;j++){
                        k2[j]=f(a[1]+k1[1]*h/2,a[2]+k1[2]*h/2,a[3]+k1[3]*h/2,
                        a[4]+k1[4]*h/2,j);}
                  
                        for (j=1;j<=4;j++){
                        k3[j]=f(a[1]+k2[1]*h/2,a[2]+k2[2]*h/2,a[3]+k2[3]*h/2,
                        a[4]+k2[4]*h/2,j);}
                  
                        for (j=1;j<=4;j++){           
                        k4[j]=f(a[1]+k3[1]*h,a[2]+k3[2]*h,a[3]+k3[3]*h,
                        a[4]+k3[4]*h,j);}

                        for (j=1;j<=4;j++){
                        b[j]=a[j]+h*(k1[j]/6+k2[j]/3+k3[j]/3+k4[j]/6);}   
                  
                        for (m=1;m<=4;m++)
                        a[m]=b[m];     
                        } 
              
              }
        
        if (i % 2 == 1)
             {
              t=t+Tr-Tp;
              x=(d-dd)*(Tr-Tp);
              y=(Tr-Tp)*q22;
                                         
              b[1]=a[1]+a[2]*(1-exp(-y));
              b[2]=a[2]*exp(-y);
              b[3]=(a[3]*cos(x)-a[4]*sin(x))*exp(-0.5*y);
              b[4]=(a[3]*sin(x)+a[4]*cos(x))*exp(-0.5*y);
              
              for (m=1;m<=4;m++)
              a[m]=b[m];
              
              }     
        
        }           

        for (m=1;m<=4;m++)
          c[m]=b[m]*exp(-0.5*(dd)*(dd)/(2*Pi*LarguraDoppler*1e6)/
                                       (2*Pi*LarguraDoppler*1e6));                
 
        w=b[1]+b[2];  
        printf("\n");
        
        printf("%d %16.14f %16.14f %16.14f %16.14f", 
               n,c[1],c[2],sqrt(c[3]*c[3]+c[4]*c[4]),w);

        fprintf(arquivo,"%16.14f %16.14f %16.14f %16.14f %16.14f %16.14f",
               n*PassoDoppler,c[1],c[2],sqrt(c[3]*c[3]+c[4]*c[4]),c[3],c[4]);   
        fprintf(arquivo,"\n");
        
        for (m=1;m<=4;m++)
           c[m]=0;
        }
     
     fclose(arquivo);
     printf("\a");
}
