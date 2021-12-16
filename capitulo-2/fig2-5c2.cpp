#include <stdio.h>
#include <math.h>

double q12,q22,gama,w,alpha,N,wL,w2;
double dd,Omega,h,t,Tr,fr,phi;
double Tp,a10,a20,x,y;
int Pulsos,i,j,k,m,n,psi,PassoRKfento,g,s;
double a[4],b[4],c[4],k1[4],k2[4],k3[4],k4[4];
double Pi=3.141592654;     

double f(double t,double a1,double a2,double a12,double b12,int j)
{
    if (j==1) return            +2*Omega*(b12*cos(alpha)-a12*sin(alpha))  +q22*a2 -(a1-a10)*gama;
    if (j==2) return            -2*Omega*(b12*cos(alpha)-a12*sin(alpha))  -q22*a2       -a2*gama;
    if (j==3) return -dd*b12-Omega*(a2-a1)*sin(alpha)                    -0.5*q22*a12  -a12*gama;
    if (j==4) return  dd*a12+Omega*(a2-a1)*cos(alpha)                    -0.5*q22*b12  -b12*gama;
}      

main() 
{
    FILE *arquivo;
    arquivo=fopen("dados.dat","w");
    
    //14-10-10
    //Corrigido em 21-02-11 (erro em rho22)
    //Introdução da fase phi
    //Introdução da fase do caminho na cavidade wL*Tr
    //calculos numericos na presença do campo
    //calculos analiticos no decaimento
    //nome da pasta: 2N_tempo
    //Pi = 3.141592654
    
    gama=0.0025*1e9*0;
    PassoRKfento=100;
    q22=(2*Pi)*5e6;
    q12=0.5*q22;
    fr=100e6;
    Tp=100e-15;             
    Pulsos=100;             
    Omega=0.01*q22/(fr*Tp);
    a10=1;                  
    a20=0;                  
    w2=(2*Pi)*400e12;
    wL=(2*Pi)*400e12;
    phi=(2*Pi)*0.4*0;
              
    a[1]=a10;
    a[2]=a20;
    a[3]=0;
    a[4]=0;    
    t=-50e-15;N=-1;
    dd=w2-wL;
    Tr=1/fr;                //Intervalo entre os pulsos}
    
 
    for (i=0;i<=2*Pulsos+1;i++)
       {
        //printf("*");
        if (i % 2 == 0)
             {
              N=N+1;    
              g=PassoRKfento;
              h=1e-15;
                            
              alpha=-N*wL*Tr+N*phi;                       
              
                      for (k=1;k<=g;k++)
                       {
                        t=t+h;                
                                        
                        for (j=1;j<=4;j++){
                        k1[j]=f(t,a[1],a[2],a[3],a[4],j); }
                                                           
                        for (j=1;j<=4;j++){
                        k2[j]=f(t+h/2,a[1]+k1[1]*h/2,a[2]+k1[2]*h/2,a[3]+k1[3]*h/2,
                        a[4]+k1[4]*h/2,j);}
                  
                        for (j=1;j<=4;j++){
                        k3[j]=f(t+h/2,a[1]+k2[1]*h/2,a[2]+k2[2]*h/2,a[3]+k2[3]*h/2,
                        a[4]+k2[4]*h/2,j);}
                  
                        for (j=1;j<=4;j++){           
                        k4[j]=f(t+h,a[1]+k3[1]*h,a[2]+k3[2]*h,a[3]+k3[3]*h,
                        a[4]+k3[4]*h,j);}

                        for (j=1;j<=4;j++){
                        b[j]=a[j]+h*(k1[j]/6+k2[j]/3+k3[j]/3+k4[j]/6);}   
                        
                        //fprintf(arquivo,"%20.18f %20.18f \n", t, sech(1.763*t/Tp));
                        
                        for (m=1;m<=4;m++)
                        a[m]=b[m];    
                        
                        w=b[1]+b[2];  
                
                        printf("%d %12.10f %12.10f %12.10f %12.10f", 
                               i,b[1],b[2],sqrt(b[3]*b[3]+b[4]*b[4]),w);
                        printf("\n");       
                        fprintf(arquivo,"%12.10f %12.10f %12.10f %12.10f",
                            t*1e15 - 9999900*0,b[1],b[2],sqrt(b[3]*b[3]+b[4]*b[4]));   
                        fprintf(arquivo,"\n"); 
                         
                        } 
              
              }
        
        if (i % 2 == 1)
             {
              t=t+(Tr-Tp);
              x=dd*(Tr-Tp);
              y=(Tr-Tp)*q22;
                                         
              b[1]=a[1]+a[2]*(1-exp(-y));
              b[2]=a[2]*exp(-y);
              b[3]=(a[3]*cos(x)-a[4]*sin(x))*exp(-0.5*y);
              b[4]=(a[3]*sin(x)+a[4]*cos(x))*exp(-0.5*y);
              
              for (m=1;m<=4;m++)
              a[m]=b[m];
              
              }                   
     
        w=b[1]+b[2];  
  
        printf("%d %12.10f %12.10f %12.10f %12.10f", 
                i,b[1],b[2],sqrt(b[3]*b[3]+b[4]*b[4]),w);
        printf("\n");       
        fprintf(arquivo,"%12.10f %12.10f %12.10f %12.10f",
                t*1e15 - 9999900,b[1],b[2],sqrt(b[3]*b[3]+b[4]*b[4]));   
        fprintf(arquivo,"\n");   
              
        }
        //scanf("%d",&s);
        fclose(arquivo);              
        printf("\a");

}
