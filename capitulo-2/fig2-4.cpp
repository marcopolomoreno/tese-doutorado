#include <stdio.h>
#include <math.h>

double q12,q22,w,alpha,N,wL,w2;
double dd,Omega,h,t,Tr,fr,phi;
double Tp,x,y;
int Pulsos,i,j,k,m,PassoRKfento,g;
double a[4],b[4],k1[4],k2[4],k3[4],k4[4];
double Pi=3.141592654;

double f(double a1,double a2,double a12,double b12,int i)
{
    if (i==1) return            +2*Omega*(b12*cos(alpha)-a12*sin(alpha))  +q22*a2;
    if (i==2) return            -2*Omega*(b12*cos(alpha)-a12*sin(alpha))  -q22*a2;
    if (i==3) return -dd*b12-Omega*(a2-a1)*sin(alpha)                    -q12*a12;
    if (i==4) return  dd*a12+Omega*(a2-a1)*cos(alpha)                    -q12*b12;
    
}

main() 
{
    FILE *arquivo;
    arquivo=fopen("dados.dat","w");
    
    //14-10-10
    //Corrigido em 21-02-11 (erro em rho22)
    //Introdu??o da fase phi
    //Introdu??o da fase do caminho na cavidade wL*Tr
    //calculos numericos na presen?a do campo
    //calculos analiticos no decaimento
    //nome da pasta: 2N_tempo
    //Pi = 3.141592654
    
    PassoRKfento=10;
    q22=(2*Pi)*5e6;
    q12=0.5*q22;
    fr=100e6;
    Tp=100e-15;             
    Pulsos=50;             
    Omega=0.01*q22/(fr*Tp);                     
                
    w2=(2*Pi)*400e12;
    wL=(2*Pi)*400e12;
    phi=(2*Pi)*0.4*0;
              
    a[1]=1;
    a[2]=0;
    a[3]=0;
    a[4]=0;    
    t=0;N=-1;
    dd=w2-wL;
    Tr=1/fr;                //Intervalo entre os pulsos}
    
    //fprintf(arquivo,"%12.10f %12.10f %12.10f %12.10f",
    //            t,a[1],a[2],sqrt(a[3]*a[3]+a[4]*a[4]));   
    //fprintf(arquivo,"\n");
    
    for (i=1;i<=2*Pulsos+1;i++)
       {
        //printf("*");
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
                t*1e9-0.99999,b[1],b[2],sqrt(b[3]*b[3]+b[4]*b[4]));   
        fprintf(arquivo,"\n");   
              
        }
        
        fclose(arquivo);              
        printf("\a");
        

}
