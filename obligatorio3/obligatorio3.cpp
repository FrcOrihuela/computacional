#include <iostream>
#include <cmath>
#include <fstream>
#include <complex>
#define PI 3.141592654

using namespace std;

complex<double> calculaalfa(complex<double> alfa,double sgor, double Vj);
complex<double> calculabeta(complex<double> beta,complex<double> alfa, double sgor, double Vj,complex<double> phij);
complex<double> calculaxi(complex<double> alfa,complex<double> beta,complex<double> xi);
complex<double> phinij(double k0gor, int j, int N);
complex<double> calculaphi(complex<double> phij, complex<double> xi);
complex<double> calculaphi(complex<double> phij, complex<double> xi);





int main (void){

    int j,k,n,N,nciclos, iteraciones;
    //N nos da el ancho del pozo

    N=1000;

    //numero de ciclos

    nciclos=50;

    double k0,k0gor,sgor,V[N+1],lambda, vnorma, e,cnorm, amp;

    ofstream estados,norma;

    estados.open("estados.txt");
    norma.open("evolucionnorma.txt");


    //Utilizamos la libreria complex para usar numeros complejos tranquilamente
    complex<double> phi[N+1];
    complex<double> xi[N+1];
    complex<double> alfa[N],beta[N],gamma[N],Amenos,Acero[N],Amas,b[N+1];
    complex<double> i=(0.0,1.0);


    //Calculamos los valores de sgor k0gor

    
    k0gor=2*PI*nciclos/N;
    sgor=1/(4*k0gor*k0gor);
    lambda=1;

    //Definimos el potencial

    for(j=0;j<N+1;j++){
        if(j>(2*N/5)&&j<(3*N/5)){
            V[j]=lambda*k0gor*k0gor;
        }else{
            V[j]=0;    
        }
    }

    alfa[0]=0.+0.*i;
    alfa[N]=0.+0.*i;

    //Calculo las alfas
    Amas=1;
    Amenos=1;
    //Calculo las A0
    for ( j = 1; j< N; j++)
    {
        Acero[j]=complex<double>(-2.0-V[j], 2.0/sgor);
    }
    
    gamma[N-1]=1./Acero[N-1];
    for(j=N-2;j>=0;j--){
    alfa[j]=-Amenos*gamma[j+1];
    gamma[j]=1./(alfa[j]+Acero[j]);
    }



    //Numero de iteraciones
    iteraciones=1000;

    //Valor de beta final
    beta[N]=0;


    //Calculamos el estado inicial
    phi[0]=0;
    phi[N]=0;
    cnorm=0;

    for(j=1;j<N;j++){
    e=exp(-8.0*(4*j-N)*(4*j-N)/(1.0*N*N));
    phi[j]=complex<double>(cos(k0gor*j),sin(k0gor*j))*e;
    cnorm=cnorm+real(phi[j])*real(phi[j])+imag(phi[j])*imag(phi[j]);
    }
    cnorm=sqrt(cnorm);

    //Normalizo
    for(j=1;j<N;j++){
        phi[j]=phi[j]/cnorm;
    }


    vnorma=0;

    for(k=0;k<=N;k++){
        amp=real(phi[k])*real(phi[k])+imag(phi[k])*imag(phi[k]);
        estados<<k<<", "<<amp<<"\n";
        vnorma=vnorma+amp;
    }
    estados<<"\n";
    norma<<0<<"\t"<<sqrt(vnorma)<<"\n";
    estados.flush();
    norma.flush();

    for(n=1;n<=iteraciones;n++){
        
        //Primero calculo las betas
        for(j=1;j<N;j++){
            b[j]=complex<double>(0,4)*phi[j]/sgor;
        }

        for(j=N-2;j>=0;j--){
            beta[j]=gamma[j+1]*(b[j+1]-beta[j+1]);
        }

        //Ahora calculo las xi
        xi[0]=0;
        xi[N]=0;
        for(j=1;j<N;j++){
            xi[j]=alfa[j-1]*xi[j-1]+beta[j-1];
        }

        //Y calculo las phi
        cnorm=0;
        for(k=1;k<N;k++){
            phi[k]=xi[k]-phi[k];
         }



        //Meto la iteracion en el txt de resultados y calcula la norma
        vnorma=0;
        for(k=0;k<=N;k++){
            amp=real(phi[k])*real(phi[k])+imag(phi[k])*imag(phi[k]);
            estados<<k<<", "<<amp<<"\n";
            vnorma=vnorma+amp;
        }
        estados<<"\n";
        norma<<n<<"\t"<<sqrt(vnorma)<<"\n";
        estados.flush();
        norma.flush();
    }


    return 0;
}

//En esta función hago el cálculo del siguiente valor de xi a partir de los anteriores

complex<double> calculaalfa(complex<double> alfa,double sgor, double Vj){

    complex<double> alfamenos,Amenos,Acero,Amas, gamma;
    complex<double> i=(0.0,1.0);
    //Defino las Acosas

    Amenos=1;
    Acero=complex<float>(-2.0-Vj, 2.0/sgor);
    Amas=1;

    gamma=Acero+Amas*alfa;
    gamma=1./gamma;

    alfamenos=-Amenos*gamma;

    return alfamenos;
}

complex<double> calculabeta(complex<double> beta,complex<double> alfa,double sgor, double Vj,complex<double> phij){
    complex<double> betamenos,Amenos,Acero,Amas, gamma, bejota;
    complex<double> i=(0.0,1.0);
    //Defino las Acosas

    Amenos=1;
    Acero=complex<float>(-2.0-Vj, 2.0/sgor);
    Amas=1;

    gamma=Acero+Amas*alfa;
    gamma=1./gamma;
    bejota=complex<double>(-4*imag(phij)/sgor,4*real(phij)/sgor);

    betamenos=gamma*(bejota-beta);

    return betamenos;
}

complex<double> calculaxi(complex<double> alfa,complex<double> beta,complex<double> xi){
    complex<double> ximas;

    ximas=alfa*xi+beta;

    return ximas;
}


complex<double> calculaphi(complex<double> phij, complex<double> xi){
    complex<double> phijmas;

    phijmas=xi-phij;


    return phijmas;
}

