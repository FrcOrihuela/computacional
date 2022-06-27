#include <iostream>
#include <cmath>
#include <fstream>
#include <complex>
#include <cstdlib>
#include <ctime>
#define PI 3.141592654

using namespace std;

complex<double> calculaalfa(complex<double> alfa,double sgor, double Vj);
complex<double> calculabeta(complex<double> beta,complex<double> alfa, double sgor, double Vj,complex<double> phij);
complex<double> calculaxi(complex<double> alfa,complex<double> beta,complex<double> xi);
complex<double> phinij(double k0gor, int j, int N);
complex<double> calculaphi(complex<double> phij, complex<double> xi);
complex<double> calculaphi(complex<double> phij, complex<double> xi);





int main (void){

    int j,k,n,N,nciclos, iteraciones, h, mT, nD, contador, salir;
    //N nos da el ancho del pozo

    N=500;

    //numero de ciclos

    nciclos=50;

    double k0,k0gor,sgor,V[N+1],lambda, vnorma, e,cnorm, amp, paso;
    double PR,PL,generador, K,probmedia, desvtipica;


    ofstream estados,norma, potencial,valorK, posicionmedia, momentomedio, cineticamedia, energiamedia;

    estados.open("estados.txt");
    norma.open("revolucionnorma.txt");
    valorK.open("ValorK.txt");
    posicionmedia.open("posicionmedia.txt");
    momentomedio.open("momentomedio.txt");
    cineticamedia.open("cineticamedia.txt");
    energiamedia.open("energiamedia.txt");
    potencial.open("potencial.txt");
    
    //Utilizamos la libreria complex para usar numeros complejos tranquilamente
    complex<double> phi[N+1];
    complex<double> xi[N+1];
    complex<double> alfa[N],beta[N],gamma[N],Amenos,Acero[N],Amas,b[N+1];
    complex<double> dphi[N+1],ddphi[N+1]; //Primera y segunda derivadas, vectores complejos
    //Valores esperados
    complex<double> posicionesp, posicionerr, momentoesp, momentoerr, cineticaesp, cineticaerr, energiaesp, energiaerr;
    complex<double> i=(0.0,1.0);

    //Semilla valores aleatorios
    srand(time(NULL));


    //Calculamos los valores de sgor k0gor

    
    k0gor=2*PI*nciclos/N;
    sgor=1/(4*k0gor*k0gor);
    //calculo el paso entre puntos
    paso=2*k0gor;
    lambda=1;
    nD=500; //TIEMPO PARA QUE PUEDA APARECER UN MÁXIMO AL OTRO LADO DE LA BARRERA DE POTENCIAL


    //Definimos el potencial
    
    for(j=0;j<N+1;j++){
        if(j>(2*N/5)&&j<(3*N/5)){
            V[j]=lambda*k0gor*k0gor;
            potencial<<j<<"\t"<<V[j]<<"\n";
        }else{
            V[j]=0; 
            potencial<<j<<"\t"<<V[j]<<"\n";   
        }
    }

    potencial.flush();
    potencial.close();


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
    mT=0; //Inicializo el valor
    probmedia=0;
    desvtipica=0;
    //Valor de beta final
    beta[N]=0;


for(h=0;h<iteraciones;h++){


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
/*
    for(k=0;k<=N;k++){
        amp=real(phi[k])*real(phi[k])+imag(phi[k])*imag(phi[k]);
        estados<<k<<", "<<amp<<"\n";
        vnorma=vnorma+amp;
    }
    estados<<"\n";
    norma<<0<<"\t"<<sqrt(vnorma)<<"\n";
    estados.flush();
    norma.flush();
*/
    contador=1;
    salir=0;
    while(salir==0){
        contador++;
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

        for(k=1;k<N;k++){
            phi[k]=xi[k]-phi[k];
         }



        //Meto la iteracion en el txt de resultados y calcula la norma
        /*for(k=0;k<=N;k++){
            amp=real(phi[k])*real(phi[k])+imag(phi[k])*imag(phi[k]);
            estados<<k<<", "<<amp<<"\n";
            vnorma=vnorma+amp;
        }
        estados<<"\n";
        norma<<n<<"\t"<<sqrt(vnorma)<<"\n";
        estados.flush();
        norma.flush();*/







    //Una vez evolucionada la partícula tras t=nD medimos
    if((contador%nD)==0){
        
        //Calculo PR
        PR=0;
        for(j=0.8*N;j<N;j++){
            PR=PR+real(phi[j])*real(phi[j])+imag(phi[j])*imag(phi[j]);
            }
        probmedia=probmedia+PR;
        desvtipica=desvtipica+PR*PR;

        generador=1.0*(rand()%RAND_MAX)/RAND_MAX;
        //cout<<PR<<"\t"<<generador<<endl;
        if(generador<PR){
            mT++;
            //cout<<mT<<endl;
            salir=1;
        }else{
            for(j=0.8*N;j<N;j++){
                phi[j]=0;
            }
            vnorma=0;
            for(j=1;j<N;j++){
                vnorma=vnorma+real(phi[j])*real(phi[j])+imag(phi[j])*imag(phi[j]);
            }
            vnorma=sqrt(vnorma);
            for(j=1;j<N;j++){
                phi[j]=phi[j]/vnorma;
            }
            //CALCULO PL
            PL=0;
            for(j=1;j<=0.2*N;j++){
                PL=PL+real(phi[j])*real(phi[j])+imag(phi[j])*imag(phi[j]);
            }
            //genero otro numero aleatorio para medir a la izquierda
            generador=1.0*(rand()%RAND_MAX)/RAND_MAX;
            if(generador<PL){
                salir=1; //Que acabe el bucle pq se ha detectado la particula
            }else{
                for(j=1;j<=N/5;j++){
                phi[j]=0;
                }
            vnorma=0;
            for(j=1;j<N;j++){
                vnorma=vnorma+real(phi[j])*real(phi[j])+imag(phi[j])*imag(phi[j]);
            }
            vnorma=sqrt(vnorma);
            for(j=1;j<N;j++){
                phi[j]=phi[j]/vnorma;
            }
            }
        }

    
    }

    //Calculo de valores esperados de operadores:

    
    //Primero calculo las derivadas utilizando fórmulas de derivación numerica
    
    dphi[0]=phi[1]/e;
    dphi[N]=(phi[N]-phi[N-1])/paso;
    ddphi[0]=0;
    ddphi[N]=0;
    for ( j = 1; j < N; j++){
        dphi[j]=(phi[j+1]-phi[j-1])/(2*paso);
        ddphi[j]=(phi[j+1]-2.0*phi[j]+phi[j-1])/(paso*paso);
    }

    posicionesp=0;
    posicionerr=0;
    momentoesp=0;
    momentoerr=0;
    cineticaesp=0;
    cineticaerr=0;
    energiaesp=0;
    energiaerr=0;

    for(j=0;j<=N;j++){
        posicionesp=posicionesp+j*paso*(conj(phi[j])*phi[j]);
        momentoesp=momentoesp+(conj(phi[j])*dphi[j]);
        cineticaesp=cineticaesp+(conj(phi[j])*ddphi[j]);
        energiaesp=energiaesp+(-conj(phi[j])*ddphi[j])+complex<double>(V[j],0)*(conj(phi[j])*phi[j])/(paso*paso);

        posicionerr=posicionerr+(j*j*paso*paso*(conj(phi[j])*phi[j]));
        momentoerr=momentoerr+(conj(phi[j])*dphi[j])*(conj(phi[j])*dphi[j]);
        cineticaerr=cineticaerr+(conj(phi[j])*ddphi[j])*(conj(phi[j])*ddphi[j]);
        energiaerr=energiaerr+(-(conj(phi[j])*ddphi[j])+complex<double>(V[j],0)*(conj(phi[j])*phi[j])/(paso*paso))*(-(conj(phi[j])*ddphi[j])+complex<double>(V[j],0)*(conj(phi[j])*phi[j])/(paso*paso));
    }

    posicionerr=sqrt(posicionerr-posicionesp*posicionesp);
    momentoerr=sqrt(momentoerr-momentoesp*momentoesp);
    cineticaerr=sqrt(cineticaerr-cineticaesp*cineticaesp);
    energiaerr=sqrt(energiaerr-energiaesp*energiaesp);


    posicionmedia<<contador<<"\t"<<real(posicionesp)<<"\t"<<real(posicionerr)<<"\n";
    momentomedio<<contador<<"\t"<<imag(momentoesp)<<"\t"<<real(momentoerr)<<"\n";
    cineticamedia<<contador<<"\t"<<norm(cineticaesp)<<"\t"<<norm(cineticaerr)<<"\n";
    energiamedia<<contador<<"\t"<<norm(energiaesp)<<"\t"<<norm(energiaerr)<<"\n";

    posicionmedia.flush();
    momentomedio.flush();
    cineticamedia.flush();
    energiamedia.flush();
    


}

}

K=1.0*mT/iteraciones;
probmedia=probmedia/iteraciones;
desvtipica=sqrt((1-K)*(1-K)*K+K*K*(1-K));//Esta es la desv típica en una distribucion de Bernoulli
valorK<<"Valor de K:\t"<<K<<endl;
valorK<<"Probabilidad promedio:\t"<<probmedia<<endl;
valorK<<"Desviación Típica:\t"<<desvtipica<<endl;
desvtipica=desvtipica*2.0/sqrt(iteraciones);
valorK<<"Error:\t"<<desvtipica<<endl;
valorK.flush();
valorK.close();


posicionmedia.close();
momentomedio.close();
cineticamedia.close();
energiamedia.close();


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

