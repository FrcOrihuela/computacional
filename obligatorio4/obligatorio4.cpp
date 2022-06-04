#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <cstdlib>
#define PI 3.14159265

using namespace std;

double RungeKuttaIt(double (&nave)[4],double (&k)[4],double (&l)[4],double (&m)[4],double (&n)[4], double h, double w, double Delta, double mu, double rprima, double t);
double rpunto(double r, double phi, double pr, double pphi);
double phipunto(double r, double phi, double pr, double pphi);
double Prpunto(double r, double phi, double pr, double pphi, double w, double Delta, double mu, double rprima, double t);
double Pphipunto(double r, double phi, double pr, double pphi, double w, double Delta, double mu, double rprima, double t);


int main(void){

    double nave[4], k[4], l[4], m[4], n[4]; //Posicion y momentos de la nave y parámetros de Runge Kutta
    double Delta, mu, h, tol, t, dTL, MT, ML, G, w, RT, RL, modvel,phini,theta, r,rprima;
    bool hfija;

    int i, niteraciones;

    ofstream posicioncohete;

    posicioncohete.open("posicioncohete.txt");


    G=6.67*pow(10,-11);
    MT=5.9736*pow(10,24);
    ML=0.07349*pow(10,24);
    dTL=3.844*pow(10,8);
    w=2.6617*pow(10,-6);
    RT=6.378160*pow(10,6);
    RL=1.7374*pow(10,6);

    //Apolo11=46 768kg


    //Defino la h del paso
    h=1;

    //Voy a calcular los valores de delta y mu, necesarios para los siguientes calculos
    Delta=G*MT/(dTL*dTL*dTL);
    mu=ML/MT;


    //Valores iniciales

    modvel=1.119E4/dTL;
    phini=0.4548;
    theta=0.4548;
    r=RT/dTL;

    nave[0]=r;
    nave[1]=phini;
    nave[2]=modvel*cos(theta-phini);
    nave[3]=r*modvel*sin(theta-phini);
    t=0;
    posicioncohete<<"0, 0 \n";
    posicioncohete<<cos(w*t)<<", "<<sin(w*t)<<"\n";
    posicioncohete<<nave[0]<<", "<<nave[1]<<"\n";
    posicioncohete<<"\n";

    niteraciones=200000;

    for(i=0;i<niteraciones;i++){
        //Calculo rprima
        rprima=sqrt(1+nave[0]*nave[0]-2*nave[0]*cos(nave[1]-w*t));
        //Hago el calculo del metodo de Runge-Kutta
        t=t+RungeKuttaIt(nave,k,l,m,n,h,w,Delta,mu, rprima,t);
        
        //Escribo las posiciones de la tierra, la luna y del cohete, no pongo todas
        if(i%500==0){
        posicioncohete<<"0, 0 \n";
        posicioncohete<<cos(w*t)<<", "<<sin(w*t)<<"\n";
        posicioncohete<<nave[0]*cos(nave[1])<<", "<<nave[0]*sin(nave[1])<<"\n";
        posicioncohete<<"\n";
        }
    }



return 0;
}



double RungeKuttaIt(double (&nave)[4],double (&k)[4],double (&l)[4],double (&m)[4],double (&n)[4], double h, double w, double Delta, double mu, double rprima, double t){

    k[0]=h*rpunto(nave[0],nave[1], nave[2], nave[3]);
    l[0]=h*phipunto(nave[0],nave[1], nave[2], nave[3]);
    m[0]=h*Prpunto(nave[0],nave[1], nave[2], nave[3], w, Delta, mu, rprima, t);
    n[0]=h*Pphipunto(nave[0],nave[1], nave[2], nave[3], w, Delta, mu, rprima, t);

    k[1]=h*rpunto(nave[0]+k[0]/2,nave[1]+l[0]/2, nave[2]+m[0]/2, nave[3]+n[0]/2);
    l[1]=h*phipunto(nave[0]+k[0]/2,nave[1]+l[0]/2, nave[2]+m[0]/2, nave[3]+n[0]/2);
    m[1]=h*Prpunto(nave[0]+k[0]/2,nave[1]+l[0]/2, nave[2]+m[0]/2, nave[3]+n[0]/2, w, Delta, mu, rprima, t+h/2);
    n[1]=h*Pphipunto(nave[0]+k[0]/2,nave[1]+l[0]/2, nave[2]+m[0]/2, nave[3]+n[0]/2, w, Delta, mu, rprima, t+h/2);

    k[2]=h*rpunto(nave[0]+k[1]/2,nave[1]+l[1]/2, nave[2]+m[1]/2, nave[3]+n[1]/2);
    l[2]=h*phipunto(nave[0]+k[1]/2,nave[1]+l[1]/2, nave[2]+m[1]/2, nave[3]+n[1]/2);
    m[2]=h*Prpunto(nave[0]+k[1]/2,nave[1]+l[1]/2, nave[2]+m[1]/2, nave[3]+n[1]/2, w, Delta, mu, rprima, t+h/2);
    n[2]=h*Pphipunto(nave[0]+k[1]/2,nave[1]+l[1]/2, nave[2]+m[1]/2, nave[3]+n[1]/2, w, Delta, mu, rprima, t+h/2);

    k[3]=h*rpunto(nave[0]+k[2],nave[1]+l[2], nave[2]+m[2], nave[3]+n[2]);
    l[3]=h*phipunto(nave[0]+k[2],nave[1]+l[2], nave[2]+m[2], nave[3]+n[2]);
    m[3]=h*Prpunto(nave[0]+k[2],nave[1]+l[2], nave[2]+m[2], nave[3]+n[2], w, Delta, mu, rprima, t+h);
    n[3]=h*Pphipunto(nave[0]+k[2],nave[1]+l[2], nave[2]+m[2], nave[3]+n[2], w, Delta, mu, rprima, t+h);

    nave[0]=nave[0]+1./6*(k[0]+2*k[1]+2*k[2]+k[3]);
    nave[1]=nave[1]+1./6*(l[0]+2*l[1]+2*l[2]+l[3]);
    nave[2]=nave[2]+1./6*(m[0]+2*m[1]+2*m[2]+m[3]);
    nave[3]=nave[3]+1./6*(n[0]+2*n[1]+2*n[2]+n[3]);


    return h;
}


//Funciones de las ecuaciones diferenciales, las pongo aparte para poder calcular paso a paso
double rpunto(double r, double phi, double pr, double pphi){
    double rpunto;
    rpunto=pr;
    return rpunto;
}

double phipunto(double r, double phi, double pr, double pphi){
    double phipunto;
    phipunto=pphi/(r*r);
    return phipunto;
}

double Prpunto(double r, double phi, double pr, double pphi, double w, double Delta, double mu, double rprima, double t){
    double prpunto;
    prpunto=pphi*pphi/(r*r*r)-Delta*(1/(r*r)+mu/(rprima*rprima*rprima)*(r-cos(phi-w*t)));
    return prpunto;
}

double Pphipunto(double r, double phi, double pr, double pphi, double w, double Delta, double mu, double rprima, double t){
    double pphipunto;
    pphipunto=-Delta*mu*r*sin(phi-w*t)/(rprima*rprima*rprima);
    return pphipunto;
}

