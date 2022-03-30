//Programa obligatorio 1
//Simulacion del sistema solar utilizando el algoritmo de velvet

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

using namespace std;


float verlet2d(float (&p)[10][2],float (&v)[10][2], float (&a)[10][2],float h,int Np);
float newton(float (&p)[10][2],float (&v)[10][2], float (&a)[10][2], int Np);
float cambiomasa(float m);
float masanormal(float m);
float cambiodistancia(float d);
float distancianormal(float d);
float tiemponormal(float t,float c,float G,float Ms);


int main(void){

    float h, t, G, c, Ms, cambiosigno, energiaini[10], energia[10], masas[10], modvel, moddist, energiatot;
    float p[10][2];
    float v[10][2];
    float a[10][2];
    float yprev[10]; //Vector para guardar los valores del eje y para calcular los periodos
    float periodo[10]; //Para guardar los valores de los periodos
    int nvueltas[10]; //guarda el numero de vueltas que da cada planeta
    int n, i, j, Np; //contadores y el numero de planetas
    
    ofstream posplat;
    ofstream energias;
    ifstream condini;
    ofstream geocent;
    ofstream periodos;
    //Defino la constante gravitatoria universal, la masa del sol
    G=6.67*pow(10,-11);
    c=1.496*pow(10,11);
    Ms=2*pow(10,30);


    //Este fichero va a guardar las posiciones y las velocidades de cada planeta, que se iran escribiendo en el ciclo for

    posplat.open("posicionesplanetas.txt");
    condini.open("condicionesiniciales.txt");
    energias.open("energiasplanetas.txt");
    geocent.open("geocentrico.txt");
    periodos.open("periodos.txt");

    //Introduzco en la variable Np el numero de planetas incluyendo al sol

    Np=9;

    //Defino los valores iniciales de las posiciones y velocidades de cada planeta leyendolos del archivo condicionesiniciales.txt
    //estos valores están en el afelio de cada planeta suponiendo que estos estan alineados


    for(i=0;i<Np; i++){
        
        condini>>p[i][0];
        p[i][1]=0;
        v[i][0]=0;
        condini>>v[i][1];
        a[i][0]=0;
        a[i][1]=0;
        nvueltas[i]=0;
        periodo[i]=0;
        
    }
//valores iniciales de las masas
    masas[0]=2*pow(10,30);
    masas[1]=0.330*pow(10,24);
    masas[2]=4.87*pow(10,24);
    masas[3]=5.97*pow(10,24);
    masas[4]=0.642*pow(10,24);
    masas[5]=1898*pow(10,24);
    masas[6]=568*pow(10,24);
    masas[7]=86.8*pow(10,24);
    masas[8]=102*pow(10,24);

//Cambio las masas a unidades reducidas
    for(i=0; i<Np; i++){
        masas[i]=cambiomasa(masas[i]);
    }


    //Pasamos las distancias a las distancias relativas y las velocidades a las unidades reducidas

    for(i=0;i<Np;i++){
        p[i][0]=cambiodistancia(p[i][0]);
    }

    for(i=0;i<Np;i++){
        v[i][1]=v[i][1]*sqrt(c/(G*Ms));
    }




    //Vamos a calcular la energia inicial de cada planeta
    for(i=1;i<Np;i++){
        modvel=sqrt(pow(v[i][0],2)+pow(v[i][1],2));
        moddist=sqrt(pow(p[i][0],2)+pow(p[i][1],2));
        energiaini[i]=1/2*masas[i]*pow(modvel,2)-G*Ms*masas[i]/(pow(moddist,2));
    }



    //Realizamos 1000 iteraciones del Algoritmo de Verlet, escribiendo en cada iteracion la posicion en el archivo posicionesplanetas.txt
    h=0.01;
    t=0;
    for(n=0;n<=100000;n++){
        for(i=1;i<Np;i++){
        yprev[i]=p[i][1];
        }
       t=t+verlet2d(p,v, a,h, Np);

        //Vamos a calcular los periodos. Para ello utilizaré el vector auxiliar para saber si hay un cambio de signo en la coordenada y
        //Si hay un cambio de signo y la x es positiva entonces sumo 1 vuelta y tomo el tiempo del periodo
        for(i=1;i<Np;i++){
            cambiosigno=yprev[i]*p[i][1];
            if((cambiosigno<0)&&(p[i][0]>0)){
                nvueltas[i]++;
                periodo[i]=t;
            }
        }
        
        if(n%50==0){
           energias<<t<<"\t ";
        for(i=0;i<Np;i++){
        modvel=sqrt(pow(v[i][0],2)+pow(v[i][1],2));
        moddist=sqrt(pow(p[i][0],2)+pow(p[i][1],2));

        energia[i]=c/G*1/2*masas[i]*pow(modvel,2)-G*masas[i]/(moddist);

    //    for(j=0;j<Np;j++){
    //    if(j=!i)
    //    {
    //    moddist=sqrt(pow(p[i][0]-p[j][0],2)+pow(p[i][1]-p[j][1],2));
    //    energia[i]=-G*G*Ms*Ms*masas[i]/(c*c*c*(pow(moddist,2)));
    //    }
     //   }
        posplat<< distancianormal(p[i][0])<<", "<< distancianormal(p[i][1])<<endl;
        posplat<<p[i][0]<< ", "<< p[i][1]<<endl;
        energias<<G/c*energia[i]<<"\t ";
        geocent<<p[i][0]-p[3][0]<< ", "<< p[i][1]-p[3][1]<<endl;

        }
        posplat<<"\n";
        geocent<<"\n";
        energiatot=0;
        for (j =1; j < Np; j++)
        {
            energiatot=energiatot+G/c*energia[i];
        }
        
        energias<<energiatot<<"\n";
        }
        
    }

    //muestro los resultados
    cout<<"N vueltas \t Periodo"<<endl;
    periodos<<"N vueltas \t Periodo"<<endl;
    for(i=1;i<Np;i++){
        cout<<nvueltas[i]<<"\t"<< tiemponormal(periodo[i]/nvueltas[i], c, Ms, G)/86400<<endl;
        periodos<<nvueltas[i]<<"\t"<< tiemponormal(periodo[i]/nvueltas[i], c, Ms, G)/86400<<endl;
    }
    //cierro archivos
    posplat.close();
    condini.close();
    energias.close();
    geocent.close();
    periodos.close();

return 0;

}

float newton(float (&p)[10][2],float (&v)[10][2], float (&a)[10][2], int Np){

    float masas[9], sumac[2];
    float dist;
    int i,j;

    //Introduzco los valores de las masas de los planetas
    masas[0]=2*pow(10,30);
    masas[1]=0.330*pow(10,24);
    masas[2]=4.87*pow(10,24);
    masas[3]=5.97*pow(10,24);
    masas[4]=0.642*pow(10,24);
    masas[5]=1898*pow(10,24);
    masas[6]=568*pow(10,24);
    masas[7]=86.8*pow(10,24);
    masas[8]=102*pow(10,24);

    //Inicializo las variables utilizadas para los sumatorios
    sumac[0]=0.;
    sumac[1]=0.;

    for(i=0; i<Np; i++){
        masas[i]=cambiomasa(masas[i]);
    }
    //Este for es para calcular las aceleraciones de todos los planetas, de forma que recorremos la matriz calculando las aceleraciones
    //que estan en las posiciones 4 y 5 del vector, para ello aplicamos la formula vista en los apuntes
    for(i=1;i<Np; i++){
        //para cada planeta necesito un for que recorra la distancias en el planeta i y el resto j
        for(j=0;j<Np;j++){
            if(i==j){
                //si estamos en en la posicion del mismo planeta entonces nos lo saltamos
                sumac[0]=sumac[0];
                sumac[1]=sumac[1];
            }
            else{
                //Meto aqui la formula que viene en la ecuacion (4) de los apuntes
                dist=pow(sqrt(pow((p[i][0]-p[j][0]),2)+pow((p[i][1]-p[j][1]),2)),3);
                sumac[0]=sumac[0]+masas[j]*(p[i][0]-p[j][0])/dist;
                sumac[1]=sumac[1]+masas[j]*(p[i][1]-p[j][1])/dist;
            }
        }
    //Guardo los valores en la matriz

    a[i][0]=-sumac[0];
    a[i][1]=-sumac[1];

    //Reinicializo los valores
    sumac[0]=0.;
    sumac[1]=0.;

    

    }



    return 0;
}




float verlet2d(float (&p)[10][2],float (&v)[10][2], float (&a)[10][2],float h,int Np){
    float w[10][2];
    int i;
    
    
    //Paso 1 calcular la aceleracion

    newton(p,v,a,Np);


    //Paso 2 calcular la posicion y la w

    for(i=1;i<Np;i++){
        p[i][0]=p[i][0]+h*v[i][0]+h*h/2*a[i][0];
        p[i][1]=p[i][1]+h*v[i][1]+h*h/2*a[i][1];
        w[i][0]=v[i][0]+h/2*a[i][0];
        w[i][1]=v[i][1]+h/2*a[i][1];
    }


    //Paso 3 Reevaluar la aceleracion
    newton(p,v,a,Np);
    //Paso4 Calcular la velocidad
    for(i=1;i<Np;i++){
        v[i][0]=w[i][0]+h/2*a[i][0];
        v[i][1]=w[i][1]+h/2*a[i][1];
    }  

    return h;
}


//funcion que cambia de la masa normal a masa para utilizar 
float cambiomasa(float m){
    float mi, Ms;
    Ms=2*pow(10,30);
    mi=m/Ms;
    return mi;
}


//funcion que pasa de la masa que utilizamos en los calculos a la masa normal
float masanormal(float m){
    float mi, Ms;
    Ms=2*pow(10,30);
    mi=m*Ms;
    return mi;
}

//funcion que pasa de la masa normal a la que utilizamos

float cambiodistancia(float d)
{
    float di, c;
    c=1.496*pow(10,11);
    di=d/c;

    return di;
}

float distancianormal(float d)
{
    float di, c;
    c=1.496*pow(10,11);
    di=d*c;

    return di;
}



float tiemponormal(float t,float c,float G,float Ms){
    float ti;
    ti=t*sqrt(pow(c,3)/(G*Ms));

    return ti;
}