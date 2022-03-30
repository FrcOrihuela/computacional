//Programa que simula el movimiento de un péndulo

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

float newton(float x);
float verlet(float &x, float &v, float &a, float h);


int main (){
    ofstream posiciones, resultados;
    float x, v, w, a, b, h, t;
    int n, i;

//defino los valores iniciales de la posicion y la velocidad (y el tiempo)

    x=-2.5;
    v=0.000;   
    t=0;

//Defino el valor del salto h

    h=0.01;





    cout << "Hola";


    posiciones.open("resultados.txt");
    resultados.open("posyvel.txt");


if(!posiciones.is_open()){
    cout <<"No se ha podido abrir el archivo"<<endl;
}
else{

    posiciones<< cos(x) << ", " << sin(x) << endl;
    resultados<< x<< ", " << v<< endl;
    posiciones<< "\n";

for(i=0; i<=10000; i++)
{
    
    t=t+verlet(x, v, a, h);
    posiciones<< cos(x) << ", " << sin(x) << endl;
    resultados<< x<< ", " << v<< endl;
    resultados<< "\n";
    posiciones<< "\n";

}

}

posiciones.close();

return 0;

}


float newton(float x){
    float a, g;
    g=9.81;

    a=-g*sin(x);
    return a;
}


float verlet(float &x, float &v, float &a, float h){
    float w;
    //Paso 1 calcular la aceleracion
    a=newton(x);
    //Paso 2 calcular la posicion y la w
    x=x+h*v+h*h/2*a;
    w=v+h/2*a;
    //Paso 3 Reevaluar la aceleracion
    a=newton(x);
    //Paso4 Calcular la velocidad
    v=w+h/2*a;

    return h;
}