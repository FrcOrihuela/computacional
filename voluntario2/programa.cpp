#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <ctime>



using namespace std;

//int RandomInt(int min, int max);


int main(void){




   //Variables varias
   float deformacion, T, a, H[2], deltaH, p, e, ji;
   //float solapamiento;
   int x, y, i, j, k, l, m, contador, N, N2, generador;
   bool aleatorio;
   
   ifstream condini;
   ofstream estados;

   condini.open("condini.txt");
   estados.open("evolucionsist.txt");


   //Numero de Neuronas

   N=20;
   N2=N*N;

   //Declaro la matriz de las interacciones que depende de las coordenadas de las 2 neuronas que interactuan
   float W[N2][N2];
   //Explico aqui el funcionamiento de esta matriz. Como el compilador no me dejaba crear una matriz de dimension 4 he decido
   //Insertar todo en esta matriz. Funciona de la siguiente forma cuando se vea en el programa:
   //i*N es la fila, j la columna de forma que W[N*i+j][k*N+l] es la interacción entre la Neurona (i,j) y la (k,l)

   //La matriz de las neuronas
   int S[50][50];
   //Umbral de disparo
   float disp[50][50];
   //Patron inicial
   int P[50][50];

   //Metemos el patron inicial en P
   for(i=0;i<N;i++){
      for(j=0;j<N;j++){
         condini>>P[i][j];
      }
   }

   //Parametro aleatorio? true si aleatorio false si coger un patron

   aleatorio=false;
   srand( time(NULL) );

   //Temperatura?

   T=0.0001;

   if(aleatorio==true){
      for(i=0;i<N;i++){
         for(j=0;j<N;j++){
            generador=rand()%(1001);;
            if(generador<500){
               S[i][j]=-1;
            }else{
               S[i][j]=1;
            }
         }
      }
   }else{
      //Parametro de deformacion
      deformacion=0.01;
      for(i=0;i<N;i++){
         for(j=0;j<N;j++){
            S[i][j]=P[i][j];
         }
      }
      //el parametro de deformacion me dice el porcentaje de cuadros a los que les cambio el signo

      for ( i = 0; i < trunc(N2*deformacion); i++){
         x=rand()%(N+1);
         y=rand()%(N+1);
         //Si es 0 me da 1, si es 1 me da 0
         S[x][y]=1-S[x][y];
      }
   }


   //Metemos el estado inicial en el fichero de estados

/*   for(i=0;i<N;i++){
      for(j=0;j<N-1;j++){
         estados<<P[i][j]<<", ";
      }
      estados<<P[i][N-1]<<"\n";
   }

*/

   //Calculo el valor de a del patron introducido

   a=0;

   for(i=0;i<N;i++){
      for(j=0;j<N;j++){
         a=a+P[i][j];
      }
   }
   a=a/N2;

   //Pasamos a calcular el valor de las interacciones entre neuronas

   for(i=0;i<N;i++){
      for(j=0;j<N;j++){
         for(k=0;k<N;k++){
            for(l=0;l<N;l++){
               if (i==k&&j==l){
                  W[N*i+j][k*N+l]=0;
               }else{
                  W[N*i+j][k*N+l]=1/N2*(P[i][j]-a)*(P[k][l]-a);
               }
            }
         }
      }
   }

   //Y ahora calculamos el termino del umbral de disparo, primero inicializo los valores y hago sumatorios

   for(i=0;i<N;i++){
      for(j=0;j<N;j++){
         disp[i][j]=0;
         for(k=0;k<N;k++){
            for(l=0;l<N;l++){
               disp[i][j]=disp[i][j]+W[i*N+j][k*N+l];
            }
         }
         disp[i][j]=disp[i][j]/2;
      }
   }

   contador=0;


//Comienza el proceso
   for(m=0;m<100*N2;m++){
      //Genero las 2 coordenadas aleatorias
      x=rand()%(N+1);
      y=rand()%(N+1);

      if(m%N2 ==0){
         for(j=0;j<N;j++){
            for(i=0;i<N-1;i++){
               estados<<S[j][i]<<", ";
               }
            estados<<S[j][N-1]<<"\n";
         }estados<<"\n";
         contador++;
      }


      for(k=0;k<N;k++){
         for(l=0;l<N;l++){
            deltaH=deltaH+(W[x*N+y][k*N+l]*S[k][l]);
         }
      }
      deltaH=(1-2*S[x][y])*(disp[x][y])*(1/2+S[x][y]-deltaH);

      e=exp(-deltaH/T);
      if(e>1){
         p=1;
      }else{
         p=e;
      }
      
      //Ahora saco un numero aleatorio entre 0 y 1
      ji=1.0*(rand()%RAND_MAX)/RAND_MAX;

      if(ji<p){
         S[x][y]=1-S[x][y];
      }

   }

   condini.close();
   estados.close();


   return 0;
}


/*
int RandomInt(int min, int max){
   return min + rand() % (( max + 1 ) - min);
}
*/