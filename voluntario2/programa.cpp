#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <vector>
#include <array>


using namespace std;



//Explico aqui el funcionamiento de esta matriz. Como el compilador no me dejaba crear una matriz de dimension 4 he decido
//Insertar todo en esta matriz. Funciona de la siguiente forma cuando se vea en el programa:
//i*N es la fila, j la columna de forma que W[N*i+j][k*N+l] es la interacción entre la Neurona (i,j) y la (k,l)

static float W[4096][4096];

int main(void){




   //Variables varias
   double deformacion, T, a, deltaH, p, e, ji, solap;
   //float solapamiento;
   int x, y, i, j, k, l, m, contador, N, N2, generador;
   bool aleatorio;
   
   ifstream condini;
   ofstream sistema, solapamiento, dependetemp;

   condini.open("sans64.txt");
   sistema.open("evolucionpatron.txt");
   solapamiento.open("solapamiento.txt");
   dependetemp.open("dependenciatemp.txt");
   //Numero de Neuronas

   N=64;
   N2=N*N;

   //Declaro la matriz de las interacciones que depende de las coordenadas de las 2 neuronas que interactuan

 
   //La matriz de las neuronas
   short int S[N][N];
   //Umbral de disparo
   float disp[N][N];
   //Patron inicial
   short int P[N][N];

   //Metemos el patron inicial en P
   for(i=0;i<N;i++){
      for(j=0;j<N;j++){
         condini>>P[i][j];
      }
   }

   //Parametro aleatorio? true si aleatorio false si coger un patron

   aleatorio=true;
   srand( time(NULL) );

   //Temperatura?

   T=0.0001;


  // while(T<0.1){

   if(aleatorio==true){
      for(i=0;i<N;i++){
         for(j=0;j<N;j++){
            generador=rand()%(1001);;
            if(generador<500){
               S[i][j]=0;
            }else{
               S[i][j]=1;
            }
         }
      }
   }else{
      //Parametro de deformacion
      deformacion=0.8;
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

   for(i=0;i<N;i++){
      for(j=0;j<N-1;j++){
         sistema<<P[i][j]<<", ";
         sistema.flush();
      }
      sistema<<P[i][N-1]<<"\n";
      
   }
   sistema<<"\n";
   sistema.flush();


   //Calculo el valor de a del patron introducido

   a=0;

   for(i=0;i<N;i++){
      for(j=0;j<N;j++){
         a=a+1.0*P[i][j];
      }
   }
   a=a/N2;

   //Pasamos a calcular el valor de las interacciones entre neuronas

   for(i=0;i<N;i++){
      for(j=0;j<N;j++){
         for(k=0;k<N;k++){
            for(l=0;l<N;l++){
               if ((i==k)&&(j==l)){
                  W[N*i+j][k*N+l]=0;
               }else{
                  W[N*i+j][k*N+l]=1.0/N2*(1.0*P[i][j]-a)*(1.0*P[k][l]-a);
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
               disp[i][j]=disp[i][j]+0.5*W[i*N+j][k*N+l];
            }
         }
      }
   }

   contador=0;

if(sistema.is_open()){
//Comienza el proceso
   for(m=0;m<100*N2;m++){
      //Genero las 2 coordenadas aleatorias
      x=rand()%(N);
      y=rand()%(N);

      if(m%N2 ==0){
      solap=0;

         for(j=0;j<N;j++){
               for(i=0;i<N-1;i++){
                  sistema<<S[j][i]<<", ";
                  sistema.flush();
                  solap=solap+(P[j][i]-a)*(S[j][i]-a);
               }
               sistema<<S[j][N-1]<<"\n";
            }sistema<<"\n";
      solap=solap/(N2*a*(1-a));
      sistema.flush();

      
      solapamiento<<m<<"\t"<<solap<<"\n";
      solapamiento.flush();
      }
      deltaH=0;

      //Valor de deltaH obtenido a partir de bibliografia, consultar
      for(k=0;k<N;k++){
         for(l=0;l<N;l++){
            deltaH=deltaH+(W[x*N+y][k*N+l]*S[k][l]);
         }
      }
      deltaH=(-deltaH/2+disp[x][y])*(1-2*S[x][y]);

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
}
//dependetemp<<log10(T)<< "\t" <<solap<<"\n";
//dependetemp.flush();
//T=T+0.0005;
//}
   condini.close();
   sistema.close();
   solapamiento.close();
   dependetemp.close();

   return 0;
}


/*
int RandomInt(int min, int max){
   return min + rand() % (( max + 1 ) - min);
}
*/