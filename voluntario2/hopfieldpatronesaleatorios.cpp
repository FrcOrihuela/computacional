#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <ctime>



using namespace std;

//int RandomInt(int min, int max);

static double W[4096][4096];

int main(void){




   //Variables varias
   int x, y, i, j, k, l, m, contador, N, N2, generador, Cp,g;
   bool aleatorio;
   
   

   double deformacion, T, a[50], H[2], deltaH, p, e, ji, solap[50], solapmax;

   //ifstream patron0,patron1,patron2,patron3;
   ofstream sistema, solapamiento, memoriapatrones;

   //patron0.open("patron0.txt");
   //patron1.open("patron1.txt");
   //patron2.open("patron2.txt");
   //patron3.open("patron3.txt");

   sistema.open("randompatrones.txt");
   solapamiento.open("randomsolapamientoS.txt");
   memoriapatrones.open("Memoriapatrones.txt");
   //Numero de Neuronas

   N=64;
   N2=N*N;

   //Declaro la matriz de las interacciones que depende de las coordenadas de las 2 neuronas que interactuan
   
   //Explico aqui el funcionamiento de esta matriz. Como el compilador no me dejaba crear una matriz de dimension 4 he decido
   //Insertar todo en esta matriz. Funciona de la siguiente forma cuando se vea en el programa:
   //i*N es la fila, j la columna de forma que W[N*i+j][k*N+l] es la interacción entre la Neurona (i,j) y la (k,l)

   //La matriz de las neuronas
   int S[N][N];
   //Umbral de disparo
   double disp[N][N];
   //Patron inicial
   int P[50][N2];
   srand( time(NULL) );
   
   
   for(Cp=1;Cp<=50;Cp++){
   memoriapatrones<<Cp<<"\t";
   //Metemos los patrones aleatorios en P
    for(k=0;k<Cp;k++){
        for(i=0;i<N;i++){
            for(j=0;j<N;j++){
                generador=rand()%(1000);;
                if(generador<500){
                    P[k][N*i+j]=0;
                }else{
                    P[k][N*i+j]=1;
                }
            }
        }
   }

   //Calculo el valor de a del patron introducido
   for(m=0;m<Cp;m++){
   a[m]=0;

   for(i=0;i<N;i++){
      for(j=0;j<N;j++){
         a[m]=a[m]+1.0*P[m][i*N+j];
      }
   }
   a[m]=a[m]/N2;
   }

   //Parametro aleatorio? true si aleatorio false si coger un patron

   aleatorio=true;


   //Temperatura?

   T=0.0001;


     
   //Pasamos a calcular el valor de las interacciones entre neuronas

   for(i=0;i<N;i++){
      for(j=0;j<N;j++){
         for(k=0;k<N;k++){
            for(l=0;l<N;l++){
               W[N*i+j][k*N+l]=0;
               if ((i==k)&&(j==l)){
                  W[N*i+j][k*N+l]=0;
               }else{
                  for(m=0;m<Cp;m++){
                  W[N*i+j][k*N+l]=W[N*i+j][k*N+l]+1.0/N2*(1.0*P[m][i*N+j]-a[m])*(1.0*P[m][k*N+l]-a[m]);
                  }
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

//Genero estado inicial
for(g=0;g<25;g++){
   if(aleatorio==true){
      for(i=0;i<N;i++){
         for(j=0;j<N;j++){
            generador=rand()%(1000);;
            if(generador<500){
               S[i][j]=0;
            }else{
               S[i][j]=1;
            }
         }
      }
   }else{
      //Parametro de deformacion
      /*deformacion=0.2;
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
      }*/
   }


   //Metemos el estado inicial en el fichero de estados

   for(i=0;i<N;i++){
      for(j=0;j<N-1;j++){
         sistema<<S[i][j]<<", ";
         sistema.flush();
      }
      sistema<<S[i][N-1]<<"\n";
      
   }
   sistema<<"\n";
   sistema.flush();



   contador=0;

if(sistema.is_open()){
//Comienza el proceso
   for(m=0;m<100*N2;m++){
      //Genero las 2 coordenadas aleatorias
      x=rand()%(N);
      y=rand()%(N);

      if(m%N2 ==0){
      for(k=0;k<Cp;k++){
         solap[k]=0;
         for(j=0;j<N;j++){
            for(i=0;i<N-1;i++){
               solap[k]=solap[k]+(P[k][j*N+i]-a[k])*(S[j][i]-a[k]);
            }
         }
         
      solap[k]=solap[k]/(N2*a[k]*(1-a[k]));
      }

         for(j=0;j<N;j++){
               for(i=0;i<N-1;i++){
                  sistema<<S[j][i]<<", ";
                  sistema.flush();
                  
               }
               sistema<<S[j][N-1]<<"\n";
            }sistema<<"\n";

      sistema.flush();

      
      solapamiento<<m<<"\t"<<solap[0]<<"\t"<<solap[1]<<"\t"<<solap[2]<<"\t"<<solap[3]<<"\n";
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
solapmax=abs(solap[0]);
for(i=0;i<Cp;i++){
   if(solapmax<abs(solap[i])){
   solapmax=abs(solap[i]);
   }
}
memoriapatrones<<solapmax<<"\t";
memoriapatrones.flush();

}
memoriapatrones<<"\n";
memoriapatrones.flush();
}
   /*patron0.close();
   patron1.close();
   patron2.close();
   patron3.close();*/
   sistema.close();


   return 0;

}


/*
int RandomInt(int min, int max){
   return min + rand() % (( max + 1 ) - min);
}
*/