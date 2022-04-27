#include <iostream>
#include <cmath>
#include <fstream>
#include <cstdlib>
#include <ctime>



using namespace std;

float RandomFloat(float a, float b) ;
int RandomInt(int min, int max) ;


int main(void){
    //Declaro la matriz de spines

    int S[200][200];
    double T, p, deltaE, ji,e;
    int i,j,k,N, x, y, N2;
    int generador;
    ofstream sistema;
    sistema.open("evolucionsistema.txt");
    bool aleatorio;
    //Orden de la matriz
    N=64;
    N2=N*N;
    
    //Temperatura inicial
    T=0.001;

    //Semilla de numeros aleatorios
    srand( time(NULL) );

    //genero las condiciones iniciales de forma aleatoria, genero un numero entre -1 y 1 y si es negativo le asigno el valor -1 y si es positivo +1
    aleatorio=true;

    if(aleatorio==true){
        for(i=1;i<=N;i++){
        for(j=1;j<=N;j++){
            generador=RandomInt(0,1000);
            if(generador<500){
                S[i][j]=-1;
            }else{
                S[i][j]=1;
            }

        }
    }
    }else{
        for(i=1;i<=N;i++){
        for(j=1;j<=N;j++){
            S[i][j]=1;
            }

        }
    }
    
    //Condiciones de contorno periodicas

    for(i=1;i<=N;i++){
        S[0][i]=S[N][i];
        S[N+1][i]=S[0][i];
        S[i][0]=S[i][N];
        S[i][N+1]=S[i][1];
    }
    for(j=1;j<=N;j++){
                for(i=1;i<=N-1;i++){
                    sistema<<S[j][i]<<", ";
                }
                sistema<<S[j][N]<<"\n";
            }sistema<<"\n";

    for(k=0;k<300*N2;k++){
        //Genero las 2 coordenadas aleatorias
        x=1+rand()%(N);
        y=1+rand()%(N);

        if(k%N2 ==0){
            for(j=1;j<=N;j++){
                for(i=1;i<=N-1;i++){
                    sistema<<S[j][i]<<", ";
                }
                sistema<<S[j][N]<<"\n";
            }sistema<<"\n";
        }


        deltaE=0;
        //calculo deltaE hay que usar condiciones de contorno periodicas

       /* if(x=!N){
            if (x=!0){
                if (y=!N){
                    if(y=!0){
                        deltaE=2*S[x][y]*(S[x+1][y]+S[x-1][y]+S[x][y+1]+S[x][y-1]);
                    }else{
                        deltaE=2*S[x][y]*(S[x+1][y]+S[x-1][y]+S[x][y+1]+S[x][N]);
                    }
                }else{
                    deltaE=2*S[x][y]*(S[x+1][y]+S[x-1][y]+S[x][0]+S[x][y-1]);
                }
                
            }else{
                if (y=!N){
                    if(y=!0){
                        deltaE=2*S[x][y]*(S[x+1][y]+S[N][y]+S[x][y+1]+S[x][y-1]);
                    }else{
                        deltaE=2*S[x][y]*(S[x+1][y]+S[N][y]+S[x][y+1]+S[x][N]);
                    }
                }else{
                    deltaE=2*S[x][y]*(S[x+1][y]+S[N][y]+S[x][0]+S[x][y-1]);
                    }
                
                }    
        }else{
            if(y=!0){
                if(y=!N){
                    deltaE=2*S[x][y]*(S[0][y]+S[x-1][y]+S[x][y+1]+S[x][y-1]);
                }else{
                    deltaE=2*S[x][y]*(S[0][y]+S[x-1][y]+S[x][0]+S[x][y-1]);
                }
            }else{
                deltaE=2*S[x][y]*(S[0][y]+S[x-1][y]+S[x][y+1]+S[x][N]);
            }
        }*/
        deltaE=2*S[x][y]*(S[x+1][y]+S[x-1][y]+S[x][y+1]+S[x][y-1]);

        e=exp(-deltaE/T);
        //Evaluo p como el minimo entre 1 y exp(-deltaE/T)

        if(e>1){
            p=1;
        }else{
            p=exp(-deltaE/T);
        }


        //Ahora saco un numero aleatorio entre 0 y 1
        ji=1.0*(rand()%RAND_MAX)/RAND_MAX;

        if(ji<p){
            S[x][y]=-S[x][y];
        }
    //Condiciones de contorno periodicas

    for(i=1;i<=N;i++){
        S[0][i]=S[N][i];
        S[N+1][i]=S[0][i];
        S[i][0]=S[i][N];
        S[i][N+1]=S[i][1];
    }

    }





return 0;
}



float RandomFloat(float a, float b) {
    float random = ((float) rand()) / (float) RAND_MAX;
    float diff = b - a;
    float r = random * diff;
    return a + r;
}


int RandomInt(int min, int max) //rango : [min, max]
{
   return min + rand() % (( max + 1 ) - min);
}