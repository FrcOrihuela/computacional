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

    int S[100][100];
    double T, p, deltaE, ji,e;
    // Las siguientes son las medias de magnetizacion energia
    double Mmag, Mene, Desvene, Energia[100][100],cn,en,mn,f[100];
    int i,j,k,N, x, y, N2,l,m,medidas;
    int generador;
    ofstream sistema;
    sistema.open("evolucionsistema.txt");
    
    //Orden de la matriz
    N=16;
    N2=N*N;
    
    //Temperatura inicial
    T=0.5;

    //Semilla de numeros aleatorios
    srand( time(NULL) );

    //genero las condiciones iniciales de forma aleatoria, genero un numero entre -1 y 1 y si es negativo le asigno el valor -1 y si es positivo +1
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


    Mmag=0;
    Mene=0;
    Desvene=0;
    medidas=0;

    for(k=0;k<1000000*N2;k++){
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


    //Suma para la magnetización y energía media
    

    if(k%(100*N2)==0){
        for(i=1;i<=N;i++){
            for(j=1;j<=N; j++){
                Mmag=Mmag+S[i][j];
                Energia[i][j]=-1/2*(S[i][j])*(S[i][j+1]+S[i][j-1]+S[i+1][j]+S[i-1][j]);
                Mene=Mene+Energia[i][j];
                Desvene=Desvene+Energia[i][j]*Energia[i][j];

                medidas++;
            }
            
        }

        for(i=0;i<N;i++){
            for(j=1;j<=N;j++){
                for(k=1;k<=N;k++){
                    if(i+j>N){
                        f[i]=f[i]+S[i][k]*S[i+j-N][k];
                    }else{
                        f[i]=f[i]+S[j][k]*S[i+j][k];
                    }
                }
            }
            
        }

        Mene=Mene/N2;
        Mmag=Mmag/N2;
        Desvene=Desvene/N2;
    }
    }

mn=Mmag/medidas;
en=Mene/(medidas*2*N2);
cn=1/(N2*T)*(Desvene/medidas+(Mene/medidas)*(Mene/medidas));
for(i=1;i<=N;i++){
    f[i]=f[i]/(N2*N2);
}

cout<<"Magnetizacion promedio: "<< mn<<endl;
cout<<"Energia media: "<< en<<endl;
cout<<"Calor especifico "<< cn<<endl;
cout<<"Funcion de correlacion"<<endl;
for(i=0;i<N;i++){
    cout<<"i="<<i<<"=>"<<f[i]<<endl;
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