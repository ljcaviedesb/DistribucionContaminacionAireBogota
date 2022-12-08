#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>
using namespace std;

//-----Constantes Globales-----//
const int Lx = 180;                 //tamaño de la simulacion
const int Ly = 400;

const int Q = 5;                    //numero de direcciones
const double W0 = 1.0/3;            //pesos
const int G = 2;                    //Numero de picos gaussianos

const double C = 0.5;               //velocidad de la onda - por estabilidad numerica: C < 0.707 cells/click
const double C2 = C*C;
//const double AUX0 = 1-3*C2*(1-W0);  //la parte i = 0 de la funcion de equilibrio

const double tau = 2.5;             //valores en la funcion de colision/ para obtener un D = 0.05
const double Utau = 1.0/tau;
const double UmUtau = 1-Utau;

//-----Clase LatticeBoltzmann------
class LatticeBoltzmann{
private:
  double w[Q];        //pesos por direccion
  int Vx[Q], Vy[Q];   //vectores de velocidad
  double M[2][2*G];     //Multivector, que guarda opsición y velocidad de inicio de los picos 
  double *f, *fnew;   //funciones de distribucion - * es un apuntador que lo lleva a la direccion de memoria
public:
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  int n(int ix, int iy, int i){return (ix*Ly+iy)*Q+i;};   //indice de las celdas
  double rho(int ix, int iy, bool UseNew);
  double Jx(int ix, int iy, bool UseNew);
  double Jy(int ix, int iy, bool UseNew);
  double feq(double rho0, double Ux0, double Uy0, int i);
  void Start(double rho0, double Ux0, double Uy0);
  void Collision(void);
  void ImposeFields(int t,double Ux0,double Uy0);
  void Advection(void);
  void Print(const char * NameFile);
};

//Constructor
LatticeBoltzmann::LatticeBoltzmann(void){
  //Pesos
  w[0] = W0; w[1] = w[2] = w[3] = w[4] = (1.0-W0)/4;
  //Vectores de velocidad
  Vx[0] = 0;  Vx[1] = 1;  Vx[2] = 0;  Vx[3] = -1; Vx[4] = 0;
  Vy[0] = 0;  Vy[1] = 0;  Vy[2] = 1;  Vy[3] = 0;  Vy[4] = -1;
  //Posiciones inicio de los picos
  M[0][0]=Lx/4; M[0][1]=Ly/4; //Representa el punto inicial en (45,100) para el primer pico
  M[0][2]=3*Lx/4; M[0][3]=3*Ly/4;  //punto inicial en (135,300) para el segundo pico
  //velocidades inicio de los picos
  M[1][0]=1; M[1][1]=1; //Coef de multiplicación de las velocidades U0x y U0y para el primer pico
  M[1][2]=-1; M[1][3]=-1;
  //Crea arrays dinamicos
  int ArraySize = Lx*Ly*Q;
  f = new double [ArraySize];  fnew = new double [ArraySize];
}
//Destructor
LatticeBoltzmann::~LatticeBoltzmann(void){
  delete[] f;  delete[] fnew;
}
//Rho
double LatticeBoltzmann::rho(int ix, int iy, bool UseNew){
  double sum; int i, n0;
  for(sum=0, i=0; i<Q; i++){
    n0 = n(ix,iy,i);
    if(UseNew) sum += fnew[n0]; else sum += f[n0];
  }
  return sum;
}
//Componentes del vector J
double LatticeBoltzmann::Jx(int ix, int iy, bool UseNew){
  double sum; int i, n0;
  for(sum=0, i=0; i<Q; i++){
    n0 = n(ix,iy,i);
    if(UseNew) sum += Vx[i]*fnew[n0]; else sum += Vx[i]*f[n0];
  }
  return sum;
}
double LatticeBoltzmann::Jy(int ix, int iy, bool UseNew){
  double sum; int i, n0;
  for(sum=0, i=0; i<Q; i++){
    n0 = n(ix,iy,i);
    if(UseNew) sum += Vy[i]*fnew[n0]; else sum += Vy[i]*f[n0];
  }
  return sum;
}
//Funcion equilibrio
double  LatticeBoltzmann::feq(double rho0, double Ux0, double Uy0, int i){
  double UdotVi = Ux0*Vx[i]+Uy0*Vy[i], U2 = Ux0*Ux0+Uy0*Uy0;
  return rho0*w[i]*(1+(UdotVi/C2)+(pow(UdotVi,2)/(2*pow(C2,2)))-(U2/(2.0*C2)));
}
//Start
void LatticeBoltzmann::Start(double rho0, double Ux0, double Uy0){
  int ix,iy,i,n0;
  for(ix=0; ix<Lx; ix++)    //para cada celda
    for(iy=0; iy<Ly; iy++)
      for(i=0; i<Q; i++){   //en cada direccion
        n0 = n(ix,iy,i);
        f[n0] = feq(rho0,Ux0,Uy0,i);
      }
}
//Colision
void LatticeBoltzmann::Collision(void){
  int ix,iy,i,n0; double rho0,Ux0,Uy0;
  for(ix=0;ix<Lx;ix++)      //para cada celda
    for(iy=0;iy<Ly;iy++){
      //Calcule los campos macroscopicos en la celda
      rho0=rho(ix,iy,false); Ux0=Jx(ix,iy,false)/rho0; Uy0=Jy(ix,iy,false)/rho0;
      for(i=0;i<Q;i++){     //para cada vector de velocidad
        n0 = n(ix,iy,i);
        fnew[n0] = UmUtau*f[n0]+Utau*feq(rho0,Ux0,Uy0,i);
      }
    }
}
//Imponer campos
void LatticeBoltzmann::ImposeFields(int t,double Ux0,double Uy0){
  int i, ix, iy, ix0, iy0,ix01, iy01, n0;
  double rho_gauss, x0_bar, y0_bar,x01_bar, y01_bar;
  double xbar[G], ybar[G];  //vectores posición bar
  //Constantes Gaussian Hill
  double A = 1.0, D = C2*(tau-0.5), Sigma = 5.0, Sigma2 = Sigma*Sigma;
  double Coef = A/(1.0+(2.0*D*t/Sigma2));
  //Fuente gaussiana en el medio
  /*ix0 = Lx/4; iy0 = Ly/4; //punto inicial en (45,100) 
  ix01 = 3*Lx/4; iy01 = 3*Ly/4; //punto inicial en (135,300)*/
  //x0_bar=M[0][0]+Ux0*M[1][0]*t; y0_bar=M[0][1]+M[1][1]*Uy0*t;
  //x01_bar=M[0][2]+Ux0*M[1][2]*t; y01_bar=M[0][3]+M[1][3]*Uy0*t;
  //inicialización posiciones bar
  for(int j=0;j<G;j++){
  xbar[j]=M[0][2*j]+Ux0*M[1][2*j]*t;
  ybar[j]=M[0][(2*j)+1]+M[1][(2*j)+1]*Uy0*t;
  }
  for(ix=0; ix<Lx; ix++){      //para cada celda
    for(iy=0; iy<Ly; iy++){
      //rho0 = rho(ix, iy, false); Ux0=Jx(ix,iy,false)/rho0; Uy0=Jy(ix,iy,false)/rho0;
      //rho_gauss = 0.05+Coef*exp(-(pow(ix-xbar[0],2)+pow(iy-ybar[0],2))/(2.0*(Sigma2+2.0*D*t)))+Coef*exp(-(pow(ix-xbar[1],2)+pow(iy-ybar[1],2))/(2.0*(Sigma2+2.0*D*t))); //K=0.05 es usada para evitar NaN.
      rho_gauss=0;
      for(int k=0;k<G;k++){
      rho_gauss+=Coef*exp(-(pow(ix-xbar[k],2)+pow(iy-ybar[k],2))/(2.0*(Sigma2+2.0*D*t)));
      }
      //Ux0=Jx(ix,iy,false)/rho_gauss; Uy0=Jy(ix,iy,false)/rho_gauss;
      for(i=0; i<Q; i++){
        n0=n(ix,iy,i);
        fnew[n0]=feq(rho_gauss,Ux0,Uy0,i);
      }
    }
  }
}
//Adveccion
void LatticeBoltzmann::Advection(void){
  int ix,iy,i,ixnext,iynext,n0,n0next;
  for(ix=0;ix<Lx;ix++)      //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){     //en cada direccion
        ixnext=(ix+Vx[i]+Lx)%Lx; iynext=(iy+Vy[i]+Ly)%Ly;
        n0=n(ix,iy,i); n0next=n(ixnext,iynext,i);
        f[n0next]=fnew[n0];     //fronteras periodicas
      }
}
//Print
void LatticeBoltzmann::Print(const char * NameFile){
  ofstream MyFile(NameFile); double rho0; int ix, iy;
  for(ix=0; ix<Lx; ix++){
    for(iy=0; iy<Ly; iy++){
      rho0=rho(ix,iy,true);
      MyFile<<ix<<" "<<iy<<" "<<rho0<<endl;
    }
    MyFile<<endl;
  }
  MyFile.close();
}

//-----Programa Principal-----

int main(void){
  LatticeBoltzmann Ondas;
  int t, tmax = 4;
  double rho0 = 1.0, Ux0 = 10.5, Uy0 = 10.5;

  Ondas.Start(rho0, Ux0, Uy0);
  //Evolucione
  for(t=0; t<tmax; t++){
    Ondas.Collision();
    Ondas.ImposeFields(t,Ux0,Uy0);
    Ondas.Advection();
  }
  //Imprima
  Ondas.Print("Ondas.dat");

  return 0;
}