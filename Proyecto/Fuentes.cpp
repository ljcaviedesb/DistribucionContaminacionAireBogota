/*
  Proyecto Final:
  Lattice Boltzmann para la ecuaci�n AD en coordenadas cartesianas con una fuente
*/

#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

//-----Constantes Globales-----//
const int Lx = 180;                 //tama�o de la simulacion
const int Ly = 400;

const int Q = 5;                    //numero de direcciones
const double W0 = 1.0/3.0;          //cte que define los pesos

const double C = 0.5;               //velocidad de la onda - por estabilidad numerica: C < 0.707 cells/click
const double C2 = C*C;

const double CoefDiff = 0.05;
const double tau = 0.7;             //Tesis de Juliana
const double Utau = 1.0/tau;
const double UmUtau = 1.0-Utau;

//-----Clase LatticeBoltzmann------
class LatticeBoltzmann{
private:
  double w[Q];        //pesos por direccion
  int Vx[Q], Vy[Q];   //vectores de velocidad
  double *f, *fnew;   //funciones de distribucion - * es un apuntador que lo lleva a la direccion de memoria
  double Snew[Lx][Ly][Q], Sold[Lx][Ly][Q]; 
public:
  LatticeBoltzmann(void);
  ~LatticeBoltzmann(void);
  int n(int ix, int iy, int i){return (ix*Ly+iy)*Q+i;};         //indice de las celdas
  double rho(int ix, int iy, bool UseNew);
  double Jx(int ix, int iy, bool UseNew);
  double Jy(int ix, int iy, bool UseNew);
  double S(int ix, int iy, int t);
  double Si(int ix, int iy, double Ux0, double Uy0, int t, int i);           //forzamiento
  double feq(double rho0, double Ux0, double Uy0, int i);
  void Start(double rho0, double Ux0, double Uy0);
  void Collision(double Ux0, double Uy0);
  void ImposeFields();
  void Advection(double Ux0, double Uy0, int t);
  void Print(const char * NameFile);
  double Varianza(void);
};
//Constructor
LatticeBoltzmann::LatticeBoltzmann(void){
  //Pesos
  w[0] = W0; w[1] = w[2] = w[3] = w[4] = (1.0-W0)/4.0;
  //Vectores de velocidad
  Vx[0] = 0;  Vx[1] = 1;  Vx[2] = 0;  Vx[3] = -1; Vx[4] = 0;
  Vy[0] = 0;  Vy[1] = 0;  Vy[2] = 1;  Vy[3] = 0;  Vy[4] = -1;
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

//FUENTES
double LatticeBoltzmann::S(int ix, int iy, int t){
  
 if(iy==60 && (ix >45 && ix<90)){
    return 0.5; // Auto Sur
  }
  if(iy==120 && (ix >0 && ix<132)){
    return 0.5;  // Cll 13
  }
  if(ix==55  && (iy>0 && iy<250)){
    return 0.5;  // Av Boyaca
  }
  if(iy==48 && (ix >130 && ix<180)){
    return 0.3; // cra 16 (parece mas calle)
  }
  if(iy==100 && (ix >130 && ix<180)){
    return 0.4; // cll 1
  }
  if(iy==160 && (ix >20 && ix<180)){
    return 0.3  ; // Av cll 26
  }
  if(iy==250 && (ix >80 && ix<180)){
    return 0.4;   // Cll 80
  }
  if(ix==150 && (iy >120 && iy<170)){
    return 0.5;  // cra 30
  }
  if(ix==80 && (iy >45 && iy<310)){
    return 0.5; // cra 68
  }
  if(iy==250 && (ix >80 && ix<180)){
    return 0.3; // cll 100
  }
  if(ix==90 && (iy >0 && iy<60)){
    return 0.4; //cll 51 sur
  }
  if(ix==130 && (iy >10 && ix<100)){
    return 0.4; // cra 10
  }
  if(ix==140 && (iy >120 && iy<200)){
    return 0.3;  //caracas
  }
  if(ix==145 && (iy >250 && ix<400)){
    return 0.4;     //autonorte
  }
  if(ix==30 && (iy >270 && iy<400)){
    return 0.5;     // Cali
  }
  if(iy==100 && (ix >10 && ix<150)){
    return 0.3; //suba
  }
  
  
  //Fabricas
  if((ix>=138 && ix<=140) && (iy>=270 && iy<=272)){
    return 0.5;            // Bavaria cra 53 # 127
  }
  if((iy>=177 && iy<=179) && (ix>=50 && ix<=52)){
    return 0.5;            // Industria Nacional de Gaseosas cll 25 # 95
  }
  if((iy>=100 && iy<=102) && (ix>=70 && ix<=72)){
    return 0.5;            // General motors cll 56sur # 36
  }
  if((iy>=150 && iy<=152) && (ix>=160 && ix<=162)){
    return 0.5;            // Diana cra 13 #93
  }
  if((iy>=245 && iy<=247) && (ix>=145 && ix<=147)){
    return 0.5;            // Nestle diag 92 # (cra)19
  }
  if((ix>=90 && ix<=92) && (iy>=157 && iy<=159)){
    return 0.5;            // Cavisan Grupo Sas Calle 2A # 53
  }
  if((ix>=29 && ix<=31) && (iy>=154 && iy<=156)){
    return 0.5;            //Mezcladores Industriales Cra 2 # 5
  }
  if((ix>=38 && ix<=40) && (iy>=233 && iy<=235)){
    return 0.5;            // Fabrica Maquinarias  cl 64 # 110 38,233
  }
  if((ix>=134 && ix<=136) && (iy>=369 && iy<=371)){
    return 0.5;            // Fabrica Productos Limpieza  cra 12-187 134,369
  }
  
  
  //FUENTES INDUSTRIALES 
  
  //TUNJUELITO
  if((ix>55 && iy>0) && (ix<90 && iy<60)){
    return 1;
  }
  //PUENTE ARANDA
  if((ix>70 && iy>120) && (ix<110 && iy<170)){
    return 0.49; //16.51
  }
  //Kenedy
  if((ix>0 && iy>60) && (ix<80 && iy<120)){
    return 0.61; //20.60
  }
  //Fontibon
  if((ix>0 && iy>120) && (ix<80 && iy<160)){
    return 0.57;//19.05
  }
  //Engativa
  if((ix>0 && iy>160) && (ix<80 && iy<250)){
    return 0.46; //15.38
  }
  //Usaquen 
  if((ix>145 && iy>250) && (ix<180 && iy<400)){
    return 0.37; //12.49
  }
 //Suba
 if((ix>0 && iy>250) && (ix<145 && iy<400)){
   return 0.45; //15.33
 }
 // Colina
 if((ix>0 && iy>250) && (ix<145 && iy<400)){
   return 0.45; //15.33
 }
 // cIUDAD bOLIVAR
 if((ix>0 && iy>0) && (ix<55 && iy<60)){
   return 0.57; //19.30
  }
 // Carvajal -Sevillana
 if((ix>040 && iy>50) && (ix<65 && iy<70)){
   return 1; //33.40
  }
 // San Cristobal
 if (ix>130 && ix <180 && iy>48 && ix<100){
   return 0.42;
 }
 // Colina
 if (ix>100 && ix<120 && iy>270 && iy<340){
   return 0.32; //10.94
 }
 // Barrios Unidos
 if(ix>80 && ix<145 && ix>250 && ix<250){
   return 0.42; //14.92
 }
 // Chapinero
 if(ix>140 && ix<180 && iy>140 && iy<250){
   return 0.54; //18.31
 }
 // Centro
 if(ix>140 && ix<180 && iy>100 && iy<160){
   return 0.42; //14.92
 }
 // Tunal
 if(ix>90 && ix<120 && iy>50 && iy<100){
   return 0.58; //19.39
 }
  
 else
   return 0;
} 




//Forzamiento LBGK
double LatticeBoltzmann::Si(int ix, int iy, double Ux0, double Uy0, int t, int i){
  double UdotVi = Ux0*Vx[i]+Uy0*Vy[i];
  return w[i]*S(ix,iy,t)*(1+(UdotVi/C2));
}
//Funcion equilibrio
double  LatticeBoltzmann::feq(double rho0, double Ux0, double Uy0, int i){
  double UdotVi = Ux0*Vx[i]+Uy0*Vy[i], U2 = Ux0*Ux0+Uy0*Uy0;
  return rho0*w[i]*(1.0+(UdotVi/C2)+(pow(UdotVi,2)/(2*pow(C2,2)))-(U2/(2*C2)));  
} 
//Start
void LatticeBoltzmann::Start(double rho0, double Ux0, double Uy0){
  int ix, iy, i, n0; 
  for(ix=0; ix<Lx; ix++)    //para cada celda
    for(iy=0; iy<Ly; iy++)
      for(i=0; i<Q; i++){   //en cada direccion
        n0 = n(ix,iy,i);
        f[n0] = feq(rho0,Ux0,Uy0,i);
        Sold[ix][iy][i] = Si(ix,iy,Ux0,Uy0,0,i);
	Snew[ix][iy][i] = Sold[ix][iy][i];
      }
}  
//Colision
void LatticeBoltzmann::Collision(double Ux0, double Uy0){
  int ix, iy, i, n0; double rho0;
  for(ix=0;ix<Lx;ix++)      //para cada celda
    for(iy=0;iy<Ly;iy++){
      //Calcule los campos macroscopicos en la celda
      rho0 = rho(ix,iy,false);
      for(i=0;i<Q;i++){     //para cada vector de velocidad
        n0 = n(ix,iy,i);
        fnew[n0] =UmUtau*f[n0]+Utau*feq(rho0,Ux0,Uy0,i)+(3.0*Snew[ix][iy][i])/2.0-Sold[ix][iy][i]/2.0; // con forzamiento
      }
    }  
}

//Imponer campos

void LatticeBoltzmann::ImposeFields(void){
  int ix,iy,i,n0;
  for(ix=0,iy=0;iy<Ly;iy++) //Pared izquierda
    for(i=0;i<Q;i++){
      n0=n(ix,iy,i);
      fnew[n0]=feq(1.0,0,0,i);
    }
  for(ix=Lx-1,iy=0;iy<Ly;iy++) //Pared derecha
    for(i=0;i<Q;i++){
      n0=n(ix,iy,i);
      fnew[n0]=feq(1.0,0,0,i);
    }
  for(ix=0,iy=Ly-1;ix<Lx;ix++) //Pared arriba
    for(i=0;i<Q;i++){
      n0=n(ix,iy,i);
      fnew[n0]=feq(1.0,0,0,i);
    }
  for(ix=0,iy=0;ix<Lx;ix++) //Pared abajo
    for(i=0;i<Q;i++){
      n0=n(ix,iy,i);
      fnew[n0]=feq(1.0,0,0,i);
    }
}


//Adveccion
void LatticeBoltzmann::Advection(double Ux0, double Uy0, int t){
  int ix, iy, i, ixnext, iynext, ixback, iyback, n0, n0next;
  for(ix=0;ix<Lx;ix++)      //para cada celda
    for(iy=0;iy<Ly;iy++)
      for(i=0;i<Q;i++){     //en cada direccion
	if( (ix==0) || (ix==Lx-1) || (iy==0) || (iy==Ly-1) ){ continue; }
        ixnext = (ix+Vx[i]+Lx)%Lx; iynext = (iy+Vy[i]+Ly)%Ly;
        //un if que corte las paredes
        n0 = n(ix,iy,i); n0next = n(ixnext,iynext,i);
        f[n0next] = fnew[n0];     //fronteras periodicas
        //Sold = Snew
	Sold[ix][iy][i] = Snew[ix][iy][i];
        Snew[ix][iy][i] = Si(ix,iy,Ux0,Uy0,t+1,i);
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

double LatticeBoltzmann::Varianza(void){
  int ix, iy;
  double N = 0;
  double xprom = 0;
  double yprom = 0;
  double sum = 0;
  double var = 0;

  
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      N+=rho(ix,iy,true);
    }
  }
  
  N=N-(Lx*Ly);
  
  for(ix=0;ix<Lx;ix++){
    for(iy=0;iy<Ly;iy++){
      yprom+=iy*(rho(ix,iy,true)-1.0); 
    }
  }

  yprom/=N;
  
  for(iy=0;iy<Ly;iy++){
    for(ix=0;ix<Lx;ix++){
      xprom+=ix*(rho(ix,iy,true)-1.0); 
    }
  }
  
  xprom/=N;

for(ix=0; ix<Lx; ix++){
    for(iy=0; iy<Ly; iy++){
      sum += (pow(iy-yprom,2) + pow(ix-xprom,2.0))*(rho(ix,iy,true)-1.0);
    }
  }
 
 var = sum/(N-1);

  return var;
}

//-----Programa Principal-----

int main(void){
  LatticeBoltzmann Ondas;
  int t, tmax = 50;
  double rho0 = 1.0, Ux0 = 0.3, Uy0 = 0.0; 

  Ondas.Start(rho0, Ux0, Uy0);
  //Evolucione
  for(t=1; t<tmax; t++){
    //if(t>3){
    //cout << t << '\t' << Ondas.Varianza() << endl;
    //}
    Ondas.Collision(Ux0,Uy0);
    Ondas.ImposeFields();
    Ondas.Advection(Ux0, Uy0, t);
    
  }
  //Imprima
  Ondas.Print("datos.dat");

  return 0;
}
