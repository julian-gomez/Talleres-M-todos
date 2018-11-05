#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
using namespace std;
double c = 300;
double deltaX = 0.01;
double A0 = 0.1;
double L = 1.0;
int N = int(L/deltaX);
double funCos(double t);
int main()
{
	ofstream archivo_salida("datosExtremosFijos.dat");
	int i,j,k,l,m,n,contador;
	double deltaT = 0.5*deltaX/c;
	double arregloPasado[N];
	double arregloPresente[N];
	double arregloFuturo[N];
	double intervaloTiempo = 0.1;
	int M = int(intervaloTiempo/deltaT);
	double r = c*deltaT/deltaX;
	//Condiciones iniciales
	arregloPasado[0] = 0.0;
	arregloPasado[N-1] = 0.0;
	arregloPresente[0] = 0.0;
	arregloPresente[N-1] = 0.0;
	arregloFuturo[0] = 0.0;
	arregloFuturo[N-1] = 0.0;
	for(i = 0; i < N; i++)
	{	
		if(i < int(N/2.0))
		{
			arregloPasado[i] = 2.0*A0/L*deltaX*i;
		}
		else
		{
			arregloPasado[i] = -2.0*A0/L*deltaX*i + 2*A0;
		}
		archivo_salida << arregloPasado[i] << endl;
	}
	
	
	for(j = 1; j < N-1; j ++)
	{
		arregloPresente[j] = pow(r,2)/2.0*(arregloPasado[j+1] + arregloPasado[j-1] - 2.0*arregloPasado[j]) + arregloPasado[j];
	}

	//Ahora haremos para todo t
	
	contador = 0;

	for(k = 2; k < M; k++)
	{
		for(l = 1; l < N-1; l++)
		{
			arregloFuturo[l] = pow(r,2)*(arregloPresente[l+1] + arregloPresente[l-1] - 2.0*arregloPresente[l]) + 2.0*arregloPresente[l] - arregloPasado[l];
		}
		for(m = 0; m < N; m++)
		{
			arregloPasado[m] = arregloPresente[m];
			arregloPresente[m] = arregloFuturo[m];
		}
		contador ++;
		if(contador == 300)
		{
			for(n = 0; n < N; n++)
			{
				archivo_salida << arregloPresente[n] << endl;
			}
			contador = 0;
		}
	
	}

	//Ahora, consideramos un extremo libre
	//Establecemos de nuevo las condiciones iniciales
	
	arregloPasado[0] = 0.0;
	arregloPasado[N-1] = 0.0;
	arregloPresente[0] = 0.0;
	arregloFuturo[0] = 0.0;
	ofstream archivo_salida2("datosExtremoLibre.dat");

	for(i = 0; i < N; i++)
	{	
		if(i < int(N/2.0))
		{
			arregloPasado[i] = 2.0*A0/L*deltaX*i;
		}
		else
		{
			arregloPasado[i] = -2.0*A0/L*deltaX*i + 2*A0;
		}
		archivo_salida2 << arregloPasado[i] << endl;
	}
	for(j = 1; j < N-1; j ++)
	{
		arregloPresente[j] = pow(r,2)/2.0*(arregloPasado[j+1] + arregloPasado[j-1] - 2.0*arregloPasado[j]) + arregloPasado[j];
	}
	arregloPresente[N-1] = arregloPresente[N-2];

	//Ahora haremos para todo t
	
	contador = 0;

	for(k = 2; k < M; k++)
	{
		for(l = 1; l < N-1; l++)
		{
			arregloFuturo[l] = pow(r,2)*(arregloPresente[l+1] + arregloPresente[l-1] - 2.0*arregloPresente[l]) + 2.0*arregloPresente[l] - arregloPasado[l];
		}
		arregloFuturo[N-1] = arregloFuturo[N-2];
		for(m = 0; m < N-1; m++)
		{
			arregloPasado[m] = arregloPresente[m];
			
		}
		arregloPasado[N-1] = arregloPasado[N-2];
		for(int h = 0; h < N-1; h++)
		{
			arregloPresente[h] = arregloFuturo[h];
		}
		arregloPresente[N-1] = arregloPresente[N-2];
		
		contador ++;
		if(contador == 300)
		{
			for(n = 0; n < N; n++)
			{
				archivo_salida2 << arregloPresente[n] << endl;
			}
			contador = 0;
		}
	
	}
	
	//Ahora consideremos el caso 3: un extremos fijo y el otro bajo un movimiento cosenosoidal. Consideraremos el extremo fijo al primer punto

	arregloPasado[0] = 0.0;
	arregloPasado[N-1] = 0.0;
	arregloPresente[0] = 0.0;
	arregloFuturo[0] = 0.0;
	ofstream archivo_salida3("datosExtremoCoseno.dat");

	for(i = 0; i < N-1; i++)
	{	
		arregloPasado[i] = 0.0;
		archivo_salida3 << arregloPasado[i] << endl;
	}
	arregloPasado[N-1] = funCos(0.0);
	archivo_salida3 << arregloPasado[N-1] << endl;
	
	
	for(j = 1; j < N-1; j ++)
	{
		arregloPresente[j] = pow(r,2)/2.0*(arregloPasado[j+1] + arregloPasado[j-1] - 2.0*arregloPasado[j]) + arregloPasado[j];
	}
	arregloPresente[N-1] = funCos(deltaT);

	//Ahora haremos para todo t
	
	contador = 0;

	for(k = 2; k < M; k++)
	{
		for(l = 1; l < N-1; l++)
		{
			arregloFuturo[l] = pow(r,2)*(arregloPresente[l+1] + arregloPresente[l-1] - 2.0*arregloPresente[l]) + 2.0*arregloPresente[l] - arregloPasado[l];
		}
		arregloFuturo[N-1] = funCos(k*deltaT);
		for(m = 0; m < N; m++)
		{
			arregloPasado[m] = arregloPresente[m];
		}
		for(int h = 0; h < N; h++)
		{
			arregloPresente[h] = arregloFuturo[h];
		}	
		contador ++;
		if(contador == 300)
		{
			for(n = 0; n < N; n++)
			{
				archivo_salida3 << arregloPresente[n] << endl;
			}
			contador = 0;
		}
	
	}

	return 0;
}

double funCos(double t)
{
	return A0*cos(3.0*c*t*M_PI/L);
}










