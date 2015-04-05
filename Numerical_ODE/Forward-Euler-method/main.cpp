//Forward Euler method
//Chak-Pong CHUNG  
//email : chakpongchung@gmail.com

/* Euler for a set of first order differential equations */

#include <stdio.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>

using namespace std;
double dist=1.0/20;
//#define dist=0.25              /* stepsize in t*/
#define MAX 1.0                /* max for t */
const int N=1;                // N equation(s) solver

double  f(double x, double y[], int i)
	//double f(double x , vector<double> y, int i)
{
	//return y[0]*(1-y[0]);                /* derivative of first equation */
	return y[0]; 

} 

void euler(double t, double y[], double step) /* Euler function */
	//void euler(double x , vector<double> y, double step)
{
	double  s[N];      /* for euler */

	int i;
	{
		for (i=0;i<N;i++)
		{     s[i]=step*f(t, y, i);
		}
	}

	{
		for (i=0;i<N;i++) 
			y[i]+=s[i];
	}
}


double analyticSol(double t)
{
	//return 1/(9+exp(t)) * exp(t);
	return exp(t);
}

int main()
{


	double t, y[N]={};
	//	double t;
	//vector<double> y;
	y[0]=1;                                       /* initial condition */

	ofstream myfile;
	double error=0;

	myfile.open (to_string(dist));
	for (int j=0; j*dist<=MAX ;j++)                     /* time loop */
	{
		t=j*dist;

		error=abs(y[0]-analyticSol(t));

		printf("t = %f ,e(t) = %f\t ,log(dist)=%f\t ,log(error) = %f\n", 
			t,error,log(dist),log(error));

		myfile<<t<<" "<<y[0]<<endl;

		euler(t, y, dist);

	}
	myfile.close();
}





