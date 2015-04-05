#include <iostream>
#include <math.h>
#include <fstream>
#include <string>
//#include <vector>
#include <sstream>

#include "ivps.h"

using namespace std;

void f ( double t, double *x, double *xp ) { 
	xp[0] = -20.0*x[0] + 10.0*x[1] ;
	xp[1] = 10.0*x[0]- 20.0*x[1] + 10.0*x[2];
	xp[2] =            10.0*x[1] - 20.0*x[2] + 10.0*x[3];
	xp[3] =                        10.0*x[2] - 20.0*x[3] + 10.0*x[4];
	xp[4] =                                    10.0*x[3] - 20.0*x[4];
} 

double exactSol(int a,double t)
{
	switch (a)
	{
	case 0:
		return exp(-30*t)*exp(-10*t*(pow(3.0 , 0.5) + 2))*((exp(10*t)*exp(10*t*(pow(3.0 , 0.5) + 2)))/3 + (pow(3.0 , 0.5)*exp(30*t)*(2*pow(3.0 , 0.5) - 3))/18 + exp(30*t)*exp(10*t*(pow(3.0 , 0.5) - 2))*exp(10*t*(pow(3.0 , 0.5) + 2))*(pow(3.0 , 0.5)/6 + 1/3))
		;
		break;
	case 1:
		return -exp(-30*t)*exp(-10*t*(pow(3.0 , 0.5) + 2))*((exp(30*t)*(2*pow(3.0 , 0.5) - 3))/6 - pow(3.0 , 0.5)*exp(30*t)*exp(10*t*(pow(3.0 , 0.5) - 2))*exp(10*t*(pow(3.0 , 0.5) + 2))*(pow(3.0 , 0.5)/6 + 1/3))
			;//-exp(-30*t)*exp(-10*t*(3^(1/2) + 2))*((exp(30*t)*(2*3^(1/2) - 3))/6 - 3^(1/2)*exp(30*t)*exp(10*t*(3^(1/2) - 2))*exp(10*t*(3^(1/2) + 2))*(3^(1/2)/6 + 1/3))

		break;
	case 2:
		return exp(-20*t)*exp(-10*t*(pow(3.0 , 0.5) + 2))*((pow(3.0 , 0.5)*exp(20*t)*(2*pow(3.0 , 0.5) - 3))/9 - exp(10*t*(pow(3.0 , 0.5) + 2))/3 + 2*exp(20*t)*exp(10*t*(pow(3.0 , 0.5) - 2))*exp(10*t*(pow(3.0 , 0.5) + 2))*(pow(3.0 , 0.5)/6 + 1/3))
;
		break;
	case 3:
		return -exp(-30*t)*exp(-10*t*(pow(3.0 , 0.5) + 2))*((exp(30*t)*(2*pow(3.0 , 0.5) - 3))/6 - pow(3.0 , 0.5)*exp(30*t)*exp(10*t*(pow(3.0 , 0.5) - 2))*exp(10*t*(pow(3.0 , 0.5) + 2))*(pow(3.0 , 0.5)/6 + 1/3))
;
		break;
	case 4:
		return exp(-30*t)*exp(-10*t*(pow(3.0 , 0.5) + 2))*((exp(10*t)*exp(10*t*(pow(3.0 , 0.5) + 2)))/3 + (pow(3.0 , 0.5)*exp(30*t)*(2*pow(3.0 , 0.5) - 3))/18 + exp(30*t)*exp(10*t*(pow(3.0 , 0.5) - 2))*exp(10*t*(pow(3.0 , 0.5) + 2))*(pow(3.0 , 0.5)/6 + 1/3))
;		
		break;

	default:
		cout<<"error in switching"<<endl;
		return 1000;
	}
}

double norm2(double* x0 ,double t0 )
{
	double m_sum = 0.0;
	for (int i=0; i<5;++i)
		m_sum  +=pow( x0[i] , 2);
		
	return sqrt(m_sum);
}


int main ()

{
	double t0, x0[5], h, work[9];
	int n, j, ts, call_num;

	ofstream myfile; //
	myfile.open ("0.0270.txt");//plot
	h = 0.0270;//time step
	int endTime=10;
	n = (int)(endTime/h); //numbers of step -1
	for ( j = 0; j < 1; j++ ) {
		x0[0] = 1.0;  x0[1] = 1.0;  x0[2] = 1.0; x0[3] = 1.0;  x0[4] = 1.0;
		t0 = 0.0;
		//h = 10.0 / ( (double) n ); for integer time step;
		call_num = 1;
		myfile<<t0<<" "<<norm2(x0,t0)<<endl;

		//myfile<<t0<<" "<<exactSol(1,t0)<<" "<<x0[1]<<endl;
		//myfile<<t0<<" "<<x0[0]<<" "<<x0[1]<<" "<<x0[2]<<" "<<x0[3]
		//<<" "<<x0[4]<<endl;

		for ( ts = 1; ts <= n; ts++ ) {
			ab2 ( 5, t0, x0, h, &call_num, work, f );
			t0 += h;
			//myfile<<t0<<" "<<exactSol(1,t0)<<" "<<x0[1]<<endl;
			myfile<<t0<<" "<<norm2(x0,t0)<<endl;
			//myfile<<t0<<" "<<x0[0]<<" "<<x0[1]<<" "<<x0[2]<<" "<<x0[3]
			//<<" "<<x0[4]<<endl;			
		}

		t0=endTime; //last time step marching
		ab2 ( 5, t0, x0, h, &call_num, work, f );
		//myfile<<t0<<" "<<exactSol(1,t0)<<" "<<x0[1]<<endl;
		myfile<<t0<<" "<<norm2(x0,t0)<<endl;		
		//myfile<<t0<<" "<<x0[0]<<" "<<x0[1]<<" "<<x0[2]<<" "<<x0[3]
		//<<" "<<x0[4]<<endl;



		// ofstream myfile;
		//myfile.open ("error.txt");
	/*	cout << x0[0] << "\t" << fabs(x0[0]-exactx) << endl
			<< x0[1] << "\t" << fabs(x0[1]-exacty) << endl
			<< x0[2] << "\t" << fabs(x0[2]-exactz) << endl
			<< x0[3] << "\t" << fabs(x0[3]-exactp) << endl
			<< x0[4] << "\t" << fabs(x0[4]-exactq) << endl;
			*/
		n *= 2;
	}

}
