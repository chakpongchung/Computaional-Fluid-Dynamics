//2-stage 3rd order implicit Runge-Kutta method
//Chak-Pong CHUNG  
//email : chakpongchung@gmail.com

#include <iostream>
#include <fstream>
using namespace std;


long double xi1(double h, double y_n)
{
	return (6-4*h) * y_n / (h*h -4*h +6);
}

long double xi2(double h, double y_n)
{
	return 6 * y_n / (h*h -4*h +6);
}

int main()
{
	long double n=20;
	long double h=1/n;
	long double a11=1.0/4, a12= -1.0/4, a21=1.0/4,a22=5.0/12;
	double b1=1.0/4, b2=3.0/4;

	
	long double y_n=1;
	ofstream myfile;
	string outFile="output.dat";

	myfile.open (outFile);
	
	for (int i= 0; i<n+1;i++)
	{
		
		myfile<<i*h<<" "<<y_n<<endl;
		cout<<i<<"\t"<<xi1(h,y_n)<<"\t"<<xi2(h,y_n)<<"\t"<<y_n<<"\t"<<exp(i*h)<<"\t"<<abs(exp(i*h)-y_n)<<endl;

		y_n= y_n+ h*b1*xi1(h,y_n)+ h*b2*xi2(h,y_n);

		
	}
	


}