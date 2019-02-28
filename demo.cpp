// demo.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"

#include <stdio.h>
#include <windows.h>
#include "LinAlg.h"
#include <stdlib.h>
#include <math.h>
#include <random>
#include "optimization.h"
#include <fstream>

using namespace std;
using namespace alglib;


// Globals


double first_bin = 0.5;
double last_bin = 4096;
const int j_max = 53;

double T2[j_max];

int number_of_poros;					// Poros number 1, 2 or 3
double porosity[3];				        // The total porosity should be between 0pu and 60pu. Amplitude Summ P1+P2+P3 = 1 ; 
double T2_init[3];						// NMR properties of the formation defined by up to three T2 values
double T1[3] = { 2000, 2000, 2000 };	// 
										// 
int number_of_echoes = 2000;			// The number of echoes could be 1000 or 2000
double tau[2000] = { 0 };				// tau - linear scale
double F[2000] = { 0 };					// F(t) simulation 
double ratio = 1.5;						// T1 = ratio * T2  

double wait_time;						// The wait time could be  2, 6, or 12 seconds
const double inter_echo_tau = 0.6;

double noise_deviation = 3;				//Gaussian noise with a standard deviation between 1pu and 4pu
double x_arr[j_max];					// 
double bnd_min_arr[j_max];				// Minimal Pu = 0 pu
double bnd_max_arr[j_max];				// Maximal Pu = 0 pu


void fid_simulation(int noise)
{
	// HERE fid simulation
	// Noise option: 0 - off, 1,2,3,4 - on
	// Gaussian noise with a standard deviation between 1pu and 4pu  

	default_random_engine generator;
	normal_distribution<double> distribution(0.0, 1.0);

	for (int current_echo = 0; current_echo < number_of_echoes; current_echo++)
	{
		tau[current_echo] = current_echo*inter_echo_tau;

		for (int current_pore = 0; current_pore < number_of_poros; current_pore++)
		{
			F[current_echo] =
				F[current_echo] + porosity[current_pore] * (1 - exp(-wait_time / (ratio*T2_init[current_pore])))
				* exp(-tau[current_echo] / T2_init[current_pore])
				+ noise * distribution(generator);
		}
	}
}

bool save_simulation()
{
	ofstream fout;
	fout.open("echoes.dat");

	for (int i = 0; i < number_of_echoes; i++)
	{
		fout << tau[i] << ", " << F[i] << endl;
	}
	fout.close();
	return FALSE;
}

bool save_distribution_T2(real_1d_array &x)
{
	ofstream fout;
	fout.open("distr.dat");

	for (int j = 0; j < j_max; j++)
	{
		fout << T2[j] << ", " << x(j) << endl;
	}
	fout.close();
	return 0;
}

void init_variables()
{
	cout << "Enter number of T2 values (1, 2, 3): ";
	cin >> number_of_poros;

	for (int i = 0; i < number_of_poros; i++)
	{
		cout << endl << endl << "Enter value of T2[" << i + 1 << "] in msec from 1 to 3000 msec: ";
		cin >> T2_init[i];
	}

	for (int i = 0; i < number_of_poros; i++)
	{
		cout << endl << endl << "Enter value of  Porosity [" << i + 1 << "] in pu from 0 to 50pu: ";
		cin >> porosity[i];
	}

	cout << endl << endl << "Enter the Wait Time in sec (2 - 12): ";
	cin >> wait_time;
	wait_time *= 1000;

	cout << endl << endl << "Enter number of Echoes (1000 - 2000): ";
	cin >> number_of_echoes;

	cout << endl << endl << "The inter-echo time is 0.6 ms";

	cout << endl << endl << "Enter the noise deviation between 0pu and 4pu: ";
	cin >> noise_deviation;
}

void T2_spread(double first_bin, double last_bin, double bin_increment)
{
	// HERE SPREADING
	//
	//	The T2 relaxation times are equally spread on a logarithmic scale
	//	The first bin is 0.5ms
	//	The last bin is  4096ms
	//	The bin increment is sqrt(2)

	for (int j = 0; j < j_max; j++)
	{
		T2[j] = first_bin*pow(2, j*bin_increment);
	}
}

void  function1_fvec(const real_1d_array &x, real_1d_array &fi, void *ptr)
{
	//	HERE callback calculates of err func
	//	Err(x1,x2,...xi) = summa(i) { f1^2 + f2^2 + ... fi^2}
	//  where   fi(xi) = F(i) - summa(j) { A(j)*exp[ -t(i)/T2(j) ] }

	for (int i = 0; i < number_of_echoes; i++)
	{
		fi[i] = 0;
		for (int j = 0; j < j_max; j++)
		{
			fi[i] = fi[i] + x[j] * exp(-tau[i] / T2[j]);
		}
		fi[i] = F[i] - fi[i];
	}
}


void inversion()
{
	// From Here inversion Start
	// Improved Levenberg-Marquardt optimizer
	// Bound constrained nonlinear least squares optimization


	real_1d_array x;
	real_1d_array bnd_min;
	real_1d_array bnd_max;

	for (int j = 0; j < j_max; j++)
	{
		x_arr[j] = 0;			// 
		bnd_min_arr[j] = 0;		// Minimal Pu = 0 pu
		bnd_max_arr[j] = 60;	// Maximal Pu = 60 pu
	}

	x.attach_to_ptr(j_max, x_arr);
	bnd_min.attach_to_ptr(j_max, bnd_min_arr);
	bnd_max.attach_to_ptr(j_max, bnd_max_arr);

	double epsx = 0.00001;    // Err function criteriom precision
	ae_int_t maxits = 0;
	minlmstate state;
	minlmreport rep;

	minlmcreatev(number_of_echoes, x, 0.0001, state);
	minlmsetbc(state, bnd_min, bnd_max);
	minlmsetcond(state, epsx, maxits);
	alglib::minlmoptimize(state, function1_fvec); // <- Здесь и будет проблема в классах!!!!!
	minlmresults(state, x, rep);

	if (!save_distribution_T2(x))
	{
		cout << endl << "Distribution saved in file: distr.dat" << endl << endl;
	}
}


int main(int argc, char **argv)
{
	init_variables();
	fid_simulation(noise_deviation);		// Noise option : 0 - off; 1, 2, 3, 4 - on
	if (!save_simulation())
	{
		cout << endl << "Echo-train saved in file: echoes.dat" << endl;
	}
	double bin_increment = log2(last_bin / first_bin) / (j_max - 1);
	T2_spread(first_bin, last_bin, bin_increment);
	inversion();

	system("pause");
	return 0;
}
