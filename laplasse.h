#pragma once

using namespace alglib;



class Echo_Train
{
	real_1d_array fid;
	real_1d_array tau;
	real_1d_array porosity;					// The total porosity should be between 0pu and 60pu. Amplitude Summ P1+P2+P3 = 1 ; 
	real_1d_array T2_init;					// NMR properties of the formation defined by up to three T2 values
	int number_of_poros;					// Inintial poros number 1, 2 or 3
	double wait_time;						// The wait time could be  2, 6, or 12 seconds
	double noise_deviation;					//Gaussian noise with a standard deviation between 1pu and 4pu
	double ratio = 1.5;						// T1 = ratio * T2 for initial ratio = 1.5
	double inter_echo_tau = 0.6;
public:
	int number_of_echoes;					// The number of echoes could be 1000 or 2000
	Echo_Train(int number_poros, int number_of_echoes);
	void get_fid(int noise_deviation);
	void get_fid();
	void put_echo(int);
	void get_init_vars();
	bool save(const char* file_name);
};

class Inversion_T2
{
	real_1d_array Amplitude;
	real_1d_array T2;
	real_1d_array F;		// Experimental Echo Train
	real_1d_array tau;
	real_1d_array bnd_min;  // Minimal Pu = 0 pu
	real_1d_array bnd_max;  // Maximal Pu = 60 pu 
	double epsx = 0.00001;    // Err function criteriom precision
	ae_int_t maxits = 0;
	minlmstate state;
	minlmreport rep;

public:
	Inversion_T2();
	void T2_spread(double first_bin, double last_bin, double bin_increment);
	void inversion();
	void function1_fvec(const alglib::real_1d_array &x, alglib::real_1d_array &fi, void *ptr);
	bool save(const char* file_name);
};


//void minlmoptimize(minlmstate &state,
//	void(*fvec)(const real_1d_array &x, real_1d_array &fi, void *ptr),
//	void(*rep)(const real_1d_array &x, double func, void *ptr) = NULL,
//	void *ptr = NULL,
//	const xparams _xparams = alglib::xdefault);

