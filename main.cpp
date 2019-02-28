// Laplasse.cpp : Defines the entry point for the console application.
// Эта программа будет считывать файл матрицу

#include "stdafx.h"
#include "LinAlg.h"
#include "laplasse.h"


using namespace std;
using namespace alglib;


int main()
{
	int ne, np;

	cout << "Enter number of T2 values (1, 2, 3): ";
	cin >> np;

	cout << endl << endl << "Enter number of Echoes (1000 - 2000): ";
	cin >> ne;
	
	Echo_Train SIM1(np,ne);
	SIM1.get_init_vars();
	SIM1.get_fid();
	SIM1.save("echo_train.dat");

	system("pause");

	return 0;
}

