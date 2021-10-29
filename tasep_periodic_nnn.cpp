/***********************************************************
* Copyright (C) 2021  
* Authors: Cade Spaulding & Hamid Teimouri
* Rice university--Department of Chemistry
* This file is distributed under the terms of the
* GNU General Public License as published by the
* Free Software Foundation; either version 3 of the
* License, or (at your option) any later version.
* http://www.gnu.org/copyleft/gpl.txt
***********************************************************/
//==============================================================================================================//
// Monte Carlo Siumulation of totally asymmetric simple exclusion process (TASEP) for periodic boundaries (ring)
// with next-nearest neighbor interacting
//==============================================================================================================//
#include <iostream>
#include <fstream>
#include <sstream>
#include <math.h>
#include <time.h>
#include <numeric>
#include <cstdlib>
#include <vector>
#include <valarray>
#include <algorithm>
#include <cstddef>
#include <iomanip>
#include <ctime>
#include "ran3.h"
#include "cpu_time.h"
#include <cmath>
#pragma hdrstop
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)
#define L 200 // number of lattice sites
#define T 10000 // number of MC steps
using std::vector;
using namespace std;
long int dum;
const int Teq=T/double(5);
double  rho=0.8;
double N=0.0;
double  J=0.0;
int     site=0;
int     nextsite=0;
long int j, num, k, i;
double  E=-2;
double q = exp(E/2);
double r = exp(-E/2); 
double rate;
double deltar=0.05;
double  dt=0.0;
double  t=0.0;
long double SiteDens;
long double AvgDens;
double densprof[L+1], lattice[L+1], possibleEnter[L+1];
//=====================================================================================================
//========================== initialize the lattice ==================================================
 void possible_enter()
{
	i = 1;
	k = 1;
	while (i < L)
	{
		if (lattice[i] == 0)
		{
			possibleEnter[k] = i;
			i++;
			k++;
		}
		else
		{
			i++;
		}	
	}
}
ofstream q1("Flux_nnn.txt");


//=====================================================================================================
//=========================== MAIN CODE ===============================================================
//=====================================================================================================
int main()
{

  dum=-time(NULL);
  ran3(&dum);

  double ctime;
  double ctime1;
  double ctime2;


 cout<< "Enter the value of particle denisty:" ;
	cout<< "\nrho = ";
	cin>> rho;
	cin.ignore();
	cout<< "      "<<endl;

 N=rho*L; 

 cout<< "Enter the value of interaction energy:";
	cout<< "\nE = ";
	cin>> E;
	cin.ignore();
cout<< "      "<<endl;

 q = exp(E/4);
 r = exp(-E/4); 

 dt=1/(1 + q*q + r*r );
//=====================================================================================================
// STEP 1 -- initialize varables
	t=0.0; // Initial time.
	J=0.0; // Initial current
	AvgDens = 0.0;
	num = 0;
	for (j = 1; j <= L; j++)
	{
		lattice[j]=0; densprof[j]=0; 
	}
	while (num < N)
	{	
		possible_enter();
		site = possibleEnter[rand()%k + 1];
		lattice[site] = 1;
		num++;
	}
//=====================================================================================================

  //update loop

for (t = 0.0; t <= T; t += dt)
	{

	for (j = 1;j <= L; j++) // Goes through as many iterations as the array's length.
	{
		site = rand()%L + 1; // Picks out a random site in the array.
		int neiborhood[6]={(site-2),(site-1),(site),(site+1),(site+2),(site+3)};
        //int S;
		for (int S = 0; S <= 4; S += 1){
			if (neiborhood[S]>L){
				neiborhood[S]-=L;
			}
			if (neiborhood[S]<1){
				neiborhood[S]+=L;
			}
        
        }
        int nextsite = neiborhood[3]; // Defines the next site. Saves the program from unecessary calculations.
        int nextnextsite=neiborhood[4];
		int nextnextnextsite=neiborhood[5];
        int prevsite=neiborhood[1];
		int prevprevsite=neiborhood[0];

        double rate=1;
        if (lattice[nextnextsite]){
        rate=rate*q;
        }
        if (lattice[prevsite]){
        rate=rate*r;
        }
		if (lattice[nextnextnextsite]){
        rate=rate*q;
        }
        if (lattice[prevprevsite]){
        rate=rate*r;
        }
		if (dt*rate>1){
			cout<<"stop"<<endl;

		}
		if (lattice[site] == 1 && lattice[nextsite] == 0 && ran3(&dum) <= dt*rate) // Bulk of the array.
		{
			lattice[site] = 0; // If the criteria are met the particle leave site.
			lattice[nextsite] = 1; // And enters the next site.

			if(site == int(L/2) && t >= Teq) // Measuring current. It's measured when a particle crosses the midpoint of the array.
			{
				J++;
			}	
		}		
	
	} // End of loop through array. J-loop.

		if (t >= Teq) // Building the density profile after the system has reached steady state. Not averaged yet.
		{
			for (j = 1; j <= L; j++) // Sweeps throught the array.
			{
				if (lattice[j] == 1) // If it finds a site which contains a particle.
				{
					densprof[j] += dt; // It adds one to its density counter. In reality this can be >1, it will be time-averaged later.
				}
			}
		}

	} // End of time loop.

 q1<<rho<<" "<<E<<" "<<J/double((T-Teq))<<endl;


 cout<<"Particle flux:"<<endl;

 cout<<"J="<<" "<<J/double((T-Teq))<<endl;


  ctime2 = cpu_time ();
  ctime = ctime2 - ctime1;
  cout << "\n";
  cout << "  Elapsed cpu time for main computation:\n";
  cout << "  " << ctime2 << " seconds.\n";

  q1.close();
  //q2.close();
  //q3.close();
  cout<<" end\n "<<endl;
  cin.get();
  
  return 0;
}
 
