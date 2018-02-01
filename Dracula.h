/*  
    DRACULA: Decay, Recombination And Collisions compUtation for Low-temperature Antihydrogen
    Copyright (C) 2013  B. Radics, D. Murtagh, Y. Yamazaki, F. Robicheaux

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    Contact: balint DOT radics AT cern.ch
*/



#ifndef _DRACULA_H_
#define _DRACULA_H_

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>

using namespace std;


//------------------------------------------------
// DRACULA: Decay, Recombination And Collisions compUtation for Low-temperature Antihydrogen 
//
// Code written by B. Radics and F. Robicheaux for simulation of Antihydrogen
// formation
//
// Dec. 2013 
// -----------------------------------------------
//
// The ODE system is defined in terms of GSL in the following way:
// - n-dimensional first-order system: dy_{i}(t)/dt = f(t, y_{1}(t), ..., y_{n}(t));
// - y is a vector of n-dimension, i=1..n;
// - the GSL stepping functions rely on the vector of derivatives, f_{i}, and the Jacobian matrix;
//  J_{ij} = df_{i}(t,y(t))/dy_{j}
// - Because the vector is indexed by the principal quantum number, n, the derivatives are non-zero only
//		   if there is a transition between the states i<->j;
  


//------------------------------------------------

const double e0 = 13.606;// eV
const double h = 6.626e-34;//m^{2}kg/s
const double me = 9.109e-31;//kg
const double mp = 1.672e-27;//kg
const double kb = 1.38e-23; // J/K
const double eV = 1.6e-19; // J
const double pi = 3.141592654;
const double eps0 = 8.854e-12; // F/m

// These guys here below must be global and static so that 
// GSL (Gnu Scientific Library) can find them (for the time being)

// Number of all states
// index 0 is for the free/continuum positrons
static int m_Nst;
// Number of active quantum states (princ. qn only), Nactst < Nst
static int m_Nactst; 

// Nst x Nst matrix for collisional excitation from ni to nf
static double ** m_kcol = NULL;
// Nst x Nst matrix bookkeeping for detailed balance check
static bool ** m_kcol_db = NULL;
// Nst x Nst matrix for radiative decay rates
static double ** m_gamrad = NULL;
// ionisation rate array, dim = Nst
static double * m_kion = NULL;
// TBR rate array, dim = Nst
static double * m_ktbr = NULL;
// RR rate array, dim = Nst
static double * m_krr = NULL;


class Dracula{

 public:
  Dracula(int Nstates);
  ~Dracula();

  // [K]
  void SetPositronTemp(double temp){m_Temp = temp;};
  // [m^{-3}]
  void SetPositronDensity(double dens){m_Density = dens;};
  // []
  void SetAntiprotonN(double n){m_NumbP = n;};
  // [us]
  void SetOverlapTime(double time){m_Time = time;};
  // [T]
  void SetBField(double bfield){m_Bfield = bfield;};

  // Use Pohl rates
  void SetUsePohlRates(bool useit){m_UsePohlRates = useit;};

  // Pick scattering rate data from database file
  // This picks up the scattering rates already
  void SetSQLDbScat(char * dbfile, bool isOverlap = true);
  void SetSQLDbRadDec(char * dbfile, bool isOverlap = true);

  // the external file containing 3 columns: ni, nf and rate(-coefficient)
  // in units of [m^{3}/s]
  void SetScatteringTable(char * scatfile);
  // the external file containing 2 columns: n, rate
  // in some units...
  void SetRadRecTable(char * rrfile);
  // the external file containing the radiative decay rates from ni to nf, 
  // 3 columns: ni, nf and rate in [1/s] units
  // or (at the moment) 5 columns: ni, li, nf, lf, rate
  void SetRadDecTable(char * rdfile);

  // For T=0 we use the Pohl scattering rates
  void UsePohlRates();

  // steering black-body radiation process for stimulated absorption and emission
  void SetBlackBodyTemp(double temp){m_UseBlackBody = true;m_TempRad = temp;};
  void AddBlackBody();
  void ApplyDetailedBalance();
  void ApplyEstabPop();

  double * GetPop(){return m_LevPop;};
  void SetPop(double * popv);

  bool Compute();// return true if success

  void PrintLevelPop(); // to screen
  void SaveLevelPop(const char * outputfile); // to ascii file
  void SaveLevelPopCSV(const char * outputfile);// to CSV file
  // for parameterscan
  void SavePop(int i, const char * outputfile);// append parameters and population in quantum state i

  // Thermal population from Saha-Boltzman distribution
  double GetThermPop(int n);

  static char * ftoa(float x);
  static char * Etoa(double x);
  static char * itoa(int x);

  // For GSL, the diff. equation
  static int func (double t, const double y[], double f[],
	  void *params);
  // For GSL the Jacobi matrix
  static int jac (double t, const double y[], double *dfdy, 
	   double dfdt[], void *params);



 public:
  // Temperature of positrons
  double m_Temp;
  // Temperature of radiation field for Black-Body
  double m_TempRad;
  // Density of positrons
  double m_Density;
  // Density of antiprotons
  double m_NumbP;
  // Antiproton-positron Overlap time
  double m_Time;
  // Magnetic field strength
  double m_Bfield; 

  bool m_UsePohlRates;//B=0T

  // whether to use black body spectrum for stimulated emission/absorption
  bool m_UseBlackBody;

  // switches
  // whether to have collisional (de)excitation and/or ionisation, tbr, rad. rec.
  bool m_isOverLap;
  // Use user level populvation
  bool m_useLevPop;


  // the level population distribution array, dim = Nst
  double * m_LevPop;
  // user level population for initial condition
  double * m_UserLevPop;
  // initial population distribution array, dim = Nst
  double * m_popinit;

};

static inline void loadbar(unsigned long long x, unsigned long long n, unsigned int w);

#endif
