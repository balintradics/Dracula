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



#include <sqlite3.h>
#include <time.h>
#include <iomanip>
#include "Dracula.h"


// Constructor
Dracula::Dracula(int Nstates){
  
  m_Nst = Nstates;
  m_Nactst = m_Nst;
  
  m_LevPop = new double[m_Nst];
  m_UserLevPop = new double[m_Nst];
  m_popinit = new double[m_Nst];
  m_kion = new double[m_Nst];
  m_ktbr = new double[m_Nst];
  m_krr = new double[m_Nst];
  
  m_kcol = new double*[m_Nst];
  m_kcol_db = new bool*[m_Nst]; // to keep track of detailed balance
  m_gamrad = new double*[m_Nst];
  for(int i = 0; i < m_Nst; i++){
    m_kcol[i] = new double[m_Nst];
    m_kcol_db[i] = new bool[m_Nst];
    m_gamrad[i] = new double[m_Nst];
    m_kion[i] = 0;
    m_ktbr[i] = 0;
    m_popinit[i] = 0;
    m_LevPop[i] = 0;
    for(int j = 0; j < m_Nst; j++){
      m_kcol[i][j] = 0;
      m_kcol_db[i][j] = false;
      m_gamrad[i][j] = 0;
    }
  }

  m_Temp = 200.0;
  m_Density = 1.0e+14;// initial e+ density
  m_NumbP = 1.0e+05;// initial antiproton density
  m_Time = 5.0;
  m_Bfield = 1.0;

  // by default we don't use Pohl's scattering rates
  m_UsePohlRates = false;

  // by default no Black-Body radiation
  m_UseBlackBody = false;
  m_TempRad = 1.0;

  // by default collisional processes are on
  m_isOverLap = true;

  // by default no initial level pop
  m_useLevPop = false;

}

Dracula::~Dracula(){}

void Dracula::SetSQLDbScat(char * dbfile, bool isOverlap){

  m_isOverLap = isOverlap;

  sqlite3 *db;
  char *zErrMsg = 0;
  int rc;
  rc = sqlite3_open(dbfile, &db);
  if( rc ){
    fprintf(stderr, "Can't open database: %s, exiting...\n", sqlite3_errmsg(db));
    sqlite3_close(db);
    exit(-1);
  }

  std::string query_base = "select RATE from ";
  int ni = 98; int nf = 99;


  if(m_isOverLap){
    // check if B==0, then we use Pohl rates
    if(m_Bfield == 0.0){
      UsePohlRates();
      
    }else{
      
      // Get data for scattering rates 
      cout << "Reading in scat.rates from sqlite3 file " << dbfile << "... " << endl;
      for(int i = 1; i < m_Nst; i++){
	for(int j = 0; j < m_Nst; j++){
	  ni = i; nf = j;
	  query_base = "select RATE from ";
	  query_base = query_base + "XRATES " + " where NI = " + string(itoa(ni)) +  " and NF = " + string(itoa(nf)) + " and BFIELD = " + string(ftoa(m_Bfield)) + " and TEMP = " + string(itoa(m_Temp));
	  char ** result = NULL;  int nrows, ncols;
	  rc = sqlite3_get_table(db, query_base.c_str(), &result, &nrows, &ncols, &zErrMsg);
	  //    cout << nrows << ", " << ncols << endl;
	  if((nrows == 1 && ncols == 1) && atof(result[1]) > 0){
	    //	cout << "Result: " << result[1] << endl;
	    if(nf == 0){
	      m_kion[i] = m_Density*atof(result[1]);
	      //	      cout << "ion[" << i << "]: " << m_kion[i] << endl;
	      //	      cout << "OK existing result for sql query: [ " << query_base << " ]  setting ionis rate to : " << m_kcol[i][j] << endl;
	    }
	    if(nf != 0){
	      m_kcol[i][j] = m_Density*atof(result[1]);
	      //cout << "col[" << i << "][" << j << "]: " << m_kcol[i][j] << endl;
	      //	      cout << "OK existing result for sql query: [ " << query_base << " ]  setting scatt rate to : " << m_kcol[i][j] << endl;
	    }
	  }else{

	    //	    cout << "Error: non existing result for sql query: [ " << query_base << " ]  setting rate to 0!" << endl;
	    // use Pohl rates
	    //--------------------
	    double temp = m_Temp; double dens = m_Density;double ei, ef, eg, el;
	    double echrg = 1.602e-19 ; double ryd = 13.606*echrg ; double melec = 9.109e-31 ;
	    double potfac = 8.988E9*echrg*echrg ; double kboltz = 1.38e-23 ; double pi = 2.0*asin(1.0) ;
	    double eps0 = 8.854e-12 ;
	    double rat = ryd/(kboltz*temp) ;
	    double x = 0;
	    double k0 = dens*potfac*potfac/(kboltz*temp*sqrt(melec*ryd)) ;
	    double lam = 6.626e-34/sqrt(2.0*pi*melec*kboltz*temp) ;
	    if(j == 0){
	      x = rat/(i*i) ;
	      m_kion[i] = k0*11.0*sqrt(rat)/(pow(x,2.333333)+4.38*pow(x,1.72)+1.32*x) ;
	      if(i > 10){
		x = log(m_kion[i]) - x;
		m_kion[i] = exp(x);
	      }else{
		m_kion[i] = 0.0;
	      }

	      //	      cout << "ion[" << i << "]: " << m_kion[i] << endl;
	      //	      cout << "Error: non existing result for sql query: [ " << query_base << " ]  setting ionis rate to Pohl et al rate: " << m_kion[i] << endl;
	    }
	    if(j != 0 && j != i && i > 10 && j > 10){
	      ei = rat/(i*i) ; ef = rat/(j*j) ;
	      if(ei > ef) {eg = ei ; el = ef ;}
	      else {eg = ef ; el = ei ;}
	      m_kcol[i][j] = k0*pow((ei/eg),2.5)*ef*sqrt(ef)*exp(el-ei)*
	  	(22.0/pow((eg+0.9),2.333333) + 4.5/(pow(eg,2.5)*pow((eg-el),1.33333))) ;
	      //cout << "col[" << i << "][" << j << "]: " << m_kcol[i][j] << endl;
	      //	      cout << "Error: non existing result for sql query: [ " << query_base << " ]  setting rate to Pohl et al rate: " << m_kcol[i][j] << endl;
	    }
	    if(j == i)m_kcol[j][j] = 0;
	    //--------------------
	    // m_kcol[i][j] = 0;
	    //	    if(i <= 50)m_kion[i] = 0;
	  }
	}
      }
      ApplyDetailedBalance();
          
      // Print out collisional rates
      // int Nini = 60;
      // for(int j = Nini - 10; j < Nini + 10; j++)
      // 	cout << "col[" << Nini << "][" << j << "]: " << m_kcol[Nini][j] << endl;
    
    
      // Now calculate TBR rates
      double temp = m_Temp; double dens = m_Density;double ei, ef, eg, el;
      double echrg = 1.602e-19 ; double ryd = 13.606*echrg ; double melec = 9.109e-31 ;
      double potfac = 8.988E9*echrg*echrg ; double kboltz = 1.38e-23 ; double pi = 2.0*asin(1.0) ;
      double eps0 = 8.854e-12 ;
      double rat = ryd/(kboltz*temp) ;
      double x = 0;
      double k0 = dens*potfac*potfac/(kboltz*temp*sqrt(melec*ryd)) ;
      double lam = 6.626e-34/sqrt(2.0*pi*melec*kboltz*temp) ;
      
      // Calculate TBR rate from ionisation rate using Saha-Boltzmann detailed balance
      double TTBR = 0;
      for(int j = 1 ; j < m_Nactst ; j++){
	m_ktbr[j] = m_kion[j]*j*j*lam*dens*lam*lam;
	x = log(m_ktbr[j]) + 13.606*echrg/(j*j*kboltz*temp);
	m_ktbr[j] = exp(x);
	//	cout << j << "\t" << m_ktbr[j] << endl;
	//	cout << "tbr[" << j << "]: " << m_ktbr[j] << endl;
	TTBR += m_ktbr[j];
	//     if(x < 39.0) dum1 += ktbr[j]*(1.0-exp(-kion[j]*2.e-4))/(kion[j]*2.e-4) ;
      }
      cout << "Total tbr rate: " << TTBR << " /s" << endl;
      
      // Extrapolate linearly the ionization rates to higher n
      for(int j = m_Nactst-1; j < m_Nactst; j++){
	m_kion[j] = m_kion[m_Nactst-2] + (j - (m_Nactst-2))*(m_kion[m_Nactst-2] - m_kion[m_Nactst-3])/((m_Nactst-2) - (m_Nactst-3));
	m_ktbr[j] = m_kion[j]*j*j*lam*dens*lam*lam*exp(13.606*echrg/(j*j*kboltz*temp));  
      }
      
    }
    
  }


  // Fill the radiative recombination rate array from file
  if(m_isOverLap){

    // for the moment we need to convert from atomic units rate to
    // SI
    double conv = 6.12e-15;//m^{3}/s
    cout << "Reading in rad. recombination rates  from sqlite3 file " << dbfile << "... " << endl;
    for(int i = 1; i < m_Nst; i++){
      nf = i;
      query_base = "select RATE from ";
      query_base = query_base + "RRATES " + " where NF = " + string(itoa(nf)) +  " and TEMP = " + string(itoa(m_Temp));
      char ** result = NULL;  int nrows, ncols;
      rc = sqlite3_get_table(db, query_base.c_str(), &result, &nrows, &ncols, &zErrMsg);
      //    cout << nrows << ", " << ncols << endl;
      if(nrows == 0 && ncols ==0){
	m_krr[nf] = 0.0;
      }
      if(nrows == 1 && ncols == 1){
	m_krr[nf] = m_Density*conv*atof(result[1]);
      }
      // for(int ii = 0; ii < nrows;ii++)
      //   cout << result[ii] << ", " << result[ii+1] << endl;
    }
  }


  // End reading sqlite3 database
  sqlite3_close(db);
  
  // Adding Black-Body radiation if needed
  // this depends on the radiation field temperature
  if(m_UseBlackBody == true)
    AddBlackBody();
  
    
}

void Dracula::SetSQLDbRadDec(char * dbfile, bool isOverlap){

  m_isOverLap = isOverlap;

  sqlite3 *db;
  char *zErrMsg = 0;
  int rc;
  rc = sqlite3_open(dbfile, &db);
  if( rc ){
    fprintf(stderr, "Can't open database: %s, exiting...\n", sqlite3_errmsg(db));
    sqlite3_close(db);
    exit(-1);
  }

  std::string query_base = "select RATE from ";
  int ni = 98; int nf = 99;


  cout << "Reading in decay rates  from sqlite3 file " << dbfile << "... " << endl;
  for(int i = 1; i < m_Nst; i++){
    for(int j = 1; j < m_Nst; j++){
      ni = i; nf = j;
      query_base = "select LI,RATE from ";
      query_base = query_base + "DRATES " + " where NI = " + string(itoa(ni)) +  " and NF = " + string(itoa(nf));
      char ** result = NULL;  int nrows, ncols;
      rc = sqlite3_get_table(db, query_base.c_str(), &result, &nrows, &ncols, &zErrMsg);
      //cout << nrows << ", " << ncols << endl;
      if(nrows == 0 && ncols == 0){
       	m_gamrad[ni][nf]+=0.0;
      }else{
       	for(int ii = 0; ii < nrows; ii = ii+2){
       	  m_gamrad[ni][nf]+=atof(result[ii+3])*(2*atoi(result[ii+2])+1)/(ni*ni);
      // 	  cout << "Result: ni, li, nf, rate: " << ni << ", " << result[ii] << ", " << nf << ", " << result[ii+1] << endl;
       	}
      }
    }
  }

 // End reading sqlite3 database
  sqlite3_close(db);

}


void Dracula::SetScatteringTable(char * scatfile){
  if(m_UsePohlRates == true){
    cout << "Error: Pohl rates have been set to use already to fill scattering matrix! Set it to false! Exiting..." << endl;
    exit(-1);
  }

  int ni, nf;
  double dum1;
  double temp = m_Temp; double dens = m_Density;double ei, ef, eg, el;
  double echrg = 1.602e-19 ; double ryd = 13.606*echrg ; double melec = 9.109e-31 ;
  double potfac = 8.988E9*echrg*echrg ; double kboltz = 1.38e-23 ; double pi = 2.0*asin(1.0) ;
  double eps0 = 8.854e-12 ;
  double rat = ryd/(kboltz*temp) ;
  double x = 0;
  double k0 = dens*potfac*potfac/(kboltz*temp*sqrt(melec*ryd)) ;
  double lam = 6.626e-34/sqrt(2.0*pi*melec*kboltz*temp) ;
  
  //Get ionization rates and scattering rates from Robicheax's CTMC data
  FILE *inscatrate ; inscatrate = fopen(scatfile,"r") ;
  cout << "Reading in scat.rate file " << scatfile << "... " << endl;
  if(inscatrate != NULL)
    {
      while(!feof(inscatrate))
	{
	  fscanf(inscatrate,"%i %i %lf",&ni,&nf,&dum1) ;
	  //	  printf("%i %i %13.6E \n", ni, nf, dum1);
	  if(ni < m_Nst && nf < m_Nst){
	    if(nf == 0){
	      m_kion[ni] = dens*dum1;
            }
	    if(nf != 0)m_kcol[ni][nf] = dens*dum1;
	  }
	  
	}
    }
    else {printf("!!!!hey where is the input scattering data file!!!!\n kill program and get the input data\n") ; getchar() ;}
    fclose(inscatrate) ;		   
    
    ApplyDetailedBalance();
    

    // Calculate TBR rate from ionisation rate using Saha-Boltzmann detailed balance
    for(int j = 10 ; j < m_Nactst ; j++){
      m_ktbr[j] = m_kion[j]*j*j*lam*dens*lam*lam*exp(13.606*echrg/(j*j*kboltz*temp));
      //     if(x < 39.0) dum1 += ktbr[j]*(1.0-exp(-kion[j]*2.e-4))/(kion[j]*2.e-4) ;
    }

    // Extrapolate linearly the ionization rates to higher n
    for(int j = m_Nactst-1; j < m_Nactst; j++){
      m_kion[j] = m_kion[m_Nactst-2] + (j - (m_Nactst-2))*(m_kion[m_Nactst-2] - m_kion[m_Nactst-3])/((m_Nactst-2) - (m_Nactst-3));
      m_ktbr[j] = m_kion[j]*j*j*lam*dens*lam*lam*exp(13.606*echrg/(j*j*kboltz*temp));
      
    }
    
    
}


void Dracula::SetRadRecTable(char * rrfile){

  // Fill the radiative recombination rate array from file
  // for the moment we need to convert from atomic units rate to
  // SI
  double conv = 6.12e-15;//m^{3}/s

  ifstream ifile;
  cout << "Reading in rad.reco. rate file " << rrfile << "... " << endl;
  ifile.open(rrfile);
  int nF;
  double rrateKramers, rrateLengthG, rrateVelG;
  while(ifile.good()){
    ifile >> nF >> rrateKramers >> rrateLengthG >>  rrateVelG;
    if(nF < m_Nst) m_krr[nF] = m_Density*conv*rrateKramers;
  }
  ifile.close();

}

void Dracula::SetRadDecTable(char * rdfile){
  // Fill decay rate matrix from file
  ifstream ifile;
  cout << "Reading in rad.dec. rate file " << rdfile << "... " << endl;
  ifile.open(rdfile);
  int nI, nF, lI, lF;
  double rrate;
  while(ifile.good()){
    ifile >> nI >> lI >> nF >> lF >> rrate;
    if(nI < m_Nst && nF < m_Nst)m_gamrad[nI][nF]+=rrate*(2*lI+1)/(nI*nI);
  }
  ifile.close();
}

// For T=0 we use the Pohl scattering rates
void Dracula::UsePohlRates(){
  m_UsePohlRates = true;
  cout << "Using Pohl et al's scattering rates for B=0 T..." << endl;

  //compute the ionization rate and the TBR rates from Pohl et al
  double temp = m_Temp; double dens = m_Density;double ei, ef, eg, el;
  double echrg = 1.602e-19 ; double ryd = 13.606*echrg ; double melec = 9.109e-31 ;
  double potfac = 8.988E9*echrg*echrg ; double kboltz = 1.38e-23 ; double pi = 2.0*asin(1.0) ;
  double eps0 = 8.854e-12 ;
  double rat = ryd/(kboltz*temp) ;
  double x = 0;
  double k0 = dens*potfac*potfac/(kboltz*temp*sqrt(melec*ryd)) ;
  double lam = 6.626e-34/sqrt(2.0*pi*melec*kboltz*temp) ;
  for(int j = 1 ; j < m_Nst ; j++){
    x = rat/(j*j) ;
    m_kion[j] = k0*11.0*sqrt(rat)/(pow(x,2.333333)+4.38*pow(x,1.72)+1.32*x) ;
    if(m_Temp > 10 && x > 40.0) m_kion[j] = 0.0 ;
    else m_kion[j] *= exp(-x) ;
    //      if(x < 39.0) dum1 += m_ktbr[j]*(1.0-exp(-m_kion[j]*2.e-4))/(m_kion[j]*2.e-4) ;
    //    printf("%i %13.6E %13.6E %13.6E %13.6E\n",int(j),m_kion[j],m_ktbr[j],m_ktbr[j]*(1.0-exp(-m_kion[j]*2.e-4))/(m_kion[j]*2.e-4), 1.0) ;
    m_kcol[j][j] = 0.0 ;
  }

  // Calculate TBR from Saha-Boltzman relation detailed balance
  for(int j=10; j < m_Nactst; j++){
    m_ktbr[j] = m_kion[j]*j*j*lam*dens*lam*lam*exp(13.606*echrg/(j*j*kboltz*temp));
    cout << j << "\t" << m_ktbr[j] << endl;
  }

  //compute the collision rate that takes n_i to n_f
  for(int ni = 1 ; ni < m_Nactst ; ni++)
    {
      //	m_kcol[ni][ni] = -m_kion[ni] ;
      for(int nf = 1 ; nf < m_Nst ; nf++)
	{
	  if(ni != nf)
	    {
	      ei = rat/(ni*ni) ; ef = rat/(nf*nf) ;
	      if(ei > ef) {eg = ei ; el = ef ;}
	      else {eg = ef ; el = ei ;}
	      m_kcol[ni][nf] = k0*pow((ei/eg),2.5)*ef*sqrt(ef)*exp(el-ei)*
		(22.0/pow((eg+0.9),2.333333) + 4.5/(pow(eg,2.5)*pow((eg-el),1.33333))) ;
	    }
	}
      
    }
  //m_kcol now contains all rates that take the atom from n_i to n_f
  

}

// Adding Black-Body radiation:
// - additional term for spontaneous radiative decay rate: for stimulated emission, nf < ni:
// New Rate(nf <- ni) = Old Rate(nf <- ni) x (1 + 1/[e^{hf/kB*T} -1])
// - additional new term responsible for stimulated absorption, nf > ni
// New Rate(nf <- ni) = Old Rate(ni <- nf) x (0 + 1/[e^{hf/kB*T} -1])
void Dracula::AddBlackBody(){
  
  cout << "Adding stimulated absorption/emission due to Black-Body radiation..." << endl;
  
  // loop over all radiative decay rates and modify them according to the
  // formulas
  double Eif = 0.0;
  for(int ni = 1; ni < m_Nst; ni++){
    for(int nf = 1; nf < m_Nst; nf++){
      Eif = e0*eV*fabs(1.0/(ni*ni) - 1.0/(nf*nf));
      if(nf < ni)
	m_gamrad[ni][nf] = m_gamrad[ni][nf]*(1.0+ 1.0/(exp(Eif/(kb*m_TempRad)) - 1));
      if(nf > ni)
	m_gamrad[ni][nf] += m_gamrad[nf][ni]*(1.0/(exp(Eif/(kb*m_TempRad)) - 1));

    }
  }


}

void Dracula::ApplyDetailedBalance(){

  double temp = m_Temp; double dens = m_Density;double ei, ef, eg, el;
  double echrg = 1.602e-19 ; double ryd = 13.606*echrg ; double melec = 9.109e-31 ;
  double potfac = 8.988E9*echrg*echrg ; double kboltz = 1.38e-23 ; double pi = 2.0*asin(1.0) ;
  double eps0 = 8.854e-12 ;
  double rat = ryd/(kboltz*temp) ;
  double x = 0;
  double k0 = dens*potfac*potfac/(kboltz*temp*sqrt(melec*ryd)) ;
  double lam = 6.626e-34/sqrt(2.0*pi*melec*kboltz*temp) ;

  // Apply detailed balance principle because
  //
  // n_i^2 exp(13.6 eV/(n_i^2 k_B T)) rate(n_i -> n_f) = n_f^2 exp(13.6 eV/(n_f^2 k_B T)) rate(n_f -> n_i)
  //
  // is only approximately correct. My suggestion is for all cases where
  // n_i -> n_f AND n_f -> n_i
  // are in the data set, define
  //
  // S(n_i -> n_f) = (1/2) [n_i^2 exp(13.6 eV/(n_i^2 k_B T)) rate(n_i -> n_f) + n_f^2 exp(13.6 eV/(n_f^2 k_B T)) rate(n_f -> n_i)]
  //
  // then do
  //
  // better rate(n_i -> n_f) = S(n_i -> n_f)/[n_i^2 exp(13.6 eV/(n_i^2 k_B T))]
  double S =0.0;
  for(int i = 1; i < m_Nst; i++){
    for(int j = 1; j < m_Nst; j++){
      //printf("m_kcol(%i, %i) : %f\n", i, j, m_kcol[i][j]);
      if(m_kcol[i][j] > 0.0 && m_kcol[j][i] > 0.0 && i != j && m_kcol_db[i][j] == false && m_kcol_db[j][i] == false){
	S = 0.5*(i*i*exp(13.606*echrg/(i*i*kboltz*temp))*m_kcol[i][j] + j*j*exp(13.606*echrg/(j*j*kboltz*temp))*m_kcol[j][i]);
	
	m_kcol[i][j] = S/(i*i*exp(13.606*echrg/(i*i*kboltz*temp)));
	m_kcol_db[i][j] = true;// detailed balance applied
	m_kcol[j][i] = S/(j*j*exp(13.606*echrg/(j*j*kboltz*temp)));
	m_kcol_db[j][i] = true;// detailed balance applied
      }
    }
  }
  
  
  // for(int i = 20; i < m_Nst ; i++){
  //   for(int j = 20; j < m_Nst ; j++){
  //     if( m_kcol[i][j] > 0 && m_kcol[j][i] == 0) {
  // 	// calculate the scattering rate from detailed balance
  // 	S = i*i*exp(13.606*echrg/(i*i*kboltz*temp))*m_kcol[i][j];
  // 	m_kcol[j][i] = S/(j*j*exp(13.606*echrg/(j*j*kboltz*temp)));
  // 	m_kcol_db[j][i] = true;// detailed balance applied
  // 	//printf("Using detailed balance: %i --> %i: %13.6E  (%13.6E)\n", j, i,m_kcol[j][i], m_kcol[i][j] );
  //     }
      
  //    }
  // }
  
}

void Dracula::ApplyEstabPop(){

  // Redefinition of collision rates to include
  // initial distribution and scattering:
  // : TBR + b-b scattering into the final state, for limited states n_i to n_f
  // k(0, nf) = k(0, nf) + \sum_ni P_ni * k(ni, nf)
  // : ionization + b-b scattering away from the state ni to any nf, for limited states n_i to n_f
  // k(ni, 0) = k(ni, 0) + \sum_nf k(ni, nf) 
  
  // We redefine the initial population to be the continuum which provides a constant
  // source of electrons, so they are never modified!
  // Only update the states below this limit (say 100)
  // But we have to sum up the contributions from the continuum for the redefinition
  
  
  double sum = 0.0;
  double prev = 0.0;
  for(int nf = 1; nf < m_Nactst; nf++){
    sum = 0.0;
    for(int ni = m_Nactst ; ni < m_Nst ; ni++){
      sum = sum + m_popinit[ni]*m_kcol[ni][nf];
    }
    prev = m_kcol[0][nf];
    m_kcol[0][nf] = m_kcol[0][nf] + sum;
    m_kcol[0][0] -= sum;
    //      printf("m_kcol(0, %i) before: %f, m_kcol(0, %i) after: %f, sum: %f\n", nf, prev, nf, m_kcol[0][nf], sum);
  }
  
  for(int ni = 1; ni < m_Nactst; ni++){
    sum = 0.0;
    for(int nf = m_Nactst ; nf < m_Nst ; nf++)
      sum = sum + m_kcol[ni][nf];
    
    prev = m_kcol[ni][0];
    m_kcol[ni][0] = m_kcol[ni][0] + sum;
    m_kcol[ni][ni] -= sum;
    //      printf("m_kcol(%i, 0) before: %f, m_kcol(%i, 0) after: %f, sum: %f\n", ni, prev, ni, m_kcol[ni][0], sum);
  }
  
}

void Dracula::SetPop(double * popv){
  m_useLevPop = true;
  for(int i = 0; i < m_Nactst; i++)
    m_UserLevPop[i] = popv[i];

}

// return true if success
bool Dracula::Compute(){

  // Turning off all collisions temporarily
  // for(int i = 0; i < m_Nst; i++){
  //   m_kion[i] = 0;
  //   for(int j = 0; j < m_Nst; j++){
  //     m_kcol[i][j] = 0;
  //   }
  // }

  // for(int j = 30; j >= 10; j--){
  //   cout << m_Temp << "\t" << m_Density << "\t" << j << "\t" << m_kcol[j][j-1] << "\t" << m_ktbr[j] << std::endl;
  // }



  // print nf = 0 
  // for(int i = 0; i < m_Nactst;i++){
  //   printf("m_kcol(%i, 0):  %f, m_kcol(0, %i): %f, m_krr: %f, m_ktbr: %f, m_gammrad(%i, 0): %f, m_gammrad(0, %i): %f,\n", i, m_kcol[i][0], i, m_kcol[0][i], m_krr[i], m_ktbr[i], i, m_gamrad[i][0], i, m_gamrad[0][i]);
  // }

  time_t start_time, finish_time;
  double seconds;
  time(&start_time); 


  double y[m_Nactst+1]; 
  // initialize array
  for(int i = 0; i < m_Nactst; i++)y[i] = 0;

  if(!m_useLevPop){
    // initialize to total number of protons
    y[m_Nactst] = m_NumbP;
  }else{
    for(int i = 0; i < m_Nactst; i++)y[i] = m_UserLevPop[i];
    y[m_Nactst] = 0;
  }

  // use thermal initial pop
  //  for(int i = m_Nactst-20; i < m_Nactst; i++)y[i] = GetThermPop(i);


  // the element m_Nactst+1 is responsible for the ion density
 const gsl_odeiv2_step_type * T 
   = gsl_odeiv2_step_bsimp;
 gsl_odeiv2_step * s 
   = gsl_odeiv2_step_alloc (T, m_Nactst+1);
 gsl_odeiv2_control * c 
   = gsl_odeiv2_control_y_new (1e-08, 1e-08);
 gsl_odeiv2_evolve * e 
   = gsl_odeiv2_evolve_alloc (m_Nactst+1);


  gsl_odeiv2_system sys;
  
  sys.function = this->func; sys.jacobian = this->jac; sys.dimension= m_Nactst+1; sys.params = NULL;
    
  // gsl_odeiv2_driver * d = 
  //   gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_bsimp,
  // 				  1e-09, 1e-06, 0.0);
  double t = 0.0, t1 = m_Time*1e-06;
  double h = 1e-09;

  cout << "Starting integration..." << endl;
  unsigned int nsteps = 100;
  unsigned int cstep = 0;
  while (t < t1){
    int status = gsl_odeiv2_evolve_apply (e, c, s,
					  &sys, 
					  &t, t1,
					  &h, y);
    
    if (status != GSL_SUCCESS)
      break;

    
     // nsteps = (unsigned long long)(t1/h);
     // cstep = (unsigned long long)(t/h);
     // loadbar(cstep, nsteps, 20);
    
    // For debugging

    // double Tot = 0;
    // for(int i = 39; i < m_Nactst;i++){
    //   if(i >0)Tot += y[i];
    // }

    //    printf ("%.5e %.5e %.5e\n", t, y[1], y[m_Nactst]);
    //    printf ("%20.6E\t %20.6E\t %20.6E\n", t, y[m_Nactst], y[1]);
    //    printf ("%20.6E\t %20.9E\t %20.6E\n", t, y[m_Nactst], Tot);
    
    printf ("%20.6E\t %20.6E\n", t, y[1]);
  }

  gsl_odeiv2_evolve_free (e);
  gsl_odeiv2_control_free (c);
  gsl_odeiv2_step_free (s);

  time(&finish_time);
  printf("Integration time: %.1f seconds\n", difftime(finish_time, start_time));

  double Tot = 0;
  for(int i = 0; i < m_Nactst;i++){
    m_LevPop[i] = y[i];
    if(i >0)Tot += y[i];
  }
  printf("Total number of bound states: %20.6E, burned antiprotons: %20.6E \n", Tot , m_NumbP - y[m_Nactst]);
  
  return true;
}

// to screen
void Dracula::PrintLevelPop(){

  for(int i = 0; i < m_Nactst-1;i++){
    cout << i+1 << "\t" << m_LevPop[i+1] << endl;
  }

}

// to ascii file
void Dracula::SaveLevelPop(const char * outputfile){

  cout << "Saving result to ascii file: " << outputfile << endl;
  ofstream outfile;
  outfile.open(outputfile);
  for(int i = 0; i<m_Nactst-1;i++)
    outfile << i+1 << "\t" << m_LevPop[i+1] << endl;

  outfile.close();
  

}

// to CSV file
void Dracula::SaveLevelPopCSV(const char * outputfile){

  cout << "Saving result to csv file: " << outputfile << endl;
  ofstream outfile;
  outfile.open(outputfile);
  // write header
  if(!m_UseBlackBody){
    outfile << "PrincQN, " << m_Temp << "K_" << m_Density << "/m3_" << m_Bfield << "T_" << m_Time <<  "us, Thermal-Equillibrium" << endl;
  }else{
    outfile << "PrincQN, " << m_Temp << "K_BB" << m_TempRad << "K_"<< m_Density << "/m3_" << m_Bfield << "T_" << m_Time <<  "us, Thermal-Equillibrium" << endl;
  }
  for(int i = 0; i<m_Nactst-1;i++){
    // this is needed for plotting in webpage
    if(i < 20 && m_Temp <= 50){
      outfile << i+1 << ", " << m_LevPop[i+1] << ", NaN" << endl;
    }else if(i < 10 && m_Temp > 50){
      outfile << i+1 << ", " << m_LevPop[i+1] << ", NaN" << endl;
    }else{
      outfile << i+1 << ", " << m_LevPop[i+1] << ", " << GetThermPop(i+1) << endl;
    }
  }
  outfile.close();
}


// for parameterscan
// append parameters and population in quantum state i
void Dracula::SavePop(int i, const char * outputfile){

 cout << "Saving result to file: " << outputfile << endl;
  ofstream outfile;
  outfile.open(outputfile, ios_base::app);

  for(int ii = 1; ii<m_Nactst;ii++){
    // this is needed for plotting in webpage
    if(i == ii){
      outfile << i << "\t" << m_LevPop[i] <<  "\t"  << m_Temp << "\t" << m_Density << "\t" << m_TempRad << "\t" << m_Bfield << "\t" << m_Time << "\t" << m_NumbP << endl;
    }
  }
  outfile.close();


}


double Dracula::GetThermPop(int n){
  double kbT_eV = kb*m_Temp/eV;
  double p = m_NumbP*n*n*m_Density*pow((h*h/(2*pi*me*kb*m_Temp)), 3.0/2.0)*exp(e0/(n*n*kbT_eV));
  return p;

}

  

int Dracula::func (double t, const double y[], double f[],
	  void *params){
  
  // First for bound states
  // and bookkeeping of losses for ions
  double ionloss = 0;
  double iongain = 0;
  for(int i = 1; i < m_Nactst;i++){
    f[i] = m_krr[i]*y[m_Nactst] + m_ktbr[i]*y[m_Nactst] - m_kion[i]*y[i];
    //f[i] = - m_kion[i]*y[i];
    ionloss += (m_krr[i]*y[m_Nactst] + m_ktbr[i]*y[m_Nactst]);
    iongain += m_kion[i]*y[i];
    double integ = 0;
    for(int j = 1; j < m_Nactst; j++){
      if(j!=i){
	integ += (m_kcol[j][i] + m_gamrad[j][i])*y[j] - (m_kcol[i][j] + m_gamrad[i][j])*y[i];
      }
    }
    f[i] += integ;
  }
  // For ion state
  f[m_Nactst] = -ionloss + iongain;
  
  return GSL_SUCCESS;
}


int Dracula::jac (double t, const double y[], double *dfdy, 
	 double dfdt[], void *params){

 gsl_matrix_view dfdy_mat 
    = gsl_matrix_view_array (dfdy, m_Nactst+1, m_Nactst+1);
  gsl_matrix * m = &dfdy_mat.matrix;

  double integ = 0.0;
  double integ2 = 0.0;
  for(int i = 0; i < m_Nactst+1; i++){
    dfdt[i] = 0.0;
    integ = 0.0;integ2 = 0.0;   
     for(int j = 0; j < m_Nactst+1; j++){
       if(i == 0 && j == 0){
	 gsl_matrix_set (m, i, j, 0.0);
       }else{
	 if(i != j){
	   if(i != m_Nactst && j == m_Nactst){
	     gsl_matrix_set (m, i, j, m_krr[i] + m_ktbr[i]);
	   }else if(i == m_Nactst && j != m_Nactst){
	     gsl_matrix_set (m, i, j, m_kion[j]);
	     integ2 += m_krr[j] + m_ktbr[j];
	   }else{
	     gsl_matrix_set (m, i, j, m_kcol[j][i] + m_gamrad[j][i]);
	     integ += m_gamrad[i][j] + m_kcol[i][j];
	   }
	 }
	 if(i == j){
	   if(i == m_Nactst){
	     gsl_matrix_set (m, i, j, -(m_krr[i] + m_ktbr[i]));
	   }else{
	     //	     gsl_matrix_set (m, i, j, -m_kcol[i][j] - m_gamrad[i][j] - m_kion[i]);
	   }
	 }
       }
     }
     if(i == m_Nactst){
       gsl_matrix_set (m, i, i, -integ2);
     }
     if(i != m_Nactst)
       gsl_matrix_set (m, i, i, -integ - m_kion[i]);     
  }


  return GSL_SUCCESS;

}

char * Dracula::ftoa(float x){
  char * res = new char[128];
  sprintf(res, "%f", x);
  return res;
}

char * Dracula::Etoa(double x){
  char * res = new char[128];
  sprintf(res, "%E", x);
  return res;
}


char * Dracula::itoa(int x){
  char * res = new char[128];
  sprintf(res, "%d", x);
  return res;
}

static inline void loadbar(unsigned long long x, unsigned long long n, unsigned int w = 50)
{

  if ( (x != n) && (x % (n/100) != 0) ) return;

  float ratio  =  x/(float)n;
  int   c      =  ratio * w;
  
  cout << setw(3) << (int)(ratio*100) << "% [";
  for (int x=0; x<c; x++) cout << "=";
  for (int x=c; x<w; x++) cout << " ";
  cout << "]\r" << flush;
}
