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



#include <algorithm>
#include <string>
#include <sstream>
#include <vector>
#include "Dracula.h"

// Example code for steering the Dracula computation from c++.
// - In this example an external input file is read in from the arguments,
// which contains the resulting level population from an earlier calculation.
// - Then letting the level population evolve with a radiative decay and stimulated 
// transitions (from from room temprature radiation)

int main(int argc, char * argv[]){

// read in level population from an external output file from a previous calculation
 char * filename = argv[1];
 char * outfilename = NULL;
 int pos_temp = 300;
 double pos_dens = 1e+14;
 double overlap_time = 10;
 double magfield = 2.0;
 double rad_temp = 300;
 double antiproton_n = 1e+06;

 // flight
 double Thbar = 300.0;//K
 double v = sqrt(2*kb*Thbar/mp);// m/s
 double L = 1.5; //m
 // double flightt = (L/v)/1.e-06;//us
 double flightt = 1000;//us
 // cout << "flight distance, time: " << L << " m, " << flightt << " us" << endl;
 cout << "flight time: "  << flightt << " us" << endl;


  ifstream infile;
  infile.open(argv[2]);
  int n;vector<double> n_v;
  double pop; vector<float> pop_v;
  vector <vector <string> > data;

  int count = 0;
  while(infile){
    string s;
    if (!getline( infile, s )) break;
    istringstream ss( s );
    vector <string> record;
    while (ss)
      {
	string s;
	if (!getline( ss, s, '\t' )) break;
	record.push_back( s );
      }
    
    data.push_back( record );
  }

 if (!infile.eof())
  {
    cerr << "Fooey!\n";
  }

 infile.close();

  for(int i = 1; i < data.size(); i++){
    vector<string> reci = data[i];
    //    cout << reci[0] << "\t" << reci[1] << "\t" << reci[2] << endl;
    //    cout << reci[0] << "\t" << reci[1] <<  endl;
    n_v.push_back(atof(reci[0].c_str()));
    pop_v.push_back(atof(reci[1].c_str()));
    
  }

  Dracula d2(pop_v.size());
  d2.SetPositronTemp(pos_temp);//K
  d2.SetPositronDensity(0.0);// /m3
  d2.SetOverlapTime(flightt);// us
  d2.SetBField(0.0);
  d2.SetAntiprotonN(0);
  d2.SetBlackBodyTemp(rad_temp);
  d2.SetSQLDbRadDec(filename, false);

  double * resOL_v = new double[pop_v.size()];
  for(int i = 0; i < pop_v.size(); i++)
    resOL_v[i] = 0;
  for(int i = 0; i < pop_v.size(); i++){
    resOL_v[(int)n_v[i]] = pop_v[i];
    //   cout << n_v[i] << "\t" << resOL_v[(int)n_v[i]] <<  endl;
  }


 // add initial level population
 d2.SetPop(resOL_v);

 // Start integartion
 d2.Compute();

 // Save result
 string oname = string(argv[2]) + "_flight.dat";
 d2.SaveLevelPop(oname.c_str());

 
 return 0;
}
