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
#include "Dracula.h"


char* getCmdOption(char ** begin, char ** end, const std::string & option)
{
    char ** itr = std::find(begin, end, option);
    if (itr != end && ++itr != end)
    {
        return *itr;
    }
    return 0;
}

bool cmdOptionExists(char** begin, char** end, const std::string& option)
{
    return std::find(begin, end, option) != end;
}

void printUsage(){

  cout << "Usage: <program> [-h] [-f sqlitefile -t positrontemp -d positrondens -an antiproton_number -b magfield -ot overlaptime -otp [ascii|csv] -of outputfile [-r blackbodytemp]] " << endl;

}

int main(int argc, char * argv[]){

  if(argc == 1){
    printUsage();
    return 0;
  }


 if(argc == 2 && cmdOptionExists(argv, argv+argc, "-h")){
   printUsage();
   return 0;
 }

 char * filename = NULL;
 char * outfilename = NULL;
 string ascii_type = string("ascii");
 string csv_type = string("csv");
 string outfiletype = ascii_type;
 int pos_temp = 200;
 double pos_dens = 1e+14;
 double overlap_time = 5;
 double magfield = 1.0;
 double rad_temp = 300;
 double antiproton_n = 1e+05;

 if(cmdOptionExists(argv, argv+argc, "-f") && 
    cmdOptionExists(argv, argv+argc, "-t") &&
    cmdOptionExists(argv, argv+argc, "-d") &&
    cmdOptionExists(argv, argv+argc, "-b") &&
    cmdOptionExists(argv, argv+argc, "-an") &&
    cmdOptionExists(argv, argv+argc, "-ot")&&
    cmdOptionExists(argv, argv+argc, "-of")){
    
   filename = getCmdOption(argv, argv + argc, "-f");
   cout << "SQLite filename : " << filename << endl;
   pos_temp = atoi(getCmdOption(argv, argv + argc, "-t"));
   cout << "positron temp : " << pos_temp << " K " << endl;
   pos_dens = atof(getCmdOption(argv, argv + argc, "-d"));
   cout << "positron dens : " << pos_dens << " /m3 " << endl;
   antiproton_n = atof(getCmdOption(argv, argv + argc, "-an"));
   cout << "number of antiprotons : " << antiproton_n << " " << endl;
   magfield = atof(getCmdOption(argv, argv + argc, "-b"));
   cout << "magnetic field : " << magfield << " T " << endl;
   overlap_time = atof(getCmdOption(argv, argv + argc, "-ot"));
   cout << "overlap time : " << overlap_time << " us " << endl;
   outfilename = getCmdOption(argv, argv + argc, "-of");
   
   if(cmdOptionExists(argv, argv+argc, "-otp")){
     outfiletype = string(getCmdOption(argv, argv + argc, "-otp"));
     if(!(outfiletype == ascii_type || outfiletype == csv_type)){
       cout << "Error: wrong outputfile type <" << outfiletype << ">, please use either \"ascii\" OR \"csv\" (default is ascii)" << endl;
       printUsage();
       return 0;
     }
   }	  


   // Instantiate class
   Dracula d(150);

   // First set all parameters
   d.SetPositronTemp(pos_temp);//K
   d.SetPositronDensity(pos_dens);// /m3
   d.SetOverlapTime(overlap_time);// us
   d.SetBField(magfield);
   d.SetAntiprotonN(antiproton_n);
   if(cmdOptionExists(argv, argv+argc, "-r")){
       rad_temp = atof(getCmdOption(argv, argv + argc, "-r"));
       cout << "black-body radiation temp : " << rad_temp << " K " << endl;
       d.SetBlackBodyTemp(rad_temp);
     }
     
     // Now read in the rates
     // and set up rate matrices
     d.SetSQLDbScat(filename);
     d.SetSQLDbRadDec(filename);
     
     // Start integration
     d.Compute();
     
     // Print/Save result
     if(outfiletype == csv_type){
       d.SaveLevelPopCSV(outfilename);
     }
     if(outfiletype == ascii_type){
       d.SaveLevelPop(outfilename);
     }
     //d.PrintLevelPop();
     d.SavePop(1, "out.txt");

 }else{
   cout << "Error: Not enough parameters!" << endl;
   printUsage();
   return 0;
 }
 
   return 0;
 }
