// To compile:
// g++ -O0 -fbounds-check corsikaReader.cpp -o corsikaReader -std=c++11 -lm
//
// Or can run command "make" in file directory (as long as accompanied Makefile is also in directory...)
//
// Usage (after compilation):
// ./corsikaReader PATH_TO_CORSIKA_OUTPUT_FILE --FILE_FLAG
//
// Output (printed to terminal):
// PRIMARY_ID PRIMARY_ENERGY ZENITH AZIMUTH MUON_NUM_TOT MUON_NUM_ECUT

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <algorithm>
#include <math.h>
#include <bitset>
#include <climits>
#include <sstream>
#include <iomanip>
using namespace std;
#include <glob.h>


bool getBinary(float g, bool thinned) {
  union
  {
    float input; // assumes sizeof(float) == sizeof(int)
    int   output;
  } data1;
  union
  {
    float input; // assumes sizeof(float) == sizeof(int)
    int   output;
  } data1a;
  union
  {
    float input; // assumes sizeof(float) == sizeof(int)
    int   output;
  } data2;

  if (thinned) {
    data1.input = 3.67252e-41; // must be this for thinned files
  } else {
    data1.input = 3.21346e-41; // must be this for regular files
  }

  data1a.input = 4.59037e-41;
  data2.input = g;

  std::bitset<sizeof(float) * CHAR_BIT> bits1(data1.output);
  std::bitset<sizeof(float) * CHAR_BIT> bits1a(data1a.output);
  std::bitset<sizeof(float) * CHAR_BIT> bits2(data2.output);
  if ((bits1 == bits2) || (bits1a == bits2)) {
    return true;
  }

  return false;
}

// Used for defining the type of corsika simulation
enum class SimType {Thinned, Standard};

/// --------------------------------------------------------------------------------------------
/// MAIN PART - READING.....
/// --------------------------------------------------------------------------------------------


int main (int argc, char *argv[]) {

    std::string primary = "iron";   // e.g. "proton"
    std::string energy = "lgE_18.0";    // e.g. "lgE_16.0"
    std::string theta = "sin2_0.4";     // e.g. "sin2_0.4"


  if (argc < 3) {
    cerr << "------------------------------------------------------\n";
    cerr << "This program counts the muons in the air shower and prints shower information to terminal:\n";
    cerr << "You must give in the input filename and type of CORSIKA file (thinned or standard)\n";
    cerr << "Usage is ./corsikaReader <InputFileName> --FILE_FLAG\n";
    cerr << "--FILE_FLAG can be: --thinned or --standard\n";
    cerr << "------------------------------------------------------\n";

    return 0;
  }

  std::string filePath = argv[1];
  std::string fileFlag = argv[argc - 1];
  //std::string fileFlag = "--thinned";
  
  SimType mode;

  if (fileFlag == "--thinned") {
    mode = SimType::Thinned;   // thinned corsika file
  } else if (fileFlag == "--standard") {
    mode = SimType::Standard;   // standard corsika file
  } else {
    cerr << "------------------------------------------------------\n";
    cerr << "Invalid file flag given!\n";
    cerr << "Usage is ./corsikaReader <InputFileName> --FILE_FLAG\n";
    cerr << "--FILE_FLAG must be either: --thinned or --standard\n";
    cerr << "------------------------------------------------------\n";
    return 0;
  }

  // Ternary operations
  // If mode is Thinned, then use "thinned corsika" record size, else use "standard corsika" record size
  const int nrecstd = (mode == SimType::Thinned) ? 26216 : 22940;
  const int nsblstd = (mode == SimType::Thinned) ? 312 : 273;

  // Other constants
  const int numbstd = nrecstd / 4;     // = 6554 for "thinned corsika", = 5735 for "standard corsika"
  float sdata[numbstd];                // to read data for a single corsika record

  // Constant for ternary operation to define particle weights in data block
  const bool isThin = (mode == SimType::Thinned) ? true : false;

  vector<string> possible_headers = {"RUNH", "EVTH", "LONG", "EVTE", "RUNE"};

  glob_t glob_result;
  glob(filePath.c_str(), GLOB_TILDE, NULL, &glob_result);

  /// init variables
  bool BROKENflag = false;
  int EVTEcnt = 0;
  int nrShow = 0;
  float primaryID, primaryEnergy; 
  double zenith, azimuth;
  primaryID = 0.;
  primaryEnergy = 0.;
  zenith = 0.;
  azimuth = 0.;

  
  cout << "RunNumber, PrimaryID, PrimaryEnergy, Zenith, Azimuth, 2ndParticle, TotalParticleNumber" << endl;


  std::vector<std::string> particles = {"mupm", "epm"};


  /// --------------------------------------------------------------------------------------------
  /// THE MAIN LOOP
  /// --------------------------------------------------------------------------------------------
  /// This reads all files one by one
  //~ for(auto t=fileList.begin(); t!=fileList.end(); ++t) {

   /// Grab only MUONS and ELECTRONS, see corsika user manual for particle id's 
  for (const auto& particle : particles) { 
  int idpa1 = -1;
  int idpa2 = -1;
      
    if (particle == "mupm") {
        idpa1 = 5;
        idpa2 = 6;
                    }
    else if (particle == "epm") {
        idpa1 = 2;
        idpa2 = 3;
                    }
      
    std::ostringstream oss2;
    oss2 << "Particles/"
    << primary << "/"
    << energy << "/"
    << theta << "/TOTAL_"
    << particle << ".dat";

    std::ofstream fout2(oss2.str());
    fout2 << "#RunNumber, PrimaryID, PrimaryEnergy, Zenith, Azimuth, 2ndParticle, TotalParticleNumber\n";
        
  for (int k = 1; k < argc - 1; ++k) {
    EVTEcnt = 0;
    BROKENflag = false;

    std::string file_ = argv[k];

     // ---- extract run number ----
    size_t pos = file_.rfind("DAT");
        if (pos == std::string::npos) {
            cerr << "Cannot find DAT in filename: " << file_ << endl;
            continue;
        }
        
        pos += 3; // move past "DAT"
        
        std::string digits;
        while (pos < file_.size() && std::isdigit(file_[pos])) {
            digits.push_back(file_[pos]);
            ++pos;
        }
        
        if (digits.empty()) {
            cerr << "No run number found in filename: " << file_ << endl;
            continue;
        }

    int runnum = std::atoi(digits.c_str());

    // ---- create output file ----
    std::ostringstream oss;
    oss << "Particles/"
    << primary << "/"
    << energy << "/"
    << theta << "/DAT"
    << std::setw(6) << std::setfill('0') << runnum
    << "_" << particle << ".dat";

    std::ofstream fout(oss.str());
    fout << "# particleID px py pz x y t weight\n";
      
    // ---- read corsika file ----
    if (!(file_.find(".long") != std::string::npos)) {

      ifstream is (file_, ifstream::binary);
      // cerr << "fileName -> " << file_ << endl;

      float nParticles = 0.;
      // float nMuons300 = 0.;

      /// Read block = record --------------------------------------------------------
      while ( is.read((char*)&sdata, sizeof(sdata)) ) { /// get full block of data at once
        if ( !getBinary( sdata[0], isThin ) ) { /// skip the first record length sdata[0]
          //cout << sdata[21 * nsblstd + 1] << endl;
          cerr << "This file is corrupted, this is not a record length - beginning of block!" << endl;
          BROKENflag = true;
          break;
        }
        /// iterate over 21 sub block inside this block
        for (int j = 0; j < 21; j++) {
          string head_word = (string) (char *) &sdata[j * nsblstd + 1];
          head_word = head_word.substr (0, 4); /// dirty hack!
          if ( find( possible_headers.begin(), possible_headers.end(), head_word ) != possible_headers.end() ) {
            if (head_word == "RUNH") {
              nrShow = sdata[j * nsblstd + 93];
            } else if (head_word == "EVTH") {
              ///  Reading primary type and energy
              primaryID = sdata[j * nsblstd + 1 + 2];
              primaryEnergy = sdata[j * nsblstd + 1 + 3];
              zenith = sdata[j * nsblstd + 11];
              azimuth = sdata[j * nsblstd + 12];
              // cout << primaryID << " " << primaryEnergy << " " << zenith << " " << azimuth << " ";
            } else if (head_word == "EVTE") {
              EVTEcnt += 1;
            }
          }
          else { /// READ DATA -> iterate every 7th position
            for (int i = j * nsblstd + 1; i <= (j * nsblstd + nsblstd); i += 8 ) {
              float particle_id = sdata[i];
              int idpa =  (int)particle_id / 1000;

             
              // // NOTE: corsika uses units of GeV for momentum and mass
                // // (see user manual, section 7, pg. 119 in version 7.7410)
              if ( idpa == idpa1 || idpa == idpa2 ) {float px = sdata[i + 1]; // momentum component in x direction (units of GeV)
                                              float py = sdata[i + 2]; // momentum component in y direction (units of GeV)
                                              float pz = sdata[i + 3]; // momentum component in z direction (units of GeV)
                                              float x  = sdata[i + 4]; // Can check these definitions in user manual
                                              float y  = sdata[i + 5];
                                              float t  = sdata[i + 6];
                                              float w  = isThin ? sdata[i+7] : 1.0;
                                              /// Count total number of muons at corsika observation level
                                              nParticles += w;

                                              // double massMu = 0.105658357;  // muon mass (units of GeV)
                                                // /// Kinetic energy (units of GeV) !!!
                                                // double ekinMu = sqrt ( px * px + py * py + pz * pz + massMu * massMu) - (massMu) ;

                                                 // /// Count number of muons at observation level w/ Kinetic Energy > 300 GeV
                                                // if ( ekinMu > 300. ) {
                                                //   nMuons300 += w;
                                                // }

                                              fout << particle_id << " "
                                                     << px << " " << py << " " << pz << " "
                                                     << x  << " " << y  << " " << t  << " "
                                                     << w  << "\n";

                
              }
            }
          }
        }
          
        /// end of the record
        if ( !getBinary( sdata[21 * nsblstd + 1], isThin ) ) {
          cerr << "This file is corrupted, this is not a record length - end of block!" << endl;
          BROKENflag = true;
          break;
        }
      }

      // Round the number of muons to the nearest integer (since weights can be fractional)
      nParticles = round(nParticles);
      // nMuons300 = round(nMuons300);

      cout << runnum << " " << primaryID << " "<< primaryEnergy << " " << zenith << " " << azimuth << " " << particle << " " << nParticles << "\n";
      fout2 << runnum << " " << primaryID << " "<< primaryEnergy << " " << zenith << " " << azimuth << " " << nParticles << "\n";
      if ( BROKENflag || !(EVTEcnt == nrShow) ) {
        cerr << "Files is broken: not enough EVTE or garbage word is wrong " << file_ << endl;
        break;
      }
      is.close();
    }
      fout.close();   
  }
  fout2.close();
  }
   return 0; 
}




