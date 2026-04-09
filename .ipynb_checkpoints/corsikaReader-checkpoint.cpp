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
#include "cnpy.h"
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


// --------------------------------
// MAIN
int main(int argc, char *argv[]) {

    if (argc < 3) {
        cerr << "Usage: ./corsikaReader <InputFileName> --thinned|--standard\n";
        return 0;
    }

    string filePath = argv[1];
    string fileFlag = argv[argc - 1];

    SimType mode;
    if      (fileFlag == "--thinned") mode = SimType::Thinned;
    else if (fileFlag == "--standard") mode = SimType::Standard;
    else {
        cerr << "Invalid file flag!\n";
        return 0;
    }

    const int nrecstd = (mode == SimType::Thinned) ? 26216 : 22940;
    const int nsblstd = (mode == SimType::Thinned) ? 312 : 273;
    const int numbstd = nrecstd / 4;
    float sdata[numbstd];
    bool isThin = (mode == SimType::Thinned);

    vector<string> possible_headers = {"RUNH", "EVTH", "LONG", "EVTE", "RUNE"};

    vector<string> particles = {"mupm", "epm"};

    string primary = "proton";
    string energy = "lgE_18.0";
    string theta = "sin2_0.2";

    cout << "corsikaReader " << primary << " " << energy << " " << theta << endl;

    // Create output directories
    string out_dir = "Particles/" + primary + "/" + energy + "/" + theta;

    // Loop over particles
    for (const auto& particle : particles) {
        int idpa1 = -1, idpa2 = -1;
        if (particle == "mupm") { idpa1 = 5; idpa2 = 6; }
        else if (particle == "epm") { idpa1 = 2; idpa2 = 3; }

        ofstream fout_total(out_dir + "/TOTAL_" + particle + ".dat");
        fout_total << "#RunNumber PrimaryID PrimaryEnergy Zenith Azimuth TotalNumber\n";

        // Loop over input files
        for (int k = 1; k < argc - 1; ++k) {

            string file_ = argv[k];

            // Extract run number
            size_t pos = file_.rfind("DAT");
            if (pos == string::npos) { cerr << "Cannot find DAT in filename: " << file_ << endl; continue; }
            pos += 3;
            string digits;
            while (pos < file_.size() && isdigit(file_[pos])) { digits.push_back(file_[pos]); ++pos; }
            if (digits.empty()) { cerr << "No run number found: " << file_ << endl; continue; }

            int runnum = stoi(digits);

            // Clear vectors
            vector<float> particle_id_v, px_v, py_v, pz_v, x_v, y_v, t_v, w_v;
            float primaryID = 0, primaryEnergy = 0;
            double zenith = 0, azimuth = 0;
            bool BROKENflag = false;
            int EVTEcnt = 0, nrShow = 0;
            float nParticles = 0;

            ifstream is(file_, ifstream::binary);
            if (!is.is_open()) { cerr << "Cannot open file: " << file_ << endl; continue; }

            while (is.read((char*)&sdata, sizeof(sdata))) {
                if (!getBinary(sdata[0], isThin)) { BROKENflag = true; break; }

                // iterate over sub-blocks
                for (int j = 0; j < 21; ++j) {
                    string head_word((char*)&sdata[j*nsblstd + 1]);
                    head_word = head_word.substr(0,4);

                    if (find(possible_headers.begin(), possible_headers.end(), head_word) != possible_headers.end()) {
                        if (head_word == "RUNH") nrShow = sdata[j*nsblstd + 93];
                        else if (head_word == "EVTH") {
                            primaryID = sdata[j*nsblstd + 3];
                            primaryEnergy = sdata[j*nsblstd + 4];
                            zenith = sdata[j*nsblstd + 11];
                            azimuth = sdata[j*nsblstd + 12];
                        }
                        else if (head_word == "EVTE") EVTEcnt++;
                    }
                    else {
                        // particle data
                        for (int i = j*nsblstd + 1; i <= j*nsblstd + nsblstd; i += 8) {
                            int idpa = (int)sdata[i]/1000;
                            if (idpa == idpa1 || idpa == idpa2) {
                                float px = sdata[i+1], py = sdata[i+2], pz = sdata[i+3];
                                float x = sdata[i+4], y = sdata[i+5], t = sdata[i+6];
                                float w = isThin ? sdata[i+7] : 1.0;

                                nParticles += w;
                                particle_id_v.push_back(sdata[i]);
                                px_v.push_back(px); py_v.push_back(py); pz_v.push_back(pz);
                                x_v.push_back(x); y_v.push_back(y); t_v.push_back(t); w_v.push_back(w);
                            }
                        }
                    }
                }

                if (!getBinary(sdata[21*nsblstd + 1], isThin)) { BROKENflag = true; break; }
            }

            is.close();
            nParticles = round(nParticles);
            fout_total << runnum << " " << primaryID << " " << primaryEnergy << " " << zenith << " " << azimuth << " " << nParticles << "\n";

            if (BROKENflag || EVTEcnt != nrShow) {
                cerr << "File broken or EVTE mismatch: " << file_ << endl;
                continue;
            }

            // Save .npz **once per file**
            string filename = out_dir + "/DAT" + string(6 - to_string(runnum).length(), '0') + to_string(runnum) + "_" + particle + ".npz";
            if (!particle_id_v.empty()) {
                cnpy::npz_save(filename, "particle_id", particle_id_v.data(), {particle_id_v.size()}, "w");
                cnpy::npz_save(filename, "px", px_v.data(), {px_v.size()}, "a");
                cnpy::npz_save(filename, "py", py_v.data(), {py_v.size()}, "a");
                cnpy::npz_save(filename, "pz", pz_v.data(), {pz_v.size()}, "a");
                cnpy::npz_save(filename, "x",  x_v.data(),  {x_v.size()},  "a");
                cnpy::npz_save(filename, "y",  y_v.data(),  {y_v.size()},  "a");
                cnpy::npz_save(filename, "t",  t_v.data(),  {t_v.size()},  "a");
                cnpy::npz_save(filename, "w",  w_v.data(),  {w_v.size()},  "a");
            }

        } // end loop over files
        fout_total.close();
    } // end loop over particles

    cout << "done." << endl;
    return 0;
}




