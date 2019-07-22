
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <boost/format.hpp>

// How to compile
// On QCL(2019): g++ -I/usr/local/gcc73.mvapich2.el75/include/boost-1_62 -L/usr/local/gcc73.mvapich2.el75/lib -L/usr/local/gcc73.mvapich2.el75/lib64 ./geninp.cc
// On Hasebe: g++ -I/home/masaaki/include/boost ./geninp2019.cc

using namespace std;

typedef std::map<double, std::vector<std::string> > len_files;  // <bond length, file_name_scf, file_name_casscf, file_name_mrci>

std::string generate_input(const double intern_dist, const std::string name_base);

//const std::string basis("def2-SVP");
const std::string basis("def2-TZVPP");

int _doCEPT2(0);
double _tCutPNO(0.0);
double _tCutPairs(0.0);

int main(int argc, char* argv[])
{
  if(argc >= 2)
    _doCEPT2   = atoi(argv[1]);
  if(argc >= 3)
    _tCutPNO   = atof(argv[2]);
  if(argc >= 4)
    _tCutPairs = atof(argv[3]);

  std::cout << " ++ Wave function Ansatz ... " << (_doCEPT2 ? "CEPT2" : "CASPT2")       << std::endl;
  std::cout << " ++ TCutPNOs used        ... " << boost::format("%16.10e") % _tCutPNO   << std::endl;
  std::cout << " ++ TCutPairs used       ... " << boost::format("%16.10e") % _tCutPairs << std::endl;

  // Base of the name of input file
  std::string n_base("R_eq");
  {
    std::ostringstream pno;
    pno << boost::format("%e") % _tCutPNO;
    std::ostringstream pairs;
    pairs << boost::format("%e") % _tCutPairs;
    n_base += (_doCEPT2 ? "_CEPT2_" : "_CASPT2_") + pno.str() + "_" + pairs.str() + "_" + basis;
  }
  
  // Bond lengths
  std::vector<double> lengths;
  //lengths.push_back(0.6 );
  ///lengths.push_back(0.7 );
  ///lengths.push_back(0.8 );
  ///lengths.push_back(0.9 );
  lengths.push_back(0.95);
  lengths.push_back(1.0 );
  lengths.push_back(1.1 );
  lengths.push_back(1.2 );
  lengths.push_back(1.3 );
  lengths.push_back(1.4 );
  lengths.push_back(1.6 );
  lengths.push_back(1.8 );
  lengths.push_back(2.0 );
  lengths.push_back(2.2 );
  lengths.push_back(2.4 );
  lengths.push_back(2.6 );
  ///lengths.push_back(2.8 );
  ///lengths.push_back(3.0 );

  len_files names;

  for(std::vector<double>::const_iterator l = lengths.begin();l != lengths.end();++l){
    std::ostringstream stm;
    stm << *l;
    std::string scf_name   ("./scf_r_"   + stm.str() + "_.out");
    std::string ccsd_name  ("./ccsd_r_"  + stm.str() + "_.out");
    std::string ccsdt_name ("./ccsdt_r_" + stm.str() + "_.out");    
    std::vector<std::string> file_names;
    file_names.push_back(scf_name);
    file_names.push_back(ccsd_name);
    file_names.push_back(ccsdt_name);
    names.insert(len_files::value_type(*l, file_names));

    std::string name(generate_input(*l, n_base));
    std::string orca_exec  ("~/orca_kofo_2qcl/orca/x86_exe/orca ./" + name + " | tee " + scf_name);
    system(orca_exec.c_str());
  }

  {
    
    std::string fileName("orcascf_summary");
    {
      std::ostringstream pno;
      pno << boost::format("%e") % _tCutPNO;
      std::ostringstream pairs;
      pairs << boost::format("%e") % _tCutPairs;
      fileName += pno.str() + "_" + pairs.str() + "_" + basis + ".out";
    }
    
    // Process output files of scf procedure
    std::ofstream outfile(fileName);

    // Extract information from the output file to compile plottable data file
    for(len_files::const_iterator n = names.begin();n != names.end();++n){
      std::ifstream infile(n->second[0].c_str());
      if(infile.fail()){
        cout << "Can't open " + n->second[0] << endl;
        abort();
      }

      // Search "total scf energy"
      std::string energy;
      std::string line;
      bool yes(false);

      bool isE0(false);
      bool isDots(false);
      bool contiguous(false);
      while(!infile.eof()){
        infile >> line;
        size_t where;
        if(yes) { energy = line; break; }
        else if(!isE0 && !isDots){
          where = line.find("E(0)", 0);
          if(where != std::string::npos)
            { isE0 = true; contiguous = true; }
          else contiguous = false;
        }
        else if(isE0 && !isDots){
          where = line.find("...", 0);
          if(where != std::string::npos)
            { isDots = true; contiguous = true; }
          else contiguous = false;
        }
        else contiguous = false;

        if(isE0 && isDots) { yes = true; /*outfile << where << " llll " << line << endl;*/ }
        else yes = false;
        if(!contiguous) { isE0 = false; isDots = false; }
      }

      if(energy.size())
        outfile << boost::format("%5.4f  %10s") % n->first % energy << endl;
      else {
        outfile << boost::format("# %5.4f  energy not found") % n->first << endl;
        cout << "scf energy can't be found for " << n->first << endl;
      }

    } // End n
  } // End scope

  {
    std::string fileName("orcaccsd_summary");
    {
      std::ostringstream pno;
      pno << boost::format("%e") % _tCutPNO;
      std::ostringstream pairs;
      pairs << boost::format("%e") % _tCutPairs;
      fileName += pno.str() + "_" + pairs.str() + "_" + basis + ".out";
    }
    
    // Process output files of casscf procedure
    std::ofstream outfile(fileName);

    // Extract information from the output file to compile plottable data file
    for(len_files::const_iterator n = names.begin();n != names.end();++n){
      std::ifstream infile(n->second[1].c_str());
      if(infile.fail()){
        cout << "Can't open " + n->second[1] << endl;
        abort();
      }

      // Search "CASSCF energy :"
      std::string energy;
      std::string line;
      bool yes(false);

      bool isECC(false);
      bool isDots(false);
      bool contiguous(false);
      while(!infile.eof()){
        infile >> line;
        size_t where;
        if(yes) { energy = line; break; }
        else if(!isECC && !isDots){
          where = line.find("E(CCSD)", 0);
          if(where != std::string::npos)
            { isECC = true; contiguous = true; }
          else contiguous = false;
        }
        else if(isECC && !isDots){
          where = line.find("...", 0);
          if(where != std::string::npos)
            { isECC = true; contiguous = true; }
          else contiguous = false;
        }
        else contiguous = false;

        if(isECC && isDots) { yes = true; /*outfile << where << " llll " << line << endl;*/ }
        else yes = false;
        if(!contiguous) { isECC = false; isDots = false; }
      }

      if(energy.size())
        outfile << boost::format("%5.4f  %10s") % n->first % energy << endl;
      else {
        outfile << boost::format("#%5.4f  energy not found") % n->first << endl;
        cout << "casscf energy can't be found for " << n->first << endl;
      }

    } // End n
  } // End scope

  {

    std::string fileName("orcaccsdt_summary");
    {
      std::ostringstream pno;
      pno << boost::format("%e") % _tCutPNO;
      std::ostringstream pairs;
      pairs << boost::format("%e") % _tCutPairs;
      fileName += pno.str() + "_" + pairs.str() + "_" + basis + ".out";
    }
    
    // Process output files of mrci procedure
    std::ofstream outfile(fileName);

    // Extract information from the output file to compile plottable data file
    for(len_files::const_iterator n = names.begin();n != names.end();++n){
      std::ifstream infile(n->second[2].c_str());
      if(infile.fail()){
        cout << "Can't open " + n->second[2] << endl;
        abort();
      }

      // Search "FINAL TOTAL ENERGY       ..."
      std::string energy;
      std::string line; 
      bool yes(false);
      
      bool isECC(false);
      bool isDots(false);
      bool contiguous(false);
      while(!infile.eof()){
        infile >> line;
        size_t where;
        if(yes) { energy = line; break; }
        else if(!isECC && !isDots){
          where = line.find("E(CCSD(T))", 0);
          if(where != std::string::npos)
            { isECC = true; contiguous = true; }
          else contiguous = false;
        }
        else if(isECC && !isDots){
          where = line.find("...", 0);
          if(where != std::string::npos)
            { isECC = true; contiguous = true; }
          else contiguous = false;
        }
        else contiguous = false;

        if(isECC && isDots) { yes = true; /*outfile << where << " llll " << line << endl;*/ }
        else yes = false;
        if(!contiguous) { isECC = false; isDots = false; }
      }

      if(energy.size())
        outfile << boost::format("%5.4f  %10s") % n->first % energy << endl;
      else {
        outfile << boost::format("#%5.4f  energy not found") % n->first << endl;
        cout << "lct energy can't be found for " << n->first << endl;
      }      

    } // End n
  } // End scope

}

std::string generate_input(const double intern_dist, const std::string name_base)
{
  std::ostringstream stm;
  stm << intern_dist;
  std::string name(name_base + "_" + stm.str() + ".inp");
  std::ofstream file(name.c_str());

  file << "! RHF " +basis + " CCSD(T) TightSCF NoAutoStart" << std::endl;
  file << "" << std::endl;

  file << "*xyz 0 1" << endl;
  // Internuclear distaces
  file << boost::format("N  0.0000000000    0.0000000000    %8.5f") % ( intern_dist/2.0) << std::endl;
  file << boost::format("N  0.0000000000    0.0000000000    %8.5f") % (-intern_dist/2.0) << std::endl;
  file << "*" << endl;
  
  return name;
}
