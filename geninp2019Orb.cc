
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
    std::string scf_name   ("./scf_r_" + stm.str() + "_.out");
    std::string casscf_name("./casscf_r_" + stm.str() + "_.out");
    std::string mrci_name  ("./lct_r_" + stm.str() + "_.out");    
    std::vector<std::string> file_names;
    file_names.push_back(scf_name);
    file_names.push_back(casscf_name);
    file_names.push_back(mrci_name);
    names.insert(len_files::value_type(*l, file_names));

    std::string name(generate_input(*l, n_base));
    std::string scf_exec   ("~/ORZSAVE3/orz/obj-opt/src/sci/scf/scf ./" + name + " | tee " + scf_name);
    std::string casscf_exec("~/ORZSAVE3/orz/obj-opt/src/sci/casscf/casscf ./" + name + " | tee " + casscf_name);
    std::string mrci_exec  ("~/ORZSAVE3/orz/obj-opt/src/sci/lct/lct ./" + name + " | tee " + mrci_name);
    system(scf_exec.c_str());
    system(casscf_exec.c_str());
    system(mrci_exec.c_str());
  }

  {
    
    std::string fileName("scf_summary");
    {
      std::ostringstream pno;
      pno << boost::format("%e") % _tCutPNO;
      std::ostringstream pairs;
      pairs << boost::format("%e") % _tCutPairs;
      fileName += (_doCEPT2 ? "_CEPT2_" : "_CASPT2_") + pno.str() + "_" + pairs.str() + "_" + basis + ".out";
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

      bool istotal(false);
      bool isscf(false);
      bool isenergy(false);
      bool contiguous(false);
      while(!infile.eof()){
        infile >> line;
        size_t where;
        if(yes) { energy = line; break; }
        else if(!istotal && !isscf && !isenergy){
          where = line.find("total", 0);
          if(where != std::string::npos)
            { istotal = true; contiguous = true; }
          else contiguous = false;
        }
        else if(istotal && !isscf && !isenergy){
          where = line.find("scf", 0);
          if(where != std::string::npos)
            { isscf = true; contiguous = true; }
          else contiguous = false;
        }
        else if(istotal && isscf && !isenergy){
          where = line.find("energy", 0);
          if(where != std::string::npos)
            { isenergy = true; contiguous = true; }
          else contiguous = false;
        }
        else contiguous = false;

        if(isscf && istotal && isenergy) { yes = true; /*outfile << where << " llll " << line << endl;*/ }
        else yes = false;
        if(!contiguous) { isscf = false; istotal = false; isenergy = false; }
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

    std::string fileName("casscf_summary");
    {
      std::ostringstream pno;
      pno << boost::format("%e") % _tCutPNO;
      std::ostringstream pairs;
      pairs << boost::format("%e") % _tCutPairs;
      fileName += (_doCEPT2 ? "_CEPT2_" : "_CASPT2_") + pno.str() + "_" + pairs.str() + "_" + basis + ".out";
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

      bool isCASSCF(false);
      bool isenergy(false);
      bool iseq(false);
      bool contiguous(false);
      while(!infile.eof()){
        infile >> line;
        size_t where;
        if(yes) { energy = line; break; }
        else if(!isCASSCF && !isenergy && !iseq){
          where = line.find("CASSCF", 0);
          if(where != std::string::npos)
            { isCASSCF = true; contiguous = true; }
          else contiguous = false;
        }
        else if(isCASSCF && !isenergy && !iseq){
          where = line.find("energy", 0);
          if(where != std::string::npos)
            { isenergy = true; contiguous = true; }
          else contiguous = false;
        }
        else if(isCASSCF && isenergy && !iseq){
          where = line.find("=", 0);
          if(where != std::string::npos)
            { iseq = true; contiguous = true; }
          else contiguous = false;
        }
        else contiguous = false;

        if(isCASSCF && isenergy && iseq) { yes = true; /*outfile << where << " llll " << line << endl;*/ }
        else yes = false;
        if(!contiguous) { isCASSCF = false; isenergy = false; iseq = false; }
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

    std::string fileName("lct_summary");
    {
      std::ostringstream pno;
      pno << boost::format("%e") % _tCutPNO;
      std::ostringstream pairs;
      pairs << boost::format("%e") % _tCutPairs;
      fileName += (_doCEPT2 ? "_CEPT2_" : "_CASPT2_") + pno.str() + "_" + pairs.str() + "_" + basis + ".out";
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
      
      bool isFINAL(false);
      bool isTOTAL(false);
      bool isENERGY(false);
      bool isDots(false);
      bool contiguous(false);
      while(!infile.eof()){
        infile >> line;
        size_t where;
        if(yes) { energy = line; break; }
        else if(!isFINAL && !isTOTAL && !isENERGY && !isDots){
          where = line.find("FINAL", 0);
          if(where != std::string::npos)
            { isFINAL = true; contiguous = true; }
          else contiguous = false;
        }
        else if(isFINAL && !isTOTAL && !isENERGY && !isDots){
          where = line.find("TOTAL", 0);
          if(where != std::string::npos)
            { isTOTAL = true; contiguous = true; }
          else contiguous = false;
        }
        else if(isFINAL &&  isTOTAL && !isENERGY && !isDots){
          where = line.find("ENERGY", 0);
          if(where != std::string::npos)
            { isENERGY = true; contiguous = true; }
          else contiguous = false;
        }
        else if(isFINAL &&  isTOTAL && isENERGY && !isDots){
          where = line.find("...", 0);
          if(where != std::string::npos)
            { isDots = true; contiguous = true; }
          else contiguous = false;
        }        
        else contiguous = false;

        if(isFINAL && isTOTAL && isENERGY && isDots) { yes = true; /*outfile << where << " llll " << line << endl;*/ }
        else yes = false;
        if(!contiguous) { isFINAL = false; isTOTAL = false; isENERGY = false; isDots = false; }
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
  std::string name(name_base + "_" + stm.str() + ".py");
  std::ofstream file(name.c_str());

  file << "from sci import *                                                           " << std::endl;
  file << "                                                                            " << std::endl;
  file << "qatom.basis = baslib(\"" + basis + "\")                                             " << std::endl;
  file << "qatom.unit = unit.angstrom                                                  " << std::endl;
  file << "                                                                            " << std::endl;

  // Internuclear distaces
  file << boost::format("N2 = [ qatom(atom.N, [    .0000000000,    .0000000000,    %8.5f     ]  ),     ") % (intern_dist/2.0) << std::endl;
  file << boost::format("       qatom(atom.N, [    .0000000000,    .0000000000,   -%8.5f     ]  ), ]   ") % (intern_dist/2.0) << std::endl;

  file << "                                                                            " << std::endl;
  file << "qmolecule = molecule( N2, charge=0, mult=1, sym=sym.C1 )                   " << std::endl;
  file << "                        " << std::endl;
  file << "nstate=1                        " << std::endl;
  file << "scf.dodirect = False                                                        " << std::endl;
  file << "#scf.calc_ao_tei = False                                                      " << std::endl;
  file << "scf.show_orbitals = True                                                    " << std::endl;
  file << "qatom.basis = baslib(\"def2-svp/jkfit\")" << std::endl;
  file << "#qatom.basis = baslib(\"aug-cc-pv5z-ri/jkfit\")" << std::endl;
  file << "" << std::endl;
  file << "acene = [ i.new_basis() for i in N2]" << std::endl;
  file << "qmolecule_RI_JK = molecule(acene, charge=0, mult=1, sym=sym.C1)" << std::endl;
  file << "casscf.no4c_RI_JK = True" << std::endl;
  file << "hint.use_RI_JK = True" << std::endl;
  file << "hint.itrf_RI_JK_memory = 2000" << std::endl;
  file << "ga_memory = 100000" << std::endl;  
  file << "                                                                            " << std::endl;
  file << "casscf.scfinp = key.scfread                                                 " << std::endl;
  file << "#casscf.orbconfig = \"symlist\"                                                " << std::endl;
  file << "casscf.orbconfig = \"orblist\"                                                " << std::endl;
  file << "" << std::endl;
  file << "#casscf.frozen = [0]" << std::endl;
  file << "#casscf.closed = [4]" << std::endl;
  file << "#casscf.occ    = [6]" << std::endl;
  file << "casscf.frozen = []" << std::endl;
  if     ( (intern_dist>=0.95) && (intern_dist< 1.3) )
    file << "casscf.closed = [0,1,2,3]" << std::endl;    
    file << "casscf.occ    = [4,5,6,7,8,10]" << std::endl;
  else if( (intern_dist>=1.3 ) && (intern_dist< 1.8) )
    file << "casscf.closed = [0,1,2,3]" << std::endl;
    file << "casscf.occ    = [4,5,6,7,8,9]" << std::endl;
  else if( (intern_dist>=1.8 ) && (intern_dist< 2.0) ){
    file << "casscf.closed = [0,1,2,4]" << std::endl;
    file << "casscf.occ    = [3,5,6,7,8,9]" << std::endl;
  }
  else if( (intern_dist>=2.0 ) && (intern_dist< 2.8) ){
    file << "casscf.closed = [0,1,2,3]" << std::endl;
    file << "casscf.occ    = [4,5,6,7,8,9]" << std::endl;
  }
///  else if( (intern_dist>=2.8 ) && (intern_dist< ) ){
///    file << "casscf.closed = [0,1,2,3]" << std::endl;
///    file << "casscf.occ    = [4,5,6,7,8,9]" << std::endl;
///  }
  
  file << "" << std::endl;
  file << "#casscf.isymstate = 0" << std::endl;
  file << "#casscf.natural_orbitals = True                                              " << std::endl;
  file << "casscf.maxMacroloop = 100                 " << std::endl;
  file << "casscf.maxStepSize  = 1.0                                                   " << std::endl;
  file << "casscf.tolMacroloop = 5e-6                  " << std::endl;
  file << "casscf.tolMicroloop = 1e-8                                                  " << std::endl;
  file << "casscf.use_check_pt = False                                                 " << std::endl;
  file << "casscf.citype = \"fullci\"                                                    " << std::endl;
  file << "casscf.roots = {1:1}                                                        " << std::endl;
  file << "                                                                            " << std::endl;
  file << "casscf.natural_orbitals = True  " << std::endl;
  file << "casscf.maxMacroloop = 100       " << std::endl;
  file << "casscf.maxStepSize  = 1.0       " << std::endl;
  file << "casscf.tolMacroloop = 5e-6      " << std::endl;
  file << "casscf.tolMicroloop = 1e-8      " << std::endl;
  file << "casscf.use_check_pt = False     " << std::endl;
  file << "casscf.citype = \"fullci\"      " << std::endl;
  file << "casscf.roots = {1:1}            " << std::endl;
  file << "fullci.rdm_newcode = True       " << std::endl;
  file << "" << std::endl;
  file << "icmr.citype = casscf.citype" << std::endl;
  file << "" << std::endl;

  file << "# ICMR" << std::endl;
  file << "icmr.skip_itrf  = False" << std::endl;
  file << "icmr.skip_casci = False" << std::endl;
  file << "icmr.scfinp = casscf.scfinp" << std::endl;
  file << "icmr.localize_type = casscf.localize_type" << std::endl;
  file << "icmr.orbconfig = casscf.orbconfig" << std::endl;
  file << "" << std::endl;
  file << "#icmr.frozen = [ 0]" << std::endl;
  file << "#icmr.closed = [ 4]" << std::endl;
  file << "#icmr.occ    = [ 6]" << std::endl;
  file << "icmr.frozen = []" << std::endl;
  file << "icmr.closed = [0,1,2,3]" << std::endl;
  file << "icmr.occ    = [4,5,6,7,8,9]" << std::endl;    
  file << "" << std::endl;
  file << "icmr.roots  = {1:nstate}" << std::endl;
  file << "icmr.citype = casscf.citype" << std::endl;
  file << "icmr.use_3RDM = True" << std::endl;
  file << "icmr.use_4RDM = True" << std::endl;
  file << "icmr.use_d4cum_of = False" << std::endl;
  file << "icmr.ampAllocType = 2 #a 1=disk, 2=GA" << std::endl;
  file << "icmr.icb_type = -2 # ICBTYPE # 0=MS-MR, 1=MS-SR, 2=SS-SR, -1=MS-SR*, -2=SS-SR*" << std::endl;
  file << "icmr.xms_type =  0 # XMSTYPE # 0=Non-XMS, 1=XMS" << std::endl;
  file << "icmr.method = \"caspt2\"" << std::endl;
  file << "icmr.pt.tolOrthExternal_1rdm = 1e-8" << std::endl;
  file << "icmr.pt.tolOrthExternal      = 1e-8" << std::endl;
  file << "icmr.pt.tolOrthInternal      = 1e-8" << std::endl;
  file << "icmr.pt.ipeashift = 0.0" << std::endl;
  file << "icmr.pt.lvl_shift = -1.0" << std::endl;
  file << "icmr.pt.toliter   = 1e-20" << std::endl;
  file << "" << std::endl;
  file << "lct.DoPNO = False" << std::endl;
  file << "" << std::endl;
  file << "# LCT" << std::endl;
  file << "lct.test        = False" << std::endl;
  file << "lct.DoPNO       = True" << std::endl;
  file << "lct.UseOrthPNOs = True" << std::endl;
  file << "lct.DoLoc       = True" << std::endl;
  file << "lct.DoCEPT2     = " << (_doCEPT2 ? "True" : "False") << std::endl;
  file << "" << std::endl;
  file << "lct.TCutPNOrho0 = " << boost::format("%14.10e") % _tCutPNO   << std::endl;
  file << "lct.TCutPNOrho1 = " << boost::format("%14.10e") % _tCutPNO   << std::endl;
  file << "lct.TCutPNOrho2 = " << boost::format("%14.10e") % _tCutPNO   << std::endl;
  file << "lct.TCutPairs   = " << boost::format("%14.10e") % _tCutPairs << std::endl;
  file << "lct.MaxIter = 500" << std::endl << std::endl;
  
  return name;
}
