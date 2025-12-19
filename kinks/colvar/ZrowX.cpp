/* +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   Copyright (c) 2014 The plumed team
   (see the PEOPLE file at the root of the distribution for a list of names)

   See http://www.plumed-code.org for more information.

   This file is part of plumed, version 2.

   plumed is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   plumed is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with plumed.  If not, see <http://www.gnu.org/licenses/>.
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ */
#include "Colvar.h"
#include "ActionRegister.h"
#include "tools/KernelFunctions.h"
#include "tools/Communicator.h"

#include <string>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstdio>

using namespace std;

namespace PLMD{
  namespace colvar{
    
//+PLUMEDOC COLVAR FEDERQ_CLUST
/*

*/
//+ENDPLUMEDOC
   
    class ZrowX : public Colvar{

    public:

      vector<AtomNumber> center;

      unsigned int nmol;   // total number of molecules involved in the calculations

      double zetal_x, zetau_x, sigmal_x, sigmau_x;
      double lb_z, ub_z, lb, ub, zetal_z, zetau_z, sigmal_z, sigmau_z;

      double f_x_i, df_x_i;
      double f_z_i, df_z_i;

      Vector com_pos;  // point of reference

      double f_lu_i, df_lu_i;
      void kernel_lu(double, double, double, double, double);

      vector<double> list;  // inded list of all molecules in defined layer interval
      unsigned int ll;      // length of list
      void layerindices(vector<Vector>&, vector<int>&);  // index list of all molecules positioned within the defined layer interval
      vector<Vector> pos;

      ZrowX(const ActionOptions&);    // CONSTRUCTOR
      ~ZrowX();                       // DESTRUCTOR
      // active methods:
      static void registerKeywords( Keywords& keys );  // KEYWORDS
      virtual void calculate();                        // CALCULATE CV

    };

    PLUMED_REGISTER_ACTION(ZrowX,"ZROWX")
    
    void ZrowX::registerKeywords( Keywords& keys ){

      Colvar::registerKeywords(keys);
      keys.add("atoms", "CENTER", "the labels of the atoms acting as center of the molecules");

      keys.add("compulsory", "ZETAL_X", "lower position of f_x_i");
      keys.add("compulsory", "ZETAU_X", "upper position of f_x_i");
      keys.add("compulsory", "SIGMAL_X", "steepness of f_x_i on lower side");
      keys.add("compulsory", "SIGMAU_X", "steepness of f_x_i on upper side");

      keys.add("compulsory", "LB_Z", "lower bound in z direction");
      keys.add("compulsory", "UB_Z", "upper bound in z direction");
      keys.add("compulsory", "ZETAL_Z", "crystal side position of f(z)");
      keys.add("compulsory", "ZETAU_Z", "liquid side position of f(z)");
      keys.add("compulsory", "SIGMAL_Z", "steepness of f(z) on liquid side");
      keys.add("compulsory", "SIGMAU_Z", "steepness of f(z) on crystal side");

      keys.remove("NOPBC");

    }
 
    ZrowX::ZrowX(const ActionOptions&ao):
      PLUMED_COLVAR_INIT(ao)
    {

      parseAtomList("CENTER", center);

      parse("ZETAL_X", zetal_x);
      parse("ZETAU_X", zetau_x);
      parse("SIGMAL_X", sigmal_x);
      parse("SIGMAU_X", sigmau_x);

      parse("LB_Z", lb_z);
      parse("UB_Z", ub_z);
      parse("ZETAL_Z", zetal_z);
      parse("ZETAU_Z", zetau_z);
      parse("SIGMAL_Z", sigmal_z);
      parse("SIGMAU_Z", sigmau_z);

      nmol=center.size();
      list.resize(nmol);  // allocate enough space for layer molecules list

      if(nmol==0) error("no molecules specified");

      //cout << "lb: " << lb << endl;
      //cout << "ub: " << ub << endl;
      //cout << "Var_z: " << Var_z << endl;
      //cout << "phi_x1: " << phi_x1 << endl;
      //cout << "phi_x2: " << phi_x2 << endl;

      addValueWithDerivatives();  // informs plumed core that we require space to store the value      
      setNotPeriodic();        // of the CV and that the CV will act on the list of atoms
      requestAtoms(center);    // named center

      checkRead();  // check that everything on the input line has been read properly,
                    // all parse command should follow before checkRead()
    }


void ZrowX::kernel_lu(double q_i, double zetal, double zetau, double sigmal, double sigmau) {

  double f_l = 0;
  double f_u = 0;

  f_l = 1/(1 + exp(-sigmal*(q_i-zetal)));
  f_u = 1/(1 + exp(-sigmau*(q_i-zetau)));

  if (zetal > zetau) {
    f_lu_i = f_l+(1-f_u);
    df_lu_i = sigmal*f_l*(1-f_l) - sigmau*f_u*(1-f_u);
  } else {
    f_lu_i = f_l*(1-f_u);
    df_lu_i = sigmal*f_l*(1-f_l)*(1-f_u) - sigmau*f_l*f_u*(1-f_u);
  }

  //cout << "q_i: " << q_i << endl;
  //cout << "sigmal: " << sigmal << endl;
  //cout << "sigmau: " << sigmau << endl;
  //cout << "zetal: " << zetal << endl;
  //cout << "zetau: " << zetau << endl;
  //cout << "f_l: " << f_l << endl;
  //cout << "f_u: " << f_u << endl;
  //cout << "f_lu_i: " << f_lu_i << endl;
  //cout << "df_lu_i: " << df_lu_i << endl;

}


void ZrowX::layerindices(vector<Vector> &pos, vector<int> &list) {

  unsigned int k = 0;

  //unsigned int stride;
  //unsigned int rank;

  //stride = comm.Get_size();
  //rank = comm.Get_rank();

  Vector pos_i;
  for (unsigned int i=0; i<nmol; i++) {
  //for (unsigned int i=rank; i<nmol; i+=stride) {

    pos_i = getPosition(i);

    //cout << "pos_i[0]: " << pos_i[0] << ", pos_i[1]: " << pos_i[1] << ", pos_i[2]: " << pos_i[2] << "\n";
    //cout << "lbound_c: " << lbound_c << ", ubound_c: " << ubound_c << ", pos_i[2]: " << pos_i[2] << endl;

    if ( (lb < pos_i[2]) && (pos_i[2] < ub) ) {
      list[k] = i;  // list is the carrier of the true index i for the molecule. k is the index of the molecules in the bias region.

      for (unsigned int ix=0; ix<3; ix++) {
        pos[k][ix] = pos_i[ix] - com_pos[ix];  // relative position to COM atom
      }
      
      if (pos[k][0] < 0) {
        pos[k][0] = pos[k][0] + getBox()[0][0];
      } else if (pos[k][0] > getBox()[0][0]) {
        pos[k][0] = pos[k][0] - getBox()[0][0];
      }
      //cout << "pos[" << k << "][1]" << pos[k][1] << endl;

      k++;

    }
  }

  //comm.Sum(k);
  ll = k;  // length of neighbour list

}


void ZrowX::calculate()
{

  com_pos[0] = getBox()[0][0]/2;
  com_pos[1] = getBox()[1][1]/2;
  com_pos[2] = getBox()[2][2]/2;

  //cout << "com_pos[0]: " << com_pos[0] << endl;
  //cout << "com_pos[1]: " << com_pos[1] << endl;
  //cout << "com_pos[2]: " << com_pos[2] << endl;

  lb = com_pos[2] + lb_z;
  ub = com_pos[2] + ub_z;
  //cout << "com_pos[2]: " << com_pos[2] << endl;
  //cout << "lb_z: " << lb_z << endl;
  //cout << "ub_z: " << ub_z << endl;
  //cout << "lb: " << lb << endl;
  //cout << "ub: " << ub << endl;

  double cv_val;   // CV
  cv_val = 0;
 
  Tensor virial;   // VIRIAL
  virial.zero();   // no virial contribution
  vector<Vector> deriv(getNumberOfAtoms());  // DERIVATIVES, vector of customized Plumed vectors

  // cv_val_i of each atom for debugging
  vector<double> cv_val_i(nmol);

  vector<Vector> pos(nmol);
  vector<int> list(nmol);
  layerindices(pos, list);  // construct list of molecules iniside the defined interval, pos, together with the list, list, which carries the true index of the molecule.

  unsigned int stride;  // SET THE PARALLELIZATION VARIABLES for the for loops
  unsigned int rank;

  stride=comm.Get_size();  // SET THE PARALLELIZATION VARIABLES for the for loops
  rank=comm.Get_rank();

  for (unsigned int i=rank; i<ll; i+=stride) {  // SUM OVER MOLECULES

    kernel_lu(pos[i][0], zetal_x, zetau_x, sigmal_x, sigmau_x);
    f_x_i = f_lu_i;
    df_x_i = df_lu_i;

    kernel_lu(pos[i][2], zetal_z, zetau_z, sigmal_z, sigmau_z);  // calculates f_z_i and df_z_i
    f_z_i = f_lu_i;
    df_z_i = df_lu_i;
    //cout << "pos[i][2]: " << pos[i][2] << endl;
    //cout << "f_z_i: " << f_z_i << endl;
    //cout << "df_z_i: " << df_z_i << endl;

    // sum the total CV
    cv_val += f_x_i*f_z_i;

    // for debugging:
    cv_val_i[list[i]] = f_x_i*f_z_i;

    // sum the total derivative in each direction
    unsigned int index_i = list[i];  // true index of molecule i
    deriv[index_i][0] += df_x_i*f_z_i;
    //deriv[index_i][1] += 0;
    deriv[index_i][2] += f_x_i*df_z_i;

  }

  comm.Sum(deriv);
  comm.Sum(cv_val);
  //comm.Sum(virial);

  for(unsigned i=0;i<nmol;i++) {
    setAtomsDerivatives(i,deriv[i]);
  }

  setBoxDerivatives(virial);
  setValue(cv_val);

  // construct MULTICOLVAR.xyz file
  Vector pos_i;
  cout << nmol << endl;
  cout << " " << getBox()[0][0] << " " << getBox()[1][1] << " " << getBox()[2][2] << endl;
  for (unsigned i=0; i<nmol; i++) {
    pos_i = getPosition(i);
    cout << "X ";
    cout << pos_i[0] << " " << pos_i[1] << " " << pos_i[2] << " ";
    cout << cv_val_i[i] << endl;
  }

}

    ZrowX::~ZrowX(){
    }

  }
}
