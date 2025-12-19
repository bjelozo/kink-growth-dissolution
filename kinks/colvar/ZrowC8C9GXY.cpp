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
   
    class ZrowC8C9GXY : public Colvar{

    public:

      vector<AtomNumber> center;

      unsigned int nmol;   // total number of molecules involved in the calculations

      double sigma_t, zeta_t;  // parameters for CV noise reduction step function

      double A_x, nu_x, phi_x1, phi_x2, eta_x;
      double zetal_x, zetau_x, sigmal_x, sigmau_x;
      double A_y, nu_y, phi_y1, phi_y2, eta_y;
      double lb_z, ub_z, lb, ub, mu_z, Var_z, Var2_z;

      double f1_x_i, df1_x_i, f2_x_i, df2_x_i;
      double f_x_i, df_x_i;
      double f1_y_i, df1_y_i, f2_y_i, df2_y_i;
      double f_xyz;

      Vector com_pos;  // point of reference

      double f_cos_i, df_cos_i; 
      void kernel_cos(double, double, double, double);  // kernel function

      double f_lu_i, df_lu_i;
      void kernel_lu(double, double, double, double, double);

      double f_z_i, df_z_i;
      void kernel_z(double, double, double);

      double f_t, df_t;  // filter for reducing noise of CV (simple logistic step function)
      void kernel_t(double);

      vector<double> list;  // inded list of all molecules in defined layer interval
      unsigned int ll;      // length of list
      void layerindices(vector<Vector>&, vector<int>&);  // index list of all molecules positioned within the defined layer interval
      vector<Vector> pos;

      ZrowC8C9GXY(const ActionOptions&);    // CONSTRUCTOR
      ~ZrowC8C9GXY();                       // DESTRUCTOR
      // active methods:
      static void registerKeywords( Keywords& keys );  // KEYWORDS
      virtual void calculate();                        // CALCULATE CV

    };

    PLUMED_REGISTER_ACTION(ZrowC8C9GXY,"ZROWC8C9GXY")
    
    void ZrowC8C9GXY::registerKeywords( Keywords& keys ){

      Colvar::registerKeywords(keys);
      keys.add("atoms", "CENTER", "the labels of the atoms acting as center of the molecules");

      keys.add("compulsory", "ZETA_T", "cutoff for CV noise reduction");
      keys.add("compulsory", "SIGMA_T", "steepness of step function");

      keys.add("compulsory", "NU_X", "number of peaks in x direction");
      keys.add("compulsory", "PHI_X1", "phase shift in x direction");
      keys.add("compulsory", "PHI_X2", "phase shift in x direction");
      keys.add("compulsory", "ETA_X", "cos exponent in x direction");
     
      keys.add("compulsory", "ZETAL_X", "lower position of f(y)");
      keys.add("compulsory", "ZETAU_X", "upper position of f(y)");
      keys.add("compulsory", "SIGMAL_X", "steepness of f(y) on lower side");
      keys.add("compulsory", "SIGMAU_X", "steepness of f(y) on upper side");

      keys.add("compulsory", "NU_Y", "number of peaks in y direction");
      keys.add("compulsory", "PHI_Y1", "phase shift in y direction");
      keys.add("compulsory", "PHI_Y2", "phase shift in y direction");
      keys.add("compulsory", "ETA_Y", "cos exponent in y direction");

      keys.add("compulsory", "LB_Z", "lower bound in z direction");
      keys.add("compulsory", "UB_Z", "upper bound in z direction");
      keys.add("compulsory", "MU_Z", "crystal side position of f(z)");
      keys.add("compulsory", "VAR_Z", "steepness of f(z) on liquid side");

      keys.remove("NOPBC");

    }
 
    ZrowC8C9GXY::ZrowC8C9GXY(const ActionOptions&ao):
      PLUMED_COLVAR_INIT(ao)
    {

      parseAtomList("CENTER", center);

      parse("ZETA_T", zeta_t);
      parse("SIGMA_T", sigma_t);

      parse("NU_X", nu_x);
      parse("PHI_X1", phi_x1);
      parse("PHI_X2", phi_x2);
      parse("ETA_X", eta_x);

      parse("ZETAL_X", zetal_x);
      parse("ZETAU_X", zetau_x);
      parse("SIGMAL_X", sigmal_x);
      parse("SIGMAU_X", sigmau_x);

      parse("NU_Y", nu_y);
      parse("PHI_Y1", phi_y1);
      parse("PHI_Y2", phi_y2);
      parse("ETA_Y", eta_y);

      parse("LB_Z", lb_z);
      parse("UB_Z", ub_z);
      parse("MU_Z", mu_z);
      parse("VAR_Z", Var_z);

      nmol=center.size();
      list.resize(nmol);  // allocate enough space for layer molecules list

      if(nmol==0) error("no molecules specified");

      //cout << "lb: " << lb << endl;
      //cout << "ub: " << ub << endl;
      //cout << "Var_z: " << Var_z << endl;
      //cout << "phi_x1: " << phi_x1 << endl;
      //cout << "phi_x2: " << phi_x2 << endl;
      //cout << "phi_y1: " << phi_y1 << endl;
      //cout << "phi_y2: " << phi_y2 << endl;
      //cout << "nu_x: " << nu_x << endl;
      //cout << "nu_y: " << nu_y << endl;
      //cout << "eta_x: " << eta_x << endl;
      //cout << "eta_y: " << eta_y << endl;

      addValueWithDerivatives();  // informs plumed core that we require space to store the value      
      setNotPeriodic();        // of the CV and that the CV will act on the list of atoms
      requestAtoms(center);    // named center

      checkRead();  // check that everything on the input line has been read properly,
                    // all parse command should follow before checkRead()
    }



void ZrowC8C9GXY::kernel_t(double f_c) {
  // f_c is crystallinity CV. Usually pretty noisy, this function filters out the noise as a simple step function

  f_t = 1/(1 + exp(-sigma_t*(f_c-zeta_t)));
  df_t = sigma_t*f_t*(1-f_t);

}


void ZrowC8C9GXY::kernel_cos(double q_i, double A, double phi, double eta) {
  // q_i is either x or y position of atom

  //double A = nu*M_PI/L;
  double val;
  val = A*(q_i-phi);

  //cout << "val: " << val << endl;
  //cout << "A: " << A << endl;

  f_cos_i = pow(cos(val), eta);
  df_cos_i = -A*sin(val)*eta*pow(cos(val), eta-1);

  //cout << "f_cos_i" << f_cos_i << endl;
  //cout << "df_cos_i" << df_cos_i << endl;
}


void ZrowC8C9GXY::kernel_lu(double q_i, double zetal, double zetau, double sigmal, double sigmau) {

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


void ZrowC8C9GXY::kernel_z(double q_i, double mu_z, double Var2) {

  double val;
  val = q_i-mu_z;

  f_z_i = exp(-(val)*(val)/(2*Var2));
  df_z_i = (-val/Var2)*f_z_i;

}


void ZrowC8C9GXY::layerindices(vector<Vector> &pos, vector<int> &list) {

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
      //cout << "pos[" << k << "][0]" << pos[k][0] << endl;

      k++;

    }
  }

  //comm.Sum(k);
  ll = k;  // length of neighbour list

}


void ZrowC8C9GXY::calculate()
{

  com_pos[0] = getBox()[0][0]/2;
  com_pos[1] = getBox()[1][1]/2;
  com_pos[2] = getBox()[2][2]/2;

  //cout << "com_pos[0]: " << com_pos[0] << endl;
  //cout << "com_pos[1]: " << com_pos[1] << endl;
  //cout << "com_pos[2]: " << com_pos[2] << endl;

  // constants for kernel_z
  Var2_z = Var_z*Var_z;
  lb = com_pos[2] + lb_z;
  ub = com_pos[2] + ub_z;
  //cout << "com_pos[2]: " << com_pos[2] << endl;
  //cout << "lb_z: " << lb_z << endl;
  //cout << "ub_z: " << ub_z << endl;
  //cout << "lb: " << lb << endl;
  //cout << "ub: " << ub << endl;

  // constants for kernel_cos
  A_x = nu_x*M_PI/getBox()[0][0];
  A_y = nu_y*M_PI/getBox()[1][1];

  double cv_val;   // CV
  cv_val = 0;
 
  Tensor virial;   // VIRIAL
  virial.zero();   // no virial contribution
  vector<Vector> deriv(getNumberOfAtoms());  // DERIVATIVES, vector of customized Plumed vectors

  // CV value of each atom for debugging
  vector<double> cv_val_i(nmol);

  vector<Vector> pos(nmol);
  vector<int> list(nmol);
  layerindices(pos, list);  // construct list of molecules iniside the defined interval, pos, together with the list, list, which carries the true index of the molecule.

  unsigned int stride;  // SET THE PARALLELIZATION VARIABLES for the for loops
  unsigned int rank;

  stride=comm.Get_size();  // SET THE PARALLELIZATION VARIABLES for the for loops
  rank=comm.Get_rank();

  for (unsigned int i=rank; i<ll; i+=stride) {  // SUM OVER MOLECULES

    kernel_cos(pos[i][0], A_x, phi_x1, eta_x);
    f1_x_i = f_cos_i;
    df1_x_i = df_cos_i;

    kernel_cos(pos[i][0], A_x, phi_x2, eta_x);
    f2_x_i = f_cos_i;
    df2_x_i = df_cos_i;

    kernel_lu(pos[i][0], zetal_x, zetau_x, sigmal_x, sigmau_x);
    f_x_i = f_lu_i;
    df_x_i = df_lu_i;

    kernel_cos(pos[i][1], A_y, phi_y1, eta_y);
    f1_y_i = f_cos_i;
    df1_y_i = df_cos_i;
    //cout << "f1_y_i: " << f1_y_i << endl;
    //cout << "df1_y_i: " << df1_y_i << endl;
    
    kernel_cos(pos[i][1], A_y, phi_y2, eta_y);
    f2_y_i = f_cos_i;
    df2_y_i = df_cos_i;
    //cout << "f2_y_i: " << f2_y_i << endl;
    //cout << "df2_y_i: " << df2_y_i << endl;

    kernel_z(pos[i][2], mu_z, Var2_z);  // calculates f_z_i and df_z_i
    //f_z_i = f_z_i;
    //df_z_i = df_z_i;

    // sum the total CV
    f_xyz = (f1_x_i*f1_y_i + f2_x_i*f2_y_i)*f_z_i;
    kernel_t(f_xyz*f_x_i);
    cv_val += f_t;

    // for debugging:
    //cv_val_i[list[i]] = f_xyz;
    cv_val_i[list[i]] = f_t;

    // sum the total derivative in each direction
    unsigned int index_i = list[i];  // true index of molecule i
    deriv[index_i][0] += df_t*((df1_x_i*f1_y_i + df2_x_i*f2_y_i)*f_z_i*f_x_i + f_xyz*df_x_i);
    deriv[index_i][1] += df_t*(f1_x_i*df1_y_i + f2_x_i*df2_y_i)*f_z_i*f_x_i;
    deriv[index_i][2] += df_t*(f1_x_i*f1_y_i + f2_x_i*f2_y_i)*df_z_i*f_x_i;

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

    ZrowC8C9GXY::~ZrowC8C9GXY(){
    }

  }
}
