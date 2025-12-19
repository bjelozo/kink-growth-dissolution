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
   
    class ZsphereG2 : public Colvar{

    public:

      vector<AtomNumber> center;

      unsigned int nmol;   // total number of molecules involved in the calculations

      double x, y, z;  // location of the center of sphere (seen from center of simulation box)
      double d_max;   // max distance along z axis
      double sigma;   // width of gaussian

      double lb_z, ub_z;    // limits of z-axis for list
      //Vector com_pos;  // center of simulation box
      Vector pos_ref;    // center of sphere (seen from com_pos)

      // ramp sphere functions
      double sigma2;
      double gauss;
      double dgauss;
      Vector dd;
      void densy(double, double, double);

      // relevant atoms list

      vector<double> list;  // inded list of all molecules in defined layer interval
      unsigned int ll;      // length of list
      void layerindices(vector<Vector>&, vector<int>&);  // index list of all molecules positioned within the defined layer interval
      vector<Vector> pos;
      ZsphereG2(const ActionOptions&);    // CONSTRUCTOR
      ~ZsphereG2();                       // DESTRUCTOR
      // active methods:
      static void registerKeywords( Keywords& keys );  // KEYWORDS
      virtual void calculate();                        // CALCULATE CV

    };

    PLUMED_REGISTER_ACTION(ZsphereG2,"ZSPHEREG2")
 
    void ZsphereG2::registerKeywords( Keywords& keys ){

      Colvar::registerKeywords(keys);
      keys.add("atoms", "CENTER", "the labels of the atoms acting as center of the molecules");
      keys.add("compulsory", "X", "x position of sphere center");
      keys.add("compulsory", "Y", "y position of sphere center");
      keys.add("compulsory", "Z", "z position of sphere center");
      keys.add("compulsory", "SIGMA", "width of gaussian");
      keys.add("compulsory", "D_MAX", "maximum distance along z");

      keys.remove("NOPBC");

    }
 
    ZsphereG2::ZsphereG2(const ActionOptions&ao):
      PLUMED_COLVAR_INIT(ao)
    {

      parseAtomList("CENTER", center);
      parse("X", x);
      parse("Y", y);
      parse("Z", z);
      parse("SIGMA", sigma);
      parse("D_MAX", d_max);

      nmol=center.size();
      list.resize(nmol);  // allocate enough space for layer molecules list

      if(nmol==0) error("no molecules specified");

      addValueWithDerivatives();  // informs plumed core that we require space to store the value      
      setNotPeriodic();        // of the CV and that the CV will act on the list of atoms
      requestAtoms(center);    // named center

      checkRead();  // check that everything on the input line has been read properly,
                    // all parse command should follow before checkRead()
    }


void ZsphereG2::densy(double x_i, double y_i, double z_i) {
 
  double d, d2;

  d2 = x_i*x_i + y_i*y_i + z_i*z_i;
  d = sqrt(d2);
  //cout << "d: " << d << endl;
  //cout << "x_i: " << x_i << "  ";
  //cout << "y_i: " << y_i << "  ";
  //cout << "z_i: " << z_i << endl;

  if (d > 0) {
    dd[0] = x_i/d;
    dd[1] = y_i/d;
    dd[2] = z_i/d;
  } else {
    dd[0] = 0;
    dd[1] = 0;
    dd[2] = 0;
  }
  //cout << "dd[0]: " << dd[0] << endl;
  //cout << "dd[1]: " << dd[1] << endl;
  //cout << "dd[2]: " << dd[2] << endl;

  gauss = exp(-0.5*d2/sigma2);
  dgauss = -(d/sigma2)*gauss;
  //cout << "d2: " << d2 << " ";
  //cout << "d: "  << d << " ";
  //cout << "gauss: " << gauss << endl;
  //cout << "dgauss: " << dgauss << endl;

}


void ZsphereG2::layerindices(vector<Vector> &pos, vector<int> &list) {

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

    if ( (lb_z < pos_i[2]) && (pos_i[2] < ub_z) ) {
      list[k] = i;  // list is the carrier of the true index i for the molecule. k is the index of the molecules in the bias region.
 
      if (pos_i[0] < pos_ref[0]-getBox()[0][0]/2) {
        pos_i[0] = pos_i[0] + getBox()[0][0];
      } else if (pos_i[0] > pos_ref[0]+getBox()[0][0]/2) {
        pos_i[0] = pos_i[0] - getBox()[0][0];
      }
      pos[k][0] = pos_i[0] - pos_ref[0];
      //cout << "pos[" << k << "][0]: " << pos[k][0] << "  ";

      if (pos_i[1] < pos_ref[1]-getBox()[1][1]/2) {
        pos_i[1] = pos_i[1] + getBox()[1][1];
      } else if (pos_i[1] > pos_ref[1]+getBox()[1][1]/2) {
        pos_i[1] = pos_i[1] - getBox()[1][1];
      }
      pos[k][1] = pos_i[1] - pos_ref[1];
      //cout << "pos[" << k << "][1]: " << pos[k][1] << "  ";

      pos[k][2] = pos_i[2] - pos_ref[2];
      //cout << "pos[" << k << "][2]: " << pos[k][2] << endl;

      k++;

    }
  }

  //comm.Sum(k);
  ll = k;  // length of neighbour list

}


void ZsphereG2::calculate()
{

  // position of sphere within simulation box
  pos_ref[0] = x;
  pos_ref[1] = y;
  pos_ref[2] = z + getBox()[2][2]/2;
  //cout << "pos_ref[0]: " << pos_ref[0] << "  ";
  //cout << "pos_ref[1]: " << pos_ref[1] << "  ";
  //cout << "pos_ref[2]: " << pos_ref[2] << endl;
  sigma2 = sigma*sigma;
  lb_z = pos_ref[2] - d_max;
  ub_z = pos_ref[2] + d_max;
  //cout << "sigma: " << sigma << endl;
  //cout << "lb_z: " << lb_z << "  ";
  //cout << "ub_z: " << ub_z << endl;

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

    densy(pos[i][0], pos[i][1], pos[i][2]);
    
    // sum the total CV
    cv_val += gauss;

    // for debugging:
    cv_val_i[list[i]] = gauss;

    // sum the total derivative in each direction
    unsigned int index_i = list[i];  // true index of molecule i
    deriv[index_i][0] += dgauss * dd[0];
    deriv[index_i][1] += dgauss * dd[1];
    deriv[index_i][2] += dgauss * dd[2];

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

    ZsphereG2::~ZsphereG2(){
    }

  }
}
