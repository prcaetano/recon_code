#ifndef __LCDM_H_
#define __LCDM_H_

#include	<cmath>
#include	<iostream>
#include	<iomanip>
#include	<fstream>
#include	<sstream>
#include	<vector>
#include	<exception>


#include	"global.h"


// This class is provded as header-only since it is so small.
//
// Author:	Martin White	(UCB/LBNL)






class   LCDM {
// A class which allows rapid distance calculations, and little else.
// Assumes "classical" LCDM, i.e. radiation and/or massive neutrino
// contributions are not included.
private:
  static const int      N=4*1024;
  std::vector<double>   zz,chi;
  std::vector<double>   zz2, chi2;
  double                omm,oml,fact,chifact,dz,dchi;
  double H(const double z) {
    const double Lhub=2997.925; // Mpc/h.
    double HofZ=sqrt(omm*pow(1+z,3)+oml);
    return(HofZ/Lhub);
  }
public:
  LCDM(const double omegam) {
    // Initialize the distance tables to be looked up by our main method.
    omm = omegam;
    oml = 1.0-omegam;
    const double zmax=10;
    try {
      zz.resize(N);
      chi.resize(N);
    } catch(std::exception& e) {myexception(e);}
    zz[0]=chi[0]=0;
    if (omegam<0) {	// This is code for using a linear Hubble law.
      for (int i=1; i<N; ++i) {
        zz[i] = i*zmax/(N-1.0);
        chi[i]= 2997.925*zz[i];
      }
    } else {		// Do a "real" calculation.
      for (int i=1; i<N; ++i) {
        zz[i] = i*zmax/(N-1.0);
        const int Nint=20;        // Must be even.
        int    wt = 4;
        double hh = (zz[i]-zz[i-1])/Nint;
        double sum= 1.0/H(zz[i-1])+1.0/H(zz[i]);
        for (int j=1; j<Nint; j++) {
          double z=zz[i-1]+j*hh;
          sum += 1.0/H(z) * wt;
          wt   = 8/wt;
        }
        sum *= hh/3.0;
        chi[i] = chi[i-1]+sum;
      }
    }
    fact = N/zmax;
    dz   = zz[1]-zz[0];

    try {
      zz2.resize(N);
      chi2.resize(N);
    } catch(std::exception& e) {myexception(e);}

    int idx=0;
    zz2[0]=chi2[0]=0;
    const double chimax=chi[N-10];
    for (int i=1; i<N; ++i) {
        chi2[i] = i*chimax/(N-1.);
        while(chi[idx+1] < chi2[i])
            idx++;
        zz2[i] = zz[idx] + dz * (chi2[i] - chi[idx]) / (chi[idx+1] - chi[idx]);
    }
    chifact = N/chimax;
    dchi    = chi2[1] - chi2[0];
  }
  double chiz(const double z) const {
    // Returns the comoving angular diameter distance, in Mpc/h.  Uses
    // a very quick linear interpolation.
    int ibin=(int)(fact*z);
    if (ibin<0 || ibin>=zz.size()-1) {
      std::cout<<"chiz: z="<<z<<" out of lookup range ["<<zz[0]<<","
               <<zz[zz.size()-1]<<")."<<std::endl;
      myexit(1);
    }
    // Note this needs one "extra" bin beyond z to do the interpolation.
    double chival=chi[ibin]+(chi[ibin+1]-chi[ibin])*(z-zz[ibin])/dz;
    return(chival);
  }
  double zchi(const double chival) const {
    // Returns the redshift given a comoving angular diameter distance, in Mpc/h.
    int ibin=(int)(chifact*chival);
    if (ibin<0 || ibin>=chi2.size()-1) {
      std::cout<<"zchi: chival="<<chival<<" out of lookup range ["<<chi2[0]<<","
               <<chi2[chi2.size()-1]<<")."<<std::endl;
      myexit(1);
    }
    // Note this needs one "extra" bin beyond z to do the interpolation.
    double zval=zz2[ibin]+(zz2[ibin+1]-zz2[ibin])*(chival-chi2[ibin])/dchi;
    return(zval);
  }
};


#endif
