#include	<cmath>
#include	<iostream>
#include	<iomanip>
#include	<fstream>
#include	<sstream>
#include	<vector>
#include	<exception>

#include	"global.h"
#include    "filehandler.h"









struct particle fill_particle(const double ra, const double dec,
                              const double  z, const double  wt,
                              const LCDM& lcdm) {
  // Converts from ra, dec and z to 3D Cartesian position.
  struct particle curp;
  if (z<0.0 || z>2.0) {
    std::cout << "z="<<z<<" out of range."<<std::endl;
    myexit(1);
  }
  if (dec<-90.0 || dec>90.0) {
    std::cout << "DEC="<<dec<<" out of range."<<std::endl;
    myexit(1);
  }
  double r    =  lcdm.chiz(z);
  double theta= (90.0 - dec)*M_PI/180.;
  double phi  = (ra        )*M_PI/180.;
  curp.pos[0] = r*sin(theta)*cos(phi);
  curp.pos[1] = r*sin(theta)*sin(phi);
  curp.pos[2] = r*cos(theta);
  curp.wt     = wt;
  return(curp);
}


struct particle fill_particle_box(const double x, const double y,
                              const double  z, const double  wt) {
  struct particle curp;
  curp.pos[0] = x;
  curp.pos[1] = y;
  curp.pos[2] = z;
  curp.wt     = wt;
  return(curp);
}




std::vector<struct particle> read_data(const char fname[], const LCDM& lcdm) {
  // Load the data from a filehandler "file".

  std::vector<long> ndims;
  std::vector<double> ra_vec  = FileHandler::read_dble(fname, "ra", ndims);
  std::vector<double> dec_vec = FileHandler::read_dble(fname, "dec", ndims);
  std::vector<double> z_vec   = FileHandler::read_dble(fname, "z", ndims);

#ifdef  READWEIGHT
  std::vector<double> w_vec   = FileHandler::read_dble(fname, "w", ndims);
#endif

  if ((ra_vec.size() != dec_vec.size())
      || (ra_vec.size() != z_vec.size())
#ifdef READWEIGHT
      || (w_vec.size() != ra_vec.size())
      || (w_vec.size() != dec_vec.size())
      || (w_vec.size() != z_vec.size())
#endif
      || (dec_vec.size() != z_vec.size())) {
      std::cerr<<"Invalid format for file "<<fname
      <<": there are columns with different sizes."<<std::endl;
      std::cerr.flush();
      myexit(1);
  }

  double ra, dec, z, w;
  std::vector<struct particle> P;
  try {P.reserve(10000000);} catch(std::exception& e) {myexception(e);}

  while (!ra_vec.empty()) {
    ra  = ra_vec.back();
    dec = dec_vec.back();
    z   = z_vec.back();

    ra_vec.pop_back();
    dec_vec.pop_back();
    z_vec.pop_back();

#ifdef READWEIGHT
    w   = w_vec.back();
    w_vec.pop_back();
#else
    w   = 1.0;
#endif
    struct particle curp = fill_particle(ra,dec,z,w,lcdm);
    try {
      P.push_back(curp);
    } catch(std::exception& e) {myexception(e);}
  }

  return P;
}


std::vector<struct particle> read_box_data(const char fname[]) {
  // Load the data from a filehandler "file".

  std::vector<long> ndims;
  std::vector<double> x_vec  = FileHandler::read_dble(fname, "x_coord", ndims);
  std::vector<double> y_vec  = FileHandler::read_dble(fname, "y_coord", ndims);
  std::vector<double> z_vec  = FileHandler::read_dble(fname, "z_coord", ndims);

#ifdef  READWEIGHT
  std::vector<double> w_vec   = FileHandler::read_dble(fname, "w", ndims);
#endif

  if ((x_vec.size() != y_vec.size())
      || (x_vec.size() != z_vec.size())
#ifdef READWEIGHT
      || (w_vec.size() != x_vec.size())
      || (w_vec.size() != y_vec.size())
      || (w_vec.size() != z_vec.size())
#endif
      || (y_vec.size() != z_vec.size())) {
      std::cerr<<"Invalid format for file "<<fname
      <<": there are columns with different sizes."<<std::endl;
      std::cerr.flush();
      myexit(1);
  }

  double x, y, z, w;
  std::vector<struct particle> P;
  try {P.reserve(10000000);} catch(std::exception& e) {myexception(e);}

  while (!x_vec.empty()) {
    x   = x_vec.back();
    y   = y_vec.back();
    z   = z_vec.back();

    x_vec.pop_back();
    y_vec.pop_back();
    z_vec.pop_back();

#ifdef READWEIGHT
    w   = w_vec.back();
    w_vec.pop_back();
#else
    w   = 1.0;
#endif
    struct particle curp = fill_particle_box(x,y,z,w);
    try {
      P.push_back(curp);
    } catch(std::exception& e) {myexception(e);}
  }

  return P;
}


void	write_and_destruct_data(std::vector<struct particle>& P, const char fname[], const LCDM& lcdm,
                                bool write_to_box) {
// Writes the normalized positions to a FH "file".
  double x, y, z, ra, dec, redshift, w;

  std::vector<double> x_vec, y_vec, z_vec;
  std::vector<double> ra_vec, dec_vec, redshift_vec;
#ifdef READWEIGHT
  std::vector<double> w_vec;
#endif

  try {
    if (write_to_box) {
      ra_vec.reserve(10000000);
      dec_vec.reserve(10000000);
      redshift_vec.reserve(10000000);
    } else {
      x_vec.reserve(10000000);
      y_vec.reserve(10000000);
      z_vec.reserve(10000000);
    }
#ifdef READWEIGHT
    w_vec.reserve(10000000);
#endif
  } catch(std::exception& e) {myexception(e);}

  while (!P.empty()) {
    struct particle p = P.back();
    P.pop_back();

    x = box.ctr[0]+box.L*(p.pos[0]-0.5);
    y = box.ctr[1]+box.L*(p.pos[1]-0.5);
    z = box.ctr[2]+box.L*(p.pos[2]-0.5);
    w = p.wt;

    ra       = 180. / M_PI * atan2(y, x);
    dec      = 90. - 180. / M_PI * atan2(sqrt(x*x + y*y), z);
    redshift = lcdm.zchi(sqrt(x*x + y*y + z*z));

    try {
      ra_vec.push_back(ra);
      dec_vec.push_back(dec);
      redshift_vec.push_back(redshift);

      if (write_to_box) {
        x_vec.push_back(x);
        y_vec.push_back(y);
        z_vec.push_back(z);
      } else {
        ra_vec.push_back(ra);
        dec_vec.push_back(dec);
        redshift_vec.push_back(redshift);
      }

#ifdef READWEIGHT
      w_vec.push_back(w);
#endif
    } catch(std::exception& e) {myexception(e);}
  }


  if (!write_to_box) {
    std::vector<long int> ndims = {(long) ra_vec.size()};
    FileHandler::write_dble(fname, "ra", ra_vec, ndims);
    FileHandler::write_dble(fname, "dec", dec_vec, ndims);
    FileHandler::write_dble(fname, "z", redshift_vec, ndims);
  } else {
    std::vector<long int> ndims = {(long) x_vec.size()};
    FileHandler::write_dble(fname, "x_coord", x_vec, ndims);
    FileHandler::write_dble(fname, "y_coord", y_vec, ndims);
    FileHandler::write_dble(fname, "z_coord", z_vec, ndims);
  }

#ifdef READWEIGHT
  FileHandler::write_dble(fname, "w", w_vec, ndims);
#endif

}

