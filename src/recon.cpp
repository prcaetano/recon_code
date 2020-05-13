#include	<cmath>
#include	<cstdlib>
#include	<iostream>
#include	<iomanip>
#include	<fstream>
#include	<sstream>
#include	<vector>
#include	<exception>

#include	"global.h"
#include	"multigrid.h"

// A code to perform reconstruction, using the lowest order algorithm.
//
// Author:	Martin White	(UCB/LBNL)
// Written:	20-Apr-2015
// Modified:	22-Apr-2015	(Take two random catalogs)
// Modified:    06-May-2020, Pedro Rangel (Optionally uses box data, takes output fnames and Om from cli)
//
//
// This code uses a standard multigrid technique to compute the
// displacements given an observed density field and then move the
// objects and "randoms" back along the displacement vector.
//
// Want to get rid of continual allocation and deallocation of arrays.
// I think we just need two arrays at each level, could call them A & B,
// or V1 and V2.
//
// Could make this use double precision throughout.  Memory isn't a
// major problem so far.
//


// Global variables.
float	bias=1,beta=0;
struct	Box box;



void    myexit(const int flag)
{
  std::cout.flush();
  std::cerr.flush();
  exit(flag);
}


void    myexception(const std::exception& e)
{
  std::cout<<": Exception: "<<e.what()<<std::endl;
  std::cout.flush();
  std::cerr.flush();
  exit(1);
}






extern "C" int	recon(char *data_file, char *randoms1_file, char *randoms2_file,
                      char *output_file, char *shifted_randoms_file,
                      float b, float f, float Rf, float Om, bool is_sim_data)
{
  bias = b;             // Sets global variable
  beta = f/b;           // Sets global variable
  LCDM lcdm(Om);

#ifdef	TESTMG
  // Make a cosine wave (only in x-direction) and solve for it with
  // beta=0.  Used to test the MG solver convergence on a simple problem.
  box.L=1.25;  box.ctr[0]=box.ctr[1]=box.ctr[2]=0.6;
  const float dt=2*M_PI/Ng;
  std::vector<float> src(Ng*Ng*Ng);
  for (int ix=0; ix<Ng; ++ix)
    for (int iy=0; iy<Ng; ++iy)
      for (int iz=0; iz<Ng; ++iz)
        src[Ng*Ng*ix+Ng*iy+iz] = cos(ix*dt);
  bias = 1; beta = 0;
  std::vector<float> ans = MultiGrid::fmg(src,Ng);
  for (int ix=0; ix<Ng; ++ix)
    std::cout<<std::fixed<<std::setprecision(6)
             <<ix<<" "<<cos(ix*dt)/(2*M_PI*2*M_PI)
             <<" "<<ans[Ng*Ng*ix+Ng*0   +Ng/4]
             <<" "<<ans[Ng*Ng*ix+Ng*Ng/2+0]
             <<" "<<ans[Ng*Ng*ix+Ng*0   +Ng/2]<<std::endl;
  return(0);
#endif

  // Read the data and figure out the 3D positions and enclosing box.
  std::vector<struct particle> D;
  std::vector<struct particle> R1;
  std::vector<struct particle> R2;
  if (is_sim_data) {
    D = read_box_data(data_file);
    R1= read_box_data(randoms1_file);
    R2= read_box_data(randoms2_file);
  }
  else {
    D = read_data(data_file,lcdm);
    R1= read_data(randoms1_file,lcdm);
    R2= read_data(randoms2_file,lcdm);
  }
  std::cout<<"# Read "<<std::setw(10)<<D.size()
           <<" objects from "<<data_file<<std::endl;
  std::cout<<"# Read "<<std::setw(10)<<R1.size()
           <<" randoms from "<<randoms1_file<<std::endl;
  std::cout<<"# Read "<<std::setw(10)<<R2.size()
           <<" randoms from "<<randoms2_file<<std::endl;
  remap_pos(D,R1,R2);
  std::cout<<"# Enclosing survey in a box of side "<<box.L<<" Mpc/h."
           <<std::endl;
  std::cout<<"# Grid/mesh size is "<<box.L/Ng<<" Mpc/h"
           <<" and filter scale is "<<Rf<<" Mpc/h."
           <<std::endl;

#ifndef	SKIPRAW
  write_data(D ,"data_raw.xyzw");
  write_data(R2,"rand_raw.xyzw");
#endif
  
  // Make the density (contrast) grid and solve for the
  // displacement potential, phi.
  std::vector<float> delta = make_grid(D,R1,Rf);
  std::vector<float> phi   = MultiGrid::fmg(delta,Ng);

  // Shift the particles and randoms back -- if you want to not enhance
  // the line-of-sight shift for the randoms you need to change beta before
  // calling shift_obj.
  shift_obj(D ,phi);
  shift_obj(R2,phi);

  write_data(D ,output_file);
  write_data(R2,shifted_randoms_file);

  return(0);
}

int	main(int argc, char **argv)
{
  if (argc!=11 && argc!=10) {
    std::cout<<"Usage: recon <data-file> <random-file> <random-file>"
             <<" <output-file> <output-shifted-randoms-file>"
             <<" <bias> <f-growth> <R-filter> <Omega m> [is-simulation]"<<std::endl;
    myexit(1);
  }
  bool is_sim_data = false;
  if (argc==11) {
    is_sim_data = bool(atoi(argv[10]));
  }
#ifdef DISTANT_OBSERVER_ZAXIS
  std::cout<<"# Distant observer in z-axis assumed for RSD corrections."<<std::endl;
#else
  std::cout<<"# Fixed observer at (x,y,z) = (0,0,0) assumed for RSD corrections."<<std::endl;
#endif
  return recon(argv[1], argv[2], argv[3],
               argv[4], argv[5],
               atof(argv[6]), atof(argv[7]), atof(argv[8]), atof(argv[9]), is_sim_data);
}

