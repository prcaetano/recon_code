#ifndef	_FILEHANDLER_H_
#define _FILEHANDLER_H_

#include	<cmath>
#include	<cstdlib>
#include	<vector>
#include	<string>
#include	<sstream>
#include	<fstream>
#include	<iostream>
#include	<iomanip>
#include	<exception>
#include	<sys/stat.h>	// For the POSIX "stat" and "mkdir" functions.

#ifndef _GLOBAL_H
#include "global.h"         // For the myexit and myexception functions
#endif


/*

A simple and lightweight system for mocking up files with 1D arrays (vectors)
of values stored by keyword.
Overloaded read/write methods do a simple version of multi-D arrays, stored
as 1D vectors but with additional information on the dimensions.
Assumes "myexit" and "myexception" are defined elsewhere.  Samples would be

inline	void	myexit(const int iflag) {
  exit(iflag);
}
inline void	myexception(std::exception& e) {
  std::cerr<<"Exception: "<<e.what()<<std::endl;
  myexit(1);
}

*/


//extern void	myexit(const int iflag);
//extern void	myexception(std::exception& e);    Templates already defined on global


class FileHandler {
// Code to implement a simple file structure which can store "keywords"
// which are arrays of values.  The code makes each "file" a directory
// on the file system with files for each "keyword" stored in the
// appropriate directory.  Files have a type appended to their name as
// one level of error checking upon re-read.
// Note that keywords are case-sensitive, and you can overload fields
// of different variables types (e.g. a long 'key1' and a float 'key1').
// To keep this stateless there is no ability to read part of a field or
// to append to a field.
// Parallel routines should divide the data across multiple directories,
// each accessed independently, to reduce the total memory per read/write
// to an acceptable value and to take maximum advantage of parallel I/O.
//
// Author:	Martin White	(UCB)
private:
  static void make_dir(const char fname[]) {
    // Checks to see whether the file (actually a directory) exists already,
    // if not it makes it, else returns.
    struct stat statbuf;
    int ret = stat(fname,&statbuf);
    if (ret==-1) { // File doesn't exist, we need to make it.
      // Just choose some semi-random permissions here...
      ret = mkdir(fname,S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      if (ret==-1) {
        std::cerr<<"Unable to make file "<<fname<<std::endl;
        myexit(1);
      }
    }
  }
public:
  static std::vector<int> read_int(const char fname[], const char field[]) {
    // Read an integer field.
    std::stringstream ss;
    ss << fname << "/" << field << "nd.i4";
    std::ifstream ifs(ss.str().c_str(),std::ios::binary);
    if (!ifs) {
      std::cerr<<"Unable to find int field "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    long nsize;
    ifs.read((char *)&nsize,sizeof(long));
    if (ifs.fail()) {
      std::cerr<<"Unable to read length for "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    std::vector<int> ret;
    try {ret.resize(nsize);} catch(std::exception& e) {myexception(e);}
    ifs.read((char *)&ret[0],nsize*sizeof(int));
    if (ifs.fail()) {
      std::cerr<<"Unable to read data for "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ifs.close();
    return(ret);
  }
  static std::vector<int> read_int(const char fname[], const char field[],
                                   std::vector<long>& ndims) {
    // This overloaded read_int will read a multidimensional array,
    // returning a vector of ints but also the dimensions so you can
    // interpret it as a multidimensional array.
    std::stringstream ss;
    ss << fname << "/" << field << ".nd.i4";
    std::ifstream ifs(ss.str().c_str(),std::ios::binary);
    if (!ifs) {
      std::cerr<<"Unable to find int field "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    int  ndim;
    ifs.read((char *)&ndim,sizeof(int));
    if (ifs.fail()) {
      std::cerr<<"Unable to read dims for "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    try{ndims.resize(ndim);}catch(std::exception& e) {myexception(e);}
    ifs.read((char *)&ndims[0],ndim*sizeof(long));
    if (ifs.fail()) {
      std::cerr<<"Unable to read ndims for "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    long nsize=1;
    for (int i=0; i<ndims.size(); ++i) nsize *= ndims[i];
    std::vector<int> ret;
    try {ret.resize(nsize);} catch(std::exception& e) {myexception(e);}
    ifs.read((char *)&ret[0],nsize*sizeof(int));
    if (ifs.fail()) {
      std::cerr<<"Unable to read data for "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ifs.close();
    return(ret);
  }
  static std::vector<float> read_float(const char fname[], const char field[]){
    // Read a float field.
    std::stringstream ss;
    ss << fname << "/" << field << "nd.f4";
    std::ifstream ifs(ss.str().c_str(),std::ios::binary);
    if (!ifs) {
      std::cerr<<"Unable to find float field "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    long nsize;
    ifs.read((char *)&nsize,sizeof(long));
    if (ifs.fail()) {
      std::cerr<<"Unable to read length for "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    std::vector<float> ret;
    try {ret.resize(nsize);} catch(std::exception& e) {myexception(e);}
    ifs.read((char *)&ret[0],nsize*sizeof(float));
    if (ifs.fail()) {
      std::cerr<<"Unable to read data for "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ifs.close();
    return(ret);
  }
  static std::vector<float> read_float(const char fname[], const char field[],
                                       std::vector<long>& ndims){
    // This overloaded read_float will read a multidimensional array,
    // returning a vector of floats but also the dimensions so you can
    // interpret it as a multidimensional array.
    std::stringstream ss;
    ss << fname << "/" << field << ".nd.f4";
    std::ifstream ifs(ss.str().c_str(),std::ios::binary);
    if (!ifs) {
      std::cerr<<"Unable to find float field "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    int ndim;
    ifs.read((char *)&ndim,sizeof(int));
    if (ifs.fail()) {
      std::cerr<<"Unable to read dims for "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    try{ndims.resize(ndim);}catch(std::exception& e) {myexception(e);}
    ifs.read((char *)&ndims[0],ndim*sizeof(long));
    if (ifs.fail()) {
      std::cerr<<"Unable to read ndims for "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    long nsize=1;
    for (int i=0; i<ndims.size(); ++i) nsize *= ndims[i];
    std::vector<float> ret;
    try {ret.resize(nsize);} catch(std::exception& e) {myexception(e);}
    ifs.read((char *)&ret[0],nsize*sizeof(float));
    if (ifs.fail()) {
      std::cerr<<"Unable to read data for "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ifs.close();
    return(ret);
  }
  static std::vector<long> read_long(const char fname[], const char field[]){
    // Read a long field.
    std::stringstream ss;
    ss << fname << "/" << field << "nd.i8";
    std::ifstream ifs(ss.str().c_str(),std::ios::binary);
    if (!ifs) {
      std::cerr<<"Unable to find long field "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    long nsize;
    ifs.read((char *)&nsize,sizeof(long));
    if (ifs.fail()) {
      std::cerr<<"Unable to read length for "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    std::vector<long> ret;
    try {ret.resize(nsize);} catch(std::exception& e) {myexception(e);}
    ifs.read((char *)&ret[0],nsize*sizeof(long));
    if (ifs.fail()) {
      std::cerr<<"Unable to read data for "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ifs.close();
    return(ret);
  }
  static std::vector<long> read_long(const char fname[], const char field[],
                                     std::vector<long>& ndims){
    // This overloaded read_long will read a multidimensional array,
    // returning a vector of longs but also the dimensions so you can
    // interpret it as a multidimensional array.
    std::stringstream ss;
    ss << fname << "/" << field << ".nd.i8";
    std::ifstream ifs(ss.str().c_str(),std::ios::binary);
    if (!ifs) {
      std::cerr<<"Unable to find long field "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    int ndim;
    ifs.read((char *)&ndim,sizeof(int));
    if (ifs.fail()) {
      std::cerr<<"Unable to read dims for "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    try{ndims.resize(ndim);}catch(std::exception& e) {myexception(e);}
    ifs.read((char *)&ndims[0],ndim*sizeof(long));
    if (ifs.fail()) {
      std::cerr<<"Unable to read ndims for "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    long nsize=1;
    for (int i=0; i<ndims.size(); ++i) nsize *= ndims[i];
    std::vector<long> ret;
    try {ret.resize(nsize);} catch(std::exception& e) {myexception(e);}
    ifs.read((char *)&ret[0],nsize*sizeof(long));
    if (ifs.fail()) {
      std::cerr<<"Unable to read data for "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ifs.close();
    return(ret);
  }
  static std::vector<double> read_dble(const char fname[], const char field[]){
    // Read a double field.
    std::stringstream ss;
    ss << fname << "/" << field << ".nd.f8";
    std::ifstream ifs(ss.str().c_str(),std::ios::binary);
    if (!ifs) {
      std::cerr<<"Unable to find dble field "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    long nsize;
    ifs.read((char *)&nsize,sizeof(long));
    if (ifs.fail()) {
      std::cerr<<"Unable to read length for "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    std::cout<<"Nsize: "<<nsize<<std::endl;
    std::vector<double> ret;
    try {ret.resize(nsize);} catch(std::exception& e) {myexception(e);}
    ifs.read((char *)&ret[0],nsize*sizeof(double));
    if (ifs.fail()) {
      std::cerr<<"Unable to read data for "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ifs.close();
    return(ret);
  }
  static std::vector<double> read_dble(const char fname[], const char field[],
                                       std::vector<long>& ndims){
    // This overloaded read_dble will read a multidimensional array,
    // returning a vector of doubles but also the dimensions so you can
    // interpret it as a multidimensional array.
    std::stringstream ss;
    ss << fname << "/" << field << ".nd.f8";
    std::ifstream ifs(ss.str().c_str(),std::ios::binary);
    if (!ifs) {
      std::cerr<<"Unable to find dble field "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    int ndim;
    ifs.read((char *)&ndim,sizeof(int));
    if (ifs.fail()) {
      std::cerr<<"Unable to read dims for "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    try{ndims.resize(ndim);}catch(std::exception& e) {myexception(e);}
    ifs.read((char *)&ndims[0],ndim*sizeof(long));
    if (ifs.fail()) {
      std::cerr<<"Unable to read ndims for "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    long nsize=1;
    for (int i=0; i<ndims.size(); ++i) nsize *= ndims[i];
    std::vector<double> ret;
    try {ret.resize(nsize);} catch(std::exception& e) {myexception(e);}
    ifs.read((char *)&ret[0],nsize*sizeof(double));
    if (ifs.fail()) {
      std::cerr<<"Unable to read data for "<<field<<" in "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ifs.close();
    return(ret);
  }
  static void write_int(const char fname[], const char field[],
                        const std::vector<int>& val) {
    // Writes an int field.
    make_dir(fname);
    std::stringstream ss;
    ss << fname << "/" << field << "nd.i4";
    std::ofstream ofs(ss.str().c_str(),std::ios::binary|std::ios::trunc);
    if (!ofs) {
      std::cerr<<ss.str().c_str()<<std::endl;
      std::cerr<<"Unable to write field "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    long nsize=val.size();
    ofs.write((char *)&nsize,sizeof(long));
    if (ofs.fail()) {
      std::cerr<<"Unable to write length for "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ofs.write((char *)&val[0],nsize*sizeof(int));
    if (ofs.fail()) {
      std::cerr<<"Unable to write data for "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ofs.close();
  }
  static void write_int(const char fname[], const char field[],
                        const std::vector<int>& val,
                        const std::vector<long>& ndims) {
    // This overloaded write_int will write a multidimensional array.
    make_dir(fname);
    std::stringstream ss;
    ss << fname << "/" << field << ".nd.i4";
    std::ofstream ofs(ss.str().c_str(),std::ios::binary|std::ios::trunc);
    if (!ofs) {
      std::cerr<<ss.str().c_str()<<std::endl;
      std::cerr<<"Unable to write field "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    int  ndim =ndims.size();
    ofs.write((char *)&ndim,sizeof(int));
    if (ofs.fail()) {
      std::cerr<<"Unable to write ndim for "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ofs.write((char *)&ndims[0],ndim*sizeof(long));
    if (ofs.fail()) {
      std::cerr<<"Unable to write ndims for "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ofs.write((char *)&val[0],val.size()*sizeof(int));
    if (ofs.fail()) {
      std::cerr<<"Unable to write data for "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ofs.close();
  }
  static void write_float(const char fname[], const char field[],
                          const std::vector<float>& val) {
    make_dir(fname);
    std::stringstream ss;
    ss << fname << "/" << field << "nd.f4";
    std::ofstream ofs(ss.str().c_str(),std::ios::binary|std::ios::trunc);
    if (!ofs) {
      std::cerr<<"Unable to write field "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    long nsize=val.size();
    ofs.write((char *)&nsize,sizeof(long));
    if (ofs.fail()) {
      std::cerr<<"Unable to write length for "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ofs.write((char *)&val[0],nsize*sizeof(float));
    if (ofs.fail()) {
      std::cerr<<"Unable to write data for "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ofs.close();
  }
  static void write_float(const char fname[], const char field[],
                          const std::vector<float>& val,
                          const std::vector<long>& ndims) {
    // This overloaded write_float will write a multidimensional array.
    make_dir(fname);
    std::stringstream ss;
    ss << fname << "/" << field << ".nd.f4";
    std::ofstream ofs(ss.str().c_str(),std::ios::binary|std::ios::trunc);
    if (!ofs) {
      std::cerr<<"Unable to write field "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    int  ndim =ndims.size();
    ofs.write((char *)&ndim,sizeof(int));
    if (ofs.fail()) {
      std::cerr<<"Unable to write ndim for "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ofs.write((char *)&ndims[0],ndim*sizeof(long));
    if (ofs.fail()) {
      std::cerr<<"Unable to write ndims for "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ofs.write((char *)&val[0],val.size()*sizeof(float));
    if (ofs.fail()) {
      std::cerr<<"Unable to write data for "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ofs.close();
  }
  static void write_long(const char fname[], const char field[],
                         const std::vector<long>& val) {
    make_dir(fname);
    std::stringstream ss;
    ss << fname << "/" << field << "nd.i8";
    std::ofstream ofs(ss.str().c_str(),std::ios::binary|std::ios::trunc);
    if (!ofs) {
      std::cerr<<"Unable to write field "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    long nsize=val.size();
    ofs.write((char *)&nsize,sizeof(long));
    if (ofs.fail()) {
      std::cerr<<"Unable to write length for "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ofs.write((char *)&val[0],nsize*sizeof(long));
    if (ofs.fail()) {
      std::cerr<<"Unable to write data for "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ofs.close();
  }
  static void write_long(const char fname[], const char field[],
                         const std::vector<long>& val,
                         const std::vector<long>& ndims) {
    make_dir(fname);
    std::stringstream ss;
    ss << fname << "/" << field << ".nd.i8";
    std::ofstream ofs(ss.str().c_str(),std::ios::binary|std::ios::trunc);
    if (!ofs) {
      std::cerr<<"Unable to write field "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    int  ndim =ndims.size();
    ofs.write((char *)&ndim,sizeof(int));
    if (ofs.fail()) {
      std::cerr<<"Unable to write ndim for "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ofs.write((char *)&ndims[0],ndim*sizeof(long));
    if (ofs.fail()) {
      std::cerr<<"Unable to write ndims for "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ofs.write((char *)&val[0],val.size()*sizeof(long));
    if (ofs.fail()) {
      std::cerr<<"Unable to write data for "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ofs.close();
  }
  static void write_dble(const char fname[], const char field[],
                         const std::vector<double>& val) {
    make_dir(fname);
    std::stringstream ss;
    ss << fname << "/" << field << "nd.f8";
    std::ofstream ofs(ss.str().c_str(),std::ios::binary|std::ios::trunc);
    if (!ofs) {
      std::cerr<<"Unable to write field "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    long nsize=val.size();
    ofs.write((char *)&nsize,sizeof(long));
    if (ofs.fail()) {
      std::cerr<<"Unable to write length for "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ofs.write((char *)&val[0],nsize*sizeof(double));
    if (ofs.fail()) {
      std::cerr<<"Unable to write data for "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ofs.close();
  }
  static void write_dble(const char fname[], const char field[],
                         const std::vector<double>& val,
                         const std::vector<long>& ndims) {
    make_dir(fname);
    std::stringstream ss;
    ss << fname << "/" << field << ".nd.f8";
    std::ofstream ofs(ss.str().c_str(),std::ios::binary|std::ios::trunc);
    if (!ofs) {
      std::cerr<<"Unable to write field "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    int  ndim =ndims.size();
    ofs.write((char *)&ndim,sizeof(int));
    if (ofs.fail()) {
      std::cerr<<"Unable to write ndim for "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ofs.write((char *)&ndims[0],ndim*sizeof(long));
    if (ofs.fail()) {
      std::cerr<<"Unable to write ndims for "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ofs.write((char *)&val[0],val.size()*sizeof(double));
    if (ofs.fail()) {
      std::cerr<<"Unable to write data for "<<field<<" to "<<fname<<std::endl;
      std::cerr.flush();
      myexit(1);
    }
    ofs.close();
  }
};
#endif
