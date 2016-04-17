#ifndef HH_Parameters_HH
#define HH_Parameters_HH
#include <iosfwd>
#include <string>
struct parameters
{
  //! max number of iteration for Gauss-Siedel
  int   itermax;
  //! Tolerance for stopping criterion
  double  toler;
  //! Bar length
   double L;
  //! First longitudinal dimension
  double a1;
 //! Second longitudinal dimension
  double a2;
  //! Dirichlet condition
  double To;
  //! External temperature 
  double Te;
  //! Conductivity
  double k;
  //! Convection coefficient
  double hc;
  //! Number of elements
  int M;
  //! Result file
  std::string Resfile;
  //! Norm Type
  int Norm; //1 for Rn, 2 for L², 3 for H1
  //! Algo choose
  int algo;
  //! Constructor takes default values
  parameters():
    itermax(1000000),
    toler(1e-8),
    L(40.),
    a1(4.),
    a2(50.),
    To(46.),
    Te(20.),
    k(0.164),
    hc(1.e-6*200.),
    M(100),
    Resfile("resultfile.dat"),
    Norm(1),
    algo(1)
  {}
 
};


//! Prints parameters
std::ostream & operator << (std::ostream &,const parameters &);
#endif
