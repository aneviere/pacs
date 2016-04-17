#include <iostream> // input output
#include <cmath> // (for sqrt)
#include <vector>
#include <tuple>
#include "readParameters.hpp"
#include "GetPot.hpp"
#include "gnuplot-iostream.hpp"// interface with gnuplot
#include <string>
/*!
  @file main.cpp
  @brief Temperature distribution in a 1D bar.

  @detail
    We solve  \f$ -T^{\prime\prime}(x)+act*(T(x)-T_e)=0, 0<x<L \f$ with 
    boundary conditions \f$ T(0)=To; T^\prime(L)=0\f$
    
    **************************************************
    Linear finite elements
    Iterative resolution by Gauss Siedel.
    **************************************************
    
    Example adapted by Luca Formaggia from  a code found in 
    "Simulation numerique an C++" di I. Danaila, F. Hecht e
    O. Pironneau.
*/
//! helper function
void printHelp()
{
  std::cout<<"USAGE: main [-h] [-v] -p parameterFile (default: parameters.pot)"<<std::endl;
  std::cout<<"-h this help"<<std::endl;
  std::cout<<"-v verbose output"<<std::endl;
}

//! main program
int main(int argc, char** argv)
{
  using namespace std; // avoid std::
  int status(0); // final program status
  GetPot   cl(argc, argv);
  if( cl.search(2, "-h", "--help") )
    {
      printHelp();
      return 0;
    }
  // check if we want verbosity
  bool verbose=cl.search(1,"-v");
  // Get file with parameter values
  string filename = cl.follow("parameters.pot","-p");
  cout<<"Reading parameters from "<<filename<<std::endl;
  // read parameters
  const parameters param=readParameters(filename,verbose);
  // Transfer parameters to local variables
  // I use references to save memory (not really an issue here, it is just
  // to show a possible  use of references)
  const int&    itermax= param.itermax;   //max number of iteration for Gauss-Siedel
  const double& toler=param.toler;   // Tolerance for stopping criterion
  // Here I use auto (remember that you need const and & if you want constant references)
  const auto& L= param.L;  // Bar length
  const auto& a1=param.a1; // First longitudinal dimension
  const auto& a2=param.a2; //  Second longitudinal dimension
  const auto& To=param.To; // Dirichlet condition
  const auto& Te=param.Te; // External temperature (Centigrades)
  const auto& k=param.k;  // Thermal conductivity
  const auto& hc=param.hc; // Convection coefficient
  const auto& M=param.M; // Number of grid elements
  const auto& file=param.Resfile; //Name of the result file
  const auto& norm=param.Norm; //Type of norm 1 for Rn, 2 for L², 3 for H1
  const auto& algo=param.algo; // GaussSeidel or Thomas algorithm

  
  //! Precomputed coefficient for adimensional form of equation
  const auto act=2.*(a1+a2)*hc*L*L/(k*a1*a2);

  // mesh size
  const auto h=1./M;
  
  // Solution vector
  std::vector<double> theta(M+1);
  
  // Gauss Siedel is initialised with a linear variation
  // of T
  
  for( int m=0;m <= M;++m)
     theta[m]=(1.-m*h)*(To-Te)/Te;
  
  // Gauss-Seidel
  // epsilon=||x^{k+1}-x^{k}||
  // Stopping criteria epsilon<=toler
  if (algo==1){
     cout<<"Gauss Seidel algorithm"<<endl;
 
  int iter=0;
  double xnew, epsilon, xold, thetaold, deriv;
     do
       { epsilon=0.;
         if (norm==3){
         	deriv=0.; //used for h1 norm
         	xold=0.; //used for h1 norm
	 	thetaold=0.; //used for h1 norm
	 }

	 // first M-1 row of linear system
         for(int m=1;m < M;m++)
         {   
	   xnew  = (theta[m-1]+theta[m+1])/(2.+h*h*act);
	   epsilon += (xnew-theta[m])*(xnew-theta[m]);
	   
         if (m>1 && norm==3) {
		  deriv+=(xnew-theta[m]-(xold-thetaold))*(xnew-theta[m]-(xold-thetaold))/h;//used for h1 norm
		  thetaold = theta[m];//used for h1 norm
		  xold=xnew; //used for h1 norm
		}
           theta[m] = xnew;
         }

	 //Last row
	 xnew = theta[M-1]; 
	 epsilon += (xnew-theta[M])*(xnew-theta[M]); //epsilon in Rn
         if (norm==2) epsilon = h*epsilon; //epsilon in L²;
         if (norm==3){ 
		deriv+=(xnew-theta[M]-(xold-thetaold))*(xnew-theta[M]-(xold-thetaold))/h; //used for h1 norm
         	epsilon = h*epsilon+deriv; //epsilon in H1
		}
	 
	 theta[M]=  xnew; 

	 iter=iter+1;     
       }while((sqrt(epsilon) > toler) && (iter < itermax) );

    if(iter<itermax)
      cout << "M="<<M<<"  Convergence in "<<iter<<" iterations"<<endl;
    else
      {
	cerr << "NOT CONVERGING in "<<itermax<<" iterations "<<
	  "||dx||="<<sqrt(epsilon)<<endl;
	status=1;
      }
   }


  // Thomas Algorithm

 
   if (algo==2) {
   cout<<"Thomas algorithm"<<endl;
 
  //Creation of a, b, c, d vectors
  std::vector<double> a(M);
  std::vector<double> b(M+1);
  std::vector<double> c(M);
  std::vector<double> d(M+1);
  
  
  // Initialization
  
  for(int m=0;m < M;++m) {

     a[m]=-1;
     c[m]=-1;
     b[m]=2+h*h*2.*act;
     d[m]=0;
    }
  b[M]=1;
  d[0]=To;
  d[M]=0;
    
  int iter=0;
  double xnew, epsilon, xold, thetaold, deriv;
     do
       { epsilon=0.;
         if (norm==3){
         	deriv=0.; //used for h1 norm
         	xold=0.; //used for h1 norm
	 	thetaold=0.; //used for h1 norm
	 }

	 //Modification of c,d
         c[0]=c[0]/b[0];
         d[0]=d[0]/b[0];
         for (int m=1; m<M;m++)
	 {
		c[m]=c[m]/(b[m]-a[m]*c[m-1]);
		d[m]=(d[m]-a[m]*d[m-1])/(b[m]-a[m]*c[m-1]);
	 }
         d[M]=(d[M]-a[M-1]*d[M-1])/(b[M]-a[M-1]*c[M-1]);

	 // New theta
	 //Last row
	 xnew = d[M]; 
	 epsilon += (xnew-theta[M])*(xnew-theta[M]); //epsilon in Rn	
	 xold=xnew; 
         thetaold = theta[M];
	 theta[M]=  xnew; 
         
         // the M-1 orthers rows
         for(int m=(M-1);m >=0;m--)
         {   
	   xnew  = d[m]-c[m]*theta[m+1];
	   epsilon += (xnew-theta[m])*(xnew-theta[m]);
	   
         if (norm==3) {
		  deriv+=(xold-thetaold-(xold-theta[m]))*(xnew-theta[m]-(xold-thetaold))/h;//used for h1 norm
		  thetaold = theta[m];//used for h1 norm
		  xold=xnew; //used for h1 norm
		}
           theta[m] = xnew;
         }

	 if (norm==2) epsilon = h*epsilon; //epsilon in L²;
	 if(norm==3) epsilon = h*epsilon+deriv; //epsilon in H1


	 iter=iter+1;     
       }while((sqrt(epsilon) > toler) && (iter < itermax) );

    if(iter<itermax)
      cout << "M="<<M<<"  Convergence in "<<iter<<" iterations"<<endl;
    else
      {
	cerr << "NOT CONVERGING in "<<itermax<<" iterations "<<
	  "||dx||="<<sqrt(epsilon)<<endl;
	status=1;
      }
   }

 // Analitic solution

    vector<double> thetaa(M+1);
     for(int m=0;m <= M;m++)
       thetaa[m]=Te+(To-Te)*cosh(sqrt(act)*(1-m*h))/cosh(sqrt(act));

     // writing results with format
     // x_i u_h(x_i) u(x_i) and lauch gnuplot 

     Gnuplot gp;
     std::vector<double> coor(M+1);
     std::vector<double> sol(M+1);
     std::vector<double> exact(M+1);


     //cout<<"Result file: resultfile.dat"<<endl;
     //ofstream f("resultfile.dat");
     cout << "Result file "<< file<< endl;
     ofstream f(file);

     for(int m = 0; m<= M; m++)
       {
	 // \t writes a tab 
         f<<m*h*L<<"\t"<<Te*(1.+theta[m])<<"\t"<<thetaa[m]<<endl;
	 // An example of use of tie and tuples!
         
	 std::tie(coor[m],sol[m],exact[m])=
	   std::make_tuple(m*h*L,Te*(1.+theta[m]),thetaa[m]);
       }
     // Using temporary files (another nice use of tie)
     gp<<"plot"<<gp.file1d(std::tie(coor,sol))<<
       "w lp title 'uh',"<< gp.file1d(std::tie(coor,exact))<<
       "w l title 'uex'"<<std::endl;
     f.close();
     return status;
}
