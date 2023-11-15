#include <TMB.hpp>

//schaefer function return biomass year 2+
template <class Type> Type schaefer(Type b1, Type r, Type K, Type m, Type u_i){
    Type u_use = exp(u_i)/(1.0+exp(u_i)); //transform exploitation rate 
    Type b2 = b1 + b1*r*(1.0-pow(b1/K,m-1.0)) - u_use*b1;//calculate biomass next year 
    return b2;
}

template<class Type>
Type objective_function<Type>::operator() ()
{
  
  // we have 3 vectors catches index years
  /* DATA */
  DATA_VECTOR(catches);
  DATA_VECTOR(index);
  DATA_IVECTOR(index_years);
  int nyrs = catches.size();
  int nobs = index.size();
  
  /* PARAMETERS */
  PARAMETER(logK); //carrying capacity
  PARAMETER(logr); //growth rate 
  PARAMETER(logm); //shape param
  PARAMETER(logSigmaC); //observation error catch
  PARAMETER(logSigmaI); //observation error index
  //PARAMETER(logq);
  PARAMETER_VECTOR(upar); //exploitation rate for each catch year 

  Type neglogL = 0.; //objective function value negative log likelihood 
  
  //first goal estimate biomass time series deterministic means there is not contribution to the likelihood
  //biomass year +1 
  
  Type m=1.0+exp(logm); // including the estimation of m shape param
  
  /* POPULATION DYNAMICS */
  vector<Type> biomass(nyrs+1);
  biomass(0) = exp(logK); //biomass first year = K 
  for (int iyr=1;iyr<=nyrs;iyr++) { //loop over years to calcutate biomass vector, now compute the predictions for our data Y hat expected value for the catch
    //schaefer biomass update goes here!
    //Type u_use = exp(upar(iyr-1))/(1.0+exp(upar(iyr-1)));
    biomass(iyr) = schaefer(biomass(iyr-1),exp(logr),exp(logK),m,upar(iyr-1));
    //biomass(iyr) = biomass(iyr-1)*(1.0+exp(logr)*(1.0-pow(biomass(iyr-1)/exp(logK),m-1.0)) - u_use);

    //std::cout << iyr << " " << biomass(iyr) << std::endl;
  }
    
  //model predictions
  
  // first predict surveys we have to predict q is a portion of our biomass (q*biomass) 

  vector<Type> index_pred(nobs);//length of our vector nobs, this goes through our observations
  //estimate q analytically
  Type logq = 0.;
  vector<Type> predbio(nobs);
  for (int iobs=0; iobs<nobs; iobs++) {
    predbio(iobs) = 0.5*(biomass(index_years(iobs))+biomass(index_years(iobs) + 1)); //index equation on the pdf Iy hat
    logq += log(index(iobs)/predbio(iobs));//Rather than treat catchability q as an estimated model parameter, we can derive and use the analytical MLE for q:   
  }
  logq = logq/nobs; //
  index_pred = exp(logq)*predbio; //*exp(logq); predict q and then apply that q to our biomass
  
  // then predict catches
  vector<Type> catch_pred(nyrs);// this one goes through our years 
  for (int iyr=0;iyr<nyrs;iyr++)
    catch_pred(iyr) = biomass(iyr)*exp(upar(iyr))/(1.0+exp(upar(iyr)));
  //std::cout << catch_pred << std::endl;
  
  // Likelihood for our observations 
   //CATCHES
   neglogL -= sum(dnorm(log(catches),log(catch_pred),exp(logSigmaC), true));
      SIMULATE {
     catches = exp(rnorm(log(catch_pred), exp(logSigmaC)));//randomly generating normal distrb values using predicted values and std dev
     REPORT(catches);
   }

   //SURVEY
   neglogL -= sum(dnorm(log(index),log(index_pred),exp(logSigmaI), true));//randomly generating normal distrb values using predicted values and std dev
     SIMULATE {
     index = exp(rnorm(log(index_pred), exp(logSigmaI)));
     REPORT(index);
   }

  ADREPORT(logq);
  ADREPORT(biomass);
  ADREPORT(m);
  ADREPORT(neglogL);
  
  return neglogL;

  /* Quick Reference
     ===============

     ** Macros to read data and declare parameters:

     _Template_Syntax_              _C++_type_                     _R_type_
     DATA_VECTOR(name)              vector<Type>                   vector
     DATA_MATRIX(name)              matrix<Type>                   matrix
     DATA_SCALAR(name)              Type                           numeric(1)
     DATA_INTEGER(name)             int                            integer(1)
     DATA_FACTOR(name)              vector<int>                    factor
     DATA_SPARSE_MATRIX(name)       Eigen::SparseMatrix<Type>      dgTMatrix
     DATA_ARRAY(name)               array<Type>                    array
     PARAMETER_MATRIX(name)         matrix<Type>                   matrix
     PARAMETER_VECTOR(name)         vector<Type>                   vector
     PARAMETER_ARRAY(name)          array<Type>                    array
     PARAMETER(name)                Type                           numeric(1)

     ** Macro to report intermediate expressions back to R:

     REPORT(x)
     ADREPORT(x)

     ** Basic constructors:

     vector<Type> v(n1);
     matrix<Type> m(n1,n2);
     array<Type> a(n1,n2,n3)

     ** Basic operations:

     v+v,v-v,v*v,v/v                Pointwise binary operations
     m*v                            Matrix-vector multiply
     a.col(i)                       R equivalent of a[,,i]
     a.col(i).col(j)                R equivalent of a[,j,i]
     a(i,j,k)                       R equivalent of a[i,j,k]
     exp(v)                         Pointwise math
     m(i,j)                         R equivalent of m[i,j]
     v.sum()                        R equivalent of sum(v)
     m.transpose()                  R equivalent of t(m)

     ** Distributions:

     Type dnbinom2(const Type &x, const Type &mu, const Type &var, int give_log=0)
     Type dpois(const Type &x, const Type &lambda, int give_log=0)
     Type dlgamma(Type y, Type shape, Type scale, int give_log=0)
     Type dnorm(Type x, Type mean, Type sd, int give_log=0)

     ** Parallel accumulator declaration (only methods "+=" and "-="):
     
     parallel_accumulator<Type> res(this);

  */

}
