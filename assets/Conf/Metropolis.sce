//==================================================================================
//
//  METROPOLIS ADJUSTED LANGEVIN ALGORITHM
//  ---------------------------------------
//
//  FIGURES 
//  0 = comparison between theoretical and empirical distributions
//  1 = trajectory x^n as a function of the step index n
//  2 = distribution of the increments x^{n+1}-x^n
//  3 = average of the observable under consideration as a function of the iteration step
//
//  WORK TO DO
//  1 -- compute the rejection rate as a function of alpha (for h = 1 and Beta = 2)
//  2 -- what is the optimal value to have the fastest sampling? 
//       For this question, it is worth considering the following observables: 
//       potential energy V(q) and position q
//  3 -- what happens as h is increased? as Beta is decreased?
//  4 -- implement another Metropolis move, for instance MALA but with a noise given by 
//       a uniform distribution between [-1/2, 1/2] rather than a Gaussian of variance 1
//
//  Last modified: Paris, July 3rd, 2017 (Gabriel STOLTZ)
// 
//==================================================================================


//-----------------------
//   POTENTIAL USED
//-----------------------

function y = V(x)
  y = h*((x-1).^2).*(x+1).^2;
endfunction

//-- derivative of the potential --
function y = nablaV(x)
  y = 2*h*( (x-1).*(x+1).^2 + (x-1).^2*(x+1) );
endfunction

//------------------------
//  VARIABLES TO DECLARE
//------------------------

h = 1;
Beta = 2;
Niter = 10000;
alpha = 0.3;
Rplot = 3;
Rmax = 5;
Nhisto = 30;
Ncomp = 100;

//------------------
// INITIALIZATIONS
//-----------------

x = zeros(Niter);
X = x;
x(1) = 1;
V_old = V(x(1));
V_new = V_old; 
E = zeros(Niter);
E(1) = V_new;

//-------------------------------------------
//  GENERATION OF THE STATISTICAL ENSEMBLE 
//-------------------------------------------

rejection = 0;
for i=2:Niter

  V_old = V(x(i-1));
 
  //--------------------
  //   PROPOSITION    
  //--------------------
  noise = rand(1,'normal');
  x_proposed = x(i-1) - nablaV(x(i-1))*alpha + sqrt(2*alpha/Beta)*noise;
  X(i) = x_proposed;
  V_new = V(x_proposed);
  pxy = exp(-0.5*noise^2); //-- this is the probability T(q^n,q^n+1)
  Noise = ( x(i-1) - x_proposed + nablaV(x_proposed)*alpha )/sqrt(2*alpha/Beta);
  pyx = exp(-0.5*Noise^2); //-- this is T(q^n+1,q^n)

  //--------------------------
  //   ACCEPTANCE / REJECTION
  //---------------------------
  p = exp(-Beta*(V_new - V_old))*pyx/pxy;
  if (p >= 1)
    x(i) = x_proposed;
  else 
    if (rand(1,'uniform') < p)
      x(i) = x_proposed;
    else
      x(i) = x(i-1);
      rejection = rejection + 1;
    end
  end

  //-- update of the observables
  E(i) = E(i-1) + x(i); //V(x(i));

end

//--------------------------------------------------------------
//  COMPARISON BETWEEN EMPIRICAL AND THEORETICAL DISTRIBUTIONS 
//--------------------------------------------------------------

//-- un-normalized Boltzmann distribution 
function z = mu(x)
  z = exp(-Beta*V(x));
endfunction

//-- plot of the empirical distribution 
histplot(Nhisto,x);

//-- plot of the theoretical distribution 
C = intg(-Rmax,Rmax,mu);
z = -Rplot:0.05:Rplot;
plot(z,mu(z)/C,'r');

//-- plot of the visited positions 
set("figure_style","new");
f = get("current_figure");
set("current_figure",1);
plot(x);

//-- difference between x^n+1 and x^n
set("figure_style","new");
f = get("current_figure");
set("current_figure",2);
histplot(Ncomp,x(2:Niter)-x(1:Niter-1)); 

//-- convergence of the observable
energie = E' ./ (1:Niter);
set("figure_style","new");
f = get("current_figure");
set("current_figure",4);
plot(energie);

//----------------------------
//     REJECTION RATE
//----------------------------

rejection = rejection/Niter
