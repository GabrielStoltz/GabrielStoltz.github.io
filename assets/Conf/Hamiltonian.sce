//==================================================================================
//
//  NUMERICAL INTEGRATION OF THE HAMILTONIAN DYNAMICS
//  -------------------------------------------------
//
//  FIGURES 
//  0 = phase-space trajectory (momenta as a function of the positions)
//  1 = energy as a function of time
//
//  WORK TO DO
//  1 -- compute trajectories of the system for explicit Euler, symplectic Euler and Heun. Is the energy preserved?
//  2 -- implement the Verlet scheme; check the preservation of energy (order \Delta t^\alpha?); 
//       determine the maximal timestep that can be used 
//  3 -- plot a map of the phase space with the Verlet scheme by playing on initial conditions
//  4 -- [OPTIONAL] implement the implicit midpoint scheme y^{n+1} = y^n + \Delta t * f((y^n+y^{n+1})/2) 
//                  and check its energy conservation 
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
dt = 0.02;
Niter = 10000;
Scheme = 3;  //-- integer giving the scheme to be used: 1 = Euler explicit, 2 = Euler symplectic, 3 = Heun, 4 = Verlet
q0 = -1;     //-- initial position  
p0 = 1;      //-- initial momentum (mass = 1 here)

x = zeros(Niter,2); //-- first colum = positions, second = momenta
x(1,1) = q0;
x(1,2) = p0;
Etot = zeros(Niter);
Epot = zeros(Niter);
Ekin = zeros(Niter);

//--------------------------------
//  INTEGRATION SCHEMES
//--------------------------------

function y = EulerExplicit(x)
  y(1) = x(1) + dt*x(2);
  y(2) = x(2) - dt*nablaV(x(1));
endfunction 

function y = EulerSymplectic(x)
  y(1) = x(1) + dt*x(2);
  y(2) = x(2) - dt*nablaV(y(1));
endfunction 

function y = Heun(x)
  z(1) = x(1) + dt*x(2);
  z(2) = x(2) - dt*nablaV(x(1));
  y(1) = x(1) + 0.5*dt*(x(2)+z(2));
  y(2) = x(2) - 0.5*dt*( nablaV(x(1))+nablaV(z(1)) );
endfunction 

function y = Verlet(x)  // TO COMPLETE
  y(1) = x(1);
  y(2) = x(2);
endfunction 

//--------------------------------
//  INTEGRATION OF THE DYNAMICS
//--------------------------------

for i=2:Niter
  //-- integrate one step with chosen scheme
  select Scheme
    case 1 then
      x(i,:) = EulerExplicit(x(i-1,:))';
    case 2 then
      x(i,:) = EulerSymplectic(x(i-1,:))';
    case 3 then
      x(i,:) = Heun(x(i-1,:))';
    case 4 then
      x(i,:) = Verlet(x(i-1,:))';
  end
  //-- compute the energies
  Ekin(i) = 0.5*x(i,2).^2;
  Epot(i) = V(x(i,1));
  Etot(i) = Ekin(i) + Epot(i);
end

//--------------------------------
//  PLOT OF THE TRAJECTORIES
//--------------------------------

//-- phase space plot: momenta as a function of the positions
set("figure_style","new");
f = get("current_figure");
set("current_figure",0);
plot(x(:,1),x(:,2));

//-- energy as a function of time
set("figure_style","new");
f = get("current_figure");
set("current_figure",1);
plot(1:Niter,Etot','k');
plot(1:Niter,Ekin','r');
plot(1:Niter,Epot','b');