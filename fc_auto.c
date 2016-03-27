#include "auto_f2c.h"
#include <math.h>

#define F2C -1

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/*   The nondimensional Yuval  model - ODE version (no diffusion terms)   */
/* ---------------------------------------------------------------------- */
/*   b_t = lamb*w*b*pow((1+eta*b),2)*(1-b) - b                            */
/*   w_t = p - nu*w*(1-rho*b) - lamb*w*b*pow((1+eta*b),2)                 */
/* ---------------------------------------------------------------------- */
int func (integer ndim, const doublereal *u, const integer *icp,
          const doublereal *par, integer ijac,
          doublereal *f, doublereal *dfdu, doublereal *dfdp)
{
  /* System generated locals */
  //integer dfdu_dim1, dfdp_dim1;
  
  /* Local variables */
  doublereal b,w, X;
  doublereal lamb,eta,p,nu,rho;
  doublereal dummy_b,dummy_w;
  doublereal pi, num_periods,L;
  pi = atan(1.) * 4.0;  
  /* Defining the parameters */
 
  lamb    = par[1+F2C];
  eta     = par[2+F2C];
  p       = par[3+F2C];
  nu      = par[4+F2C];
  rho     = par[5+F2C];
  

  num_periods = par[17+F2C];
  L       = par[18+F2C];

  // 2*pi/L = kf/2 for locking. We multiply by num_periods to make the domain larger
  
  /* Function Body */
  b      = u[0];
  w      = u[1];
 

  // f[i] = u[i]_x
  // We multiply by L because the derivatives are relative to "AUTO"s space.
  // x_real = [0,L], and x_auto = [0,1]. Therefore x_real = L * x_auto
  // d/dx_auto = L * d/dx_real
  f[0] = (lamb*w*b*pow((1+eta*b),2)*(1-b) - b);
  f[1] = (p - nu*w*(1-rho*b) - lamb*w*b*pow((1+eta*b),2));

  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int stpnt (integer ndim, doublereal x,
           doublereal *u, doublereal *par)
{
  /* local variables */
  doublereal B,W,X;
  doublereal K,E,M,N,Lambda,Gamma,R,PP;
  doublereal lamb,eta,p,nu,rho;
  doublereal dummy_b,dummy_w;
  doublereal L, amp, pi, num_periods, dm, dp,offset,tet,batch,tat;
  /*  defining the numerical value of the parameters */
  K = 0.4;
  E = 7.0;
  M = 10.5;
  N = 15;
  Lambda = 0.9;
  Gamma  = 12;
  R  = 0.7;
  PP = 20;

  lamb    = (K*Gamma)/M;
  eta     = E*K;
  p       = (Lambda*PP)/(K*Gamma*M);
  nu      = N/M;
  rho     = R;
  // With delta = 7.5 and a = 0.5 we have eps_(T) = 2.57324  (eps_T = delta / ((2-a) + 2*sqrt(1-a))
  pi = atan(1.) * 4.0;  
  num_periods = 26;
  L       = 320;
  
  /* load into internal parameters */

  par[1+F2C] = lamb;
  par[2+F2C] = eta;
  par[3+F2C] = p;
  par[4+F2C] = nu;
  par[5+F2C] = rho;

  par[17+F2C] = num_periods;
  par[18+F2C] = L;
  X     =  x*L;

  //par[12] = num_periods;
  //par[10] = L;
 
  // The exact Solution
  // the derivatives are relative to the "real" space: L * d/dx_real = d/dx_auto. That's why we divide by L.
  // The L multiplying and the L dividing cancel out, but we'd rather write it this way so we can always remember...
  //dp = 1/((kf+k0)*(kf+k0)-k0*k0)/((kf+k0)*(kf+k0)-k0*k0);
  //dm = 1/((kf-k0)*(kf-k0)-k0*k0)/((kf-k0)*(kf-k0)-k0*k0);
  //amp   =  2.0*sqrt(epsilon + gamma*gamma*epsilon*(dp+dm)/4.0)/sqrt(3.0)/1.0;
  //amp = 1.0*sqrt(epsilon/2.0)/sqrt(3.0);
  //amp = 0.0289334;
  //  amp   =  2.0*sqrt(-epsilon )/sqrt(3.0);
  
  tet=2*pi/L*num_periods;
  amp=1.02;
  //offset=0.021;
  batch=8.5*L/num_periods;
  tat=0.5;
  /*
  U  = 0;
  V  = 0;
//  U  = a/2 - 0.5 * sqrt(a*a-4*m*m);
//  V  = m / W;  

  Ux = 0;
  Vx = 0;
*/

  tet=2*pi/L*num_periods;
  amp=1.02;
  B     =  0;//0.05 - amp * cos(X*tet);
  W     =  p/nu;//-0.1 - amp * cos(X*tet);
 

/*
  U   = 0.05 + amp * cos(X*tet) * 0.5 * (1-tanh(X*tat-batch*tat)) - 0.45 * (1+tanh(X*tat-batch*tat));
  V   = -0.1 + amp * cos(X*tet) * 0.5 * (1-tanh(X*tat-batch*tat)) - 0.45 * (1+tanh(X*tat-batch*tat));
  Ux  = -amp*tet * sin(X*tet) * 0.5 * (1-tanh(X*tat-batch*tat)) - 0.5 * tat*amp*cos(X*tet)/(cosh(X*tat-batch)*cosh(X*tat-batch*tat)) - 0.45/(cosh(X*tat-batch)*cosh(X*tat-batch*tat));
  Vx  = -amp*tet * sin(X*tet) * 0.5 * (1-tanh(X*tat-batch*tat)) - 0.5 * tat*amp*cos(X*tet)/(cosh(X*tat-batch)*cosh(X*tat-batch*tat)) - 0.45/(cosh(X*tat-batch)*cosh(X*tat-batch*tat));
 */


 
  u[0] = B;  
  u[1] = W;  
 
 

  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int pvls (integer ndim, const doublereal *u,
          doublereal *par)
{

  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int bcnd (integer ndim, const doublereal *par, const integer *icp,
          integer nbc, const doublereal *u0, const doublereal *u1, integer ijac,
          doublereal *fb, doublereal *dbc)
{  

  //fb[0] = u1[2];   // Nx_right   = 0
  //fb[1] = u0[2];   // Nx_left    = 0
  //fb[2] = u1[3];   // Wx_right   = 0
  //fb[3] = u0[3];   // Wx_left    = 0

  
  return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int icnd (integer ndim, const doublereal *par, const integer *icp,
          integer nint, const doublereal *u, const doublereal *uold,
          const doublereal *udot, const doublereal *upold, integer ijac,
          doublereal *fi, doublereal *dint)
{
    return 0;
}

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
int fopt (integer ndim, const doublereal *u, const integer *icp,
          const doublereal *par, integer ijac,
          doublereal *fs, doublereal *dfdu, doublereal *dfdp)
{
    return 0;
}
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */

