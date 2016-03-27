# -*- coding: utf-8 -*-
from __future__ import print_function
import argparse
import PyDSTool as dst
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot as plt
mpl.rcParams['legend.fontsize'] = 10


def model():
    K = 0.4
    E = 7.0
    M = 10.5
    N = 15
    Lambda = 0.9
    Gamma  = 12
    R  = 0.7
    PP = 20
    
    lamb_p    = (K*Gamma)/M
    eta_p     = E*K
    p_p       = (Lambda*PP)/(K*Gamma*M)
    nu_p     = N/M
    rho_p     = R
    
    # Declare names and initial values for (symbolic) parameters
    lamb   = dst.Par(lamb_p, 'lamb')
    eta    = dst.Par(eta_p, 'eta')
    p      = dst.Par(p_p, 'p')
    nu     = dst.Par(nu_p, 'nu')
    rho    = dst.Par(rho_p, 'rho')
    
    # Compute nontrivial boundary equilibrium initial condition from parameters (see reference)
    
    
    b_0 = 0.0
    w_0 = p_p/nu_p
    # Declare symbolic variables    
    b = dst.Var('b')
    w = dst.Var('w')
    t = dst.Var('t')

    # Create Symbolic Quantity objects for definitions
    brhs = dst.Fun(lamb*w*b*((1+eta*b)**2)*(1-b) - b,[b,w],'brhs')
    wrhs = dst.Fun(p - nu*w*(1-rho*b) - lamb*w*b*((1+eta*b)**2),[b,w],'wrhs')

    F   = dst.Fun([brhs(b,w),wrhs(b,w)], [b,w], 'F')
    jac = dst.Fun(dst.Diff(F,[b,w]), [t,b,w], 'Jacobian')


    # Build Generator
    DSargs = dst.args(name='fairy_circles_ode')
    DSargs.fnspecs = [jac, brhs,wrhs]
    DSargs.varspecs = {b:brhs(b,w) ,
                       w:wrhs(b,w)}
    DSargs.pars = [lamb,eta,p,nu,rho]
    
    
    # Use eval method to get a float value from the symbolic definitions given in
    # terms of parameter values
    DSargs.ics = dst.args(b=b_0, w=w_0)
    return DSargs

def cont(model, maxnum=450,maxstep=2.0,minstep=1e-5,stepsize=2e-2,direction="forward"):
    ode  = dst.Generator.Vode_ODEsystem(model)
    # Prepare the system to start close to a steady state
#    ode.set(pars = {'p': 0.078} )       # Lower bound of the control parameter 'i'
#    ode.set(ics =  {'b': 0.0, 'w': 0.0} )   # Close to one of the steady states present for i=-220
    
    PC = dst.ContClass(ode)            # Set up continuation class
    
    PCargs = dst.args(name='EQ1', type='EP-C')     # 'EP-C' stands for Equilibrium Point Curve. The branch will be labeled 'EQ1'.
    PCargs.freepars     = ['p']                    # control parameter(s) (it should be among those specified in DSargs.pars)
    PCargs.MaxNumPoints = maxnum                      # The following 3 parameters are set after trial-and-error
    PCargs.MaxStepSize  = maxstep
    PCargs.MinStepSize  = minstep
    PCargs.StepSize     = stepsize
    PCargs.LocBifPoints = 'LP'                     # detect limit points / saddle-node bifurcations
    PCargs.SaveEigen    = True
    PC.newCurve(PCargs)
    if direction == "forward":
        PC['EQ1'].forward()
    elif direction=="backward":
        PC['EQ1'].backward()
    PC.display(['p','b'], stability=True, figure=3)        # stable and unstable branches as solid and dashed curves, resp.
    return PC


def main(args):
    global ode_model
    ode_model = model()
    return 0

def add_parser_arguments(parser):
    "Add command-line arguments."""
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        dest="verbose",
                        default=False,
                        help="Turn on debuging messages.")
    return parser

if __name__ == '__main__':
    parser = argparse.ArgumentParser(prog='PROG', usage='%(prog)s [options]')
    parser = add_parser_arguments(parser)
    args = parser.parse_args()
    main(args)