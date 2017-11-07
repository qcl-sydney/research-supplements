from __future__ import division, print_function
import numpy as np


def DPSS(N, k_max, NW, EW=False):
    '''
    Calculate the discrete prolate spheroidal sequences for all orders up to k_max using Newton's method
    and a semi-empiric expression to solve for the eigenvalues of the Toelplitz matrix and then extract the
    corresponding eigenvector. This algorithm is taken from the book "Numerical Recipes in C++ - The art of
    scientific computing" by W.H. Press and others. 
    Parameters
    ----------------------------
    
        N : integer
            Number of points in the sequence
            
        k_max : integer
            Maximum order up until which the sequences are calculated
            
        NW : integer
            The half-time bandwidth product
            
        EW : bool, optional
        Whether to return eigenvalues as well as the eigenvectors
    Returns
    ----------------------------
    
        dpss : array
        Array containing the sequences for orders up to k_max; dim (N, k_max).
        
    Examples
    ----------------------------
    
    Calculate and plot the first couple of DPSS orders
    
    >>> import matplotlib.pyplot as plt
    >>> dpss = DPSS(100, 5, 5)
    >>> dpss.shape
    (100, 5)
    >>> fig, ax = plt.subplots()
    >>> for i in range(5):
    ...    ax.plot(dpss[:, i], label=r'$k={}$'.format(i))
    >>> ax.legend(loc='lower left')
    >>> plt.show()
    
    '''
    # Define vectors to hold the diagonal and off-diagonal elements of the Toelpitz matrix as well
    # as vectors to hold intermediate results
    dg = np.zeros(N, dtype=np.float64); dgg = np.zeros(N, dtype=np.float64); gam = np.zeros(N, dtype=np.float64);
    sup = np.zeros(N-1, dtype=np.float64); sub = np.zeros(N, dtype=np.float64); 

    # Start filling the Toelplitz matrix
    n = np.arange(N)
    dg =    0.25*(N**2 - (N-1-2*n)**2 * np.cos(2*np.pi*NW/N) )
    sup = - n[1:]*(N-n[1:])/2.0; sub = - n[1:]*(N-n[1:])/2.0

    # Initial, semi-empirical guess for eigenvalues
    xx = -0.10859 - 0.068762/NW + 1.5692*NW
    xold = xx + 0.47276 + 0.20273/NW - 3.1387*NW

    u = np.zeros(N)
    dpss = np.zeros((N, k_max))
    ews = np.zeros(k_max)

    for k in range(k_max):  # Loop over orders
        u = dpss[:, k]      # Point the vector u into the dpss array

        # Find the eigenvalue for this order
        for i in range(20):             # Loop over iterations for Newton's algorithm
            pp = 1.;   p = dg[0] - xx;   dd = 0.;    d = -1.

            for j in range(1,N):
                ppp = pp; pp = p;   ddd = dd;   dd = d

                p = pp*(dg[j] - xx) - ppp*sup[j-1]**2
                d = -pp + dd*(dg[j] - xx) - ddd*sub[j-1]**2

                if abs(p) > 1.e30: 
                    p = IDexp(p, -100); pp = IDexp(pp, -100);
                    d = IDexp(d, -100); dd = IDexp(dd, -100);
               
                elif abs(p) < 1.e-30:
                    p = IDexp(p, 100); pp = IDexp(pp, 100);
                    d = IDexp(d, 100); dd = IDexp(dd, 100);

            # Adjust the eigenvalue (Newton's method)
            xnew = xx - p/d             
            
            # Stop iteration when desired precision is reached            
            if abs(xx-xnew) < 1e-10*abs(xnew): break    
                
            xx = xnew
        ews[k] = xnew
        xx = xnew - (xold - xnew)
        xold = xnew

        # Subtract eigenvalue from diagonal, then set one component
        # fixed (but save the old value) 
        dgg = dg - xnew         
        
        nl = N // 3            
        dgg[nl] = 1.
        ssup = sup[nl];  ssub = sub[nl-1]
        
        u[0] = 0.; sup[nl] = 0.;  sub[nl-1] = 0.;

        bet = dgg[0]
        
        # Find the eigenvector corresponding to the eigenvalue
        for i in range(1, N):           
            gam[i] = sup[i-1]/bet
            bet = dgg[i] - sub[i-1]*gam[i]

            s = 1. if i == nl else 0.
            u[i] = (s - sub[i-1]*u[i-1])/bet

        for i in range(N-2, -1, -1):
            u[i] -= gam[i+1]*u[i+1]

        # Restore saved values
        sup[nl] = ssup;   sub[nl-1] = ssub      

        # Renormalise and abide to the sign convention
        s = np.sum(u**2)                      
        s = np.sqrt(s) if u[3] > 0. else - np.sqrt(s)
        u /= s
    
    if EW:
        return dpss, ews   
    else:
        return dpss
        

def IDexp(x, exp):
    return x * 2**exp

        

if __name__ == '__main__':
    
    import matplotlib.pyplot as plt
    dpss = DPSS(100, 5, 5)
    
    fig, ax = plt.subplots()
    for i in range(5):
        ax.plot(dpss[:, i], label=r'$k={}$'.format(i))
    ax.legend(loc='lower left')
plt.show()
