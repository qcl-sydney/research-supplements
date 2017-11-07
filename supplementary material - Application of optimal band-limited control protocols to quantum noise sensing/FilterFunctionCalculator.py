# -*- coding: utf-8 -*-
from __future__ import print_function, division
import numpy as np


def FF_Amplitude(RabiRates, w_vec, tau, FFT=False):
    ''' Calculate the amplitude Filter function (FF) for a control sequence around a single,
        fixed rotation axis (x or y) as defined by the time-dependent Rabi rates.
        To calculate the Rabi rates, see functions RabiRates_SLEP or RabiRates_PRIM.
        
     Parameters
     ----------------------------
        RabiRates : array
                1D array containing the Rabi rates of the control sequence
                
        w_vec : array
                1D vector of angular frequencies over which the Filter function is to be calculated
         
        tau : float
                Time of the full control sequence
        
        FFT : bool, optional
            Whether to perform the calculation via a fast Fourier transform
         
         
     Returns 
     ----------------------------
            FF : array
                1D array containing the resulting Filter function
                
        
     Notes 
     ----------------------------
        For a control sequence around a single rotation axis only, the Filter function is given
        simply by the Fourier transform of the time-domain Rabi rates.
        
     Examples 
     ----------------------------
        
        Calculate and plot the Filter function for a rotary spin echo sequence with 5 sign flips, 
        1ms pi time and 1ms sequence duration.
    
        >>> import matplotlib.pyplot as plt
        >>> rabi_prim = RabiRates_PRIM(10, 1e-3, n_segments=5)
        >>> w_vec = np.arange(0, 30e3, 0.5e3)
        >>> FF_prim = FF_Amplitude(rabi_prim, w_vec, 1e-3)
        >>> fig, ax = plt.subplots()
        >>> ax.plot(w_vec, FF_prim/w_vec**2)
        >>> plt.show() 
    '''
    FF = np.zeros(len(w_vec))
    N = RabiRates.shape[0]
    pre_factor = np.sin(w_vec*tau/(2*N))**2
    t = np.arange(N)
    if FFT is True:
        n_f = len(w_vec); n_t = len(RabiRates)
        delt_f = w_vec[-1] / 4 /np.pi / n_f; delt_t = tau / n_t
        n_zeros = int(1 / (2 * delt_f * delt_t) - tau / delt_t)
        fft = np.fft.fft(RabiRates, n=n_t + n_zeros)
        FF = np.abs(fft[:n_f])**2 * pre_factor
    else:
        FF = pre_factor * np.abs(np.sum(RabiRates[:, np.newaxis]*np.exp(1j*w_vec*t[:, np.newaxis]*tau/N), axis=0))**2
    return FF


def RabiRates_PRIM(N, piTime, tau=0, theta=np.pi, n_segments=1):
    '''
        Calculate the correctly normalised Rabi rates for a flat-top control sequence around a 
        single, fixed rotation axis (x or y) with (optionally) alternating signs in the amplitude, 
        to construct a rotary echo sequence.


    Parameters
    ----------------------------
        N : integer
                Number of points per segment
                
        piTime : float
                Time required to perform a full pi rotation, in seconds
         
        tau : float, optional
                Time of the full control sequence, if different from piTime
                
        theta : float, optional
                Rotation angle to be swept out by the protocol
                
        n_segments : integer
                Number of segments with alternating amplitudes (+1/-1)

    Returns
    ----------------------------
    
        RabiRates : array
                Array containing the correctly normalised Rabi rates for a primitive pulse sequence
        
    Examples
    ----------------------------
    
      Calculate and plot the Rabi rates for a rotary spin echo sequence with 5 sign flips, 1ms pi time
    
      >>> import matplotlib.pyplot as plt
      >>> rabi_prim = RabiRates_PRIM(10, 1e-3, n_segments=5)
      >>> fig, ax = plt.subplots()
      >>> ax.plot(rabi_slep)
      >>> plt.show()
    '''

    RabiRates = -np.ones(n_segments)
    
    for i in range(0, n_segments):
        RabiRates[i] = RabiRates[i]**(i+2)
        
    RabiRates *= theta/(piTime if tau == 0 else tau)
    
    return np.repeat(RabiRates, N) 
    
def RabiRates_SLEP(dpss, piTime, tau=0, theta=np.pi, SSB=False, 
                   xfunc=None, w_0=0, normalise=False):
    '''
        Calculate the correctly normalised Rabi rates for a control sequence around a single,
        fixed rotation axis (x or y) with Slepian-modulated amplitudes. The units of Rabi rates 
        are angle (rad) / time (s), so this function needs information about the total time of 
        the sequence. In addition to the DPSS amplitude, this function provides the option to 
        apply additional amplitude modulation  in form of cosinusoidal or single-sideband 
        modulation to shift the frequency response of the control sequence.

    Parameters
    ----------------------------
        dpss : 1D array 
                DPSS envelope
                
        piTime : float
                Time required to perform a full pi rotation, in seconds
         
        tau : float, optional
                Time of the full control sequence, if different from piTime
                
        theta : float, optional
                Rotation angle to be swept out by the protocol
                
        SSB : bool, optional
                (AM) Apply single-sideband modulation to the DPSS envelope with frequency w_0
            
        xfunc : function, optional
                (AM) Apply AM with supplied function (accepts sine or cosine) with frequency w_0
                
        w_0 : float, optional
                (AM) Frequency for xfunc or SSB
          
        normalise : bool, optional
                Whether or not to normalise the Slepian pulse such that the FF areas are equal
                (becomes important under application of AM, SSB AM)
        

    Returns
    ----------------------------
    
        RabiRates : array
                Array containing the correctly normalised, optionally amplitude-modulated Rabi rates
        
    Examples
    ----------------------------
    
      Calculate and plot the Rabi rates for a 0th order DPSS with a 1ms pi time
    
      >>> import matplotlib.pyplot as plt
      >>> from slepian import DPSS
      >>> dpss = DPSS(100, 1, 3)[:, 0]
      >>> rabi_dpss = RabiRates_SLEP(dpss, 1e-3)
      >>> fig, ax = plt.subplots()
      >>> ax.plot(rabi_dpss)
      >>> plt.show()
    '''
    RabiRates = dpss
    N = RabiRates.shape[0]
    t = np.arange(N)*(tau if tau != 0 else piTime)/N

    # Apply optional amplitude modulation
    if xfunc is not None and SSB is False:
        RabiRates *= xfunc(w_0*t)
    elif SSB is True:
        arg = w_0*np.arange(N)*piTime/N;
        RabiRates_H = np.imag(hilbert(RabiRates))
        RabiRates = RabiRates*np.cos(arg) - RabiRates_H*np.sin(arg)

    if not normalise:
        RabiRates *= 1 / np.sum(np.abs(RabiRates/(N)))
    else:
        if xfunc is not None and w_0 == 0:
            RabiRates *= 2 / (1*np.sum(np.abs(RabiRates/(N))))
        elif xfunc is not None and w_0 != 0:
            RabiRates *= 1 / np.sum(np.abs(RabiRates/(N)))
        elif SSB is True:
            RabiRates *= 1 / np.sum( np.abs(RabiRates/(N)) )**2

    RabiRates *= theta/(piTime if tau == 0 else tau)

    return RabiRates
    
    
if __name__ == '__main__':

    import matplotlib    
    import matplotlib.pyplot as plt
    from Slepian import DPSS
    
    if matplotlib.__version__.startswith('1'):
        colours = [c['color'] for c in list(plt.rcParams['axes.prop_cycle'])]
    else:
        prop_cycle = plt.rcParams['axes.prop_cycle']
        colours = prop_cycle.by_key()['color']

    pitime = 1e-3
    w_vec = np.linspace(1,30000,10000)

    k_max = 6
    N = 100
    dpss = DPSS(N, k_max, k_max)
    
    fig, ax = plt.subplots()
    for i in range(0, k_max, 2):
        r_slep  = RabiRates_SLEP(dpss[:, i], pitime)
        r_prim  = RabiRates_PRIM(N, pitime, n_segments=i+1) 
        
        FF_slep = FF_Amplitude(r_slep, w_vec, pitime)
        FF_prim = FF_Amplitude(r_prim, w_vec, pitime)
    
        ax.plot(w_vec/2/np.pi, FF_slep/w_vec**2, color=colours[i])
        ax.plot(w_vec/2/np.pi, FF_prim/w_vec**2, color=colours[i], ls='--')
    
    ax.set_xlabel('Frequency (Hz)')
    ax.set_ylabel('Filter function amplitude (a.u.)')
    plt.show()
