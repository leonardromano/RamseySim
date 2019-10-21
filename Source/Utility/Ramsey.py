import numpy as np
import Source.Utility.Math_Utility as mut

from Source.Utility.Constants import *

def polarizationVector(delta_omega, t_p, omega_rabi, Omega_eff, \
                       T = T, T1 = T1, alpha_0 = alpha_0):
    down = ((alpha_0*np.exp(-T/T1)*omega_rabi**2)/(4*Omega_eff**4)) *\
    (np.sin(Omega_eff*t_p))**2 *(2*Omega_eff*np.cos(Omega_eff*t_p)*\
     np.cos(delta_omega/2 * T) - delta_omega*np.sin(Omega_eff*t_p)*\
     np.sin(delta_omega/2 * T))**2
     
    return np.array([1-down, down])

def analytical(omega, omega_larmor, T = T, T1 = T1, \
                      offset = 0, scale = 1, alpha_0 = alpha_0,): 
    
    delta_omega, t_p, omega_rabi, Omega_eff = \
    mut.calculate_effective_frequency(omega, omega_larmor)
    
    prob_up  = polarizationVector(delta_omega, t_p, omega_rabi, Omega_eff, \
                                  T, T1, alpha_0)[0]
    polarization = 2*prob_up - 1
         
    return (polarization + offset)*scale

def analytical_derivative(omega, omega_larmor, T = T, T1 = T, \
                          offset = 0, scale = 1, alpha_0 = alpha_0, ):
        
    delta_omega, t_p, omega_rabi, Omega_eff = \
    mut.calculate_effective_frequency(omega, omega_larmor)

    B = 2*Omega_eff*np.cos(Omega_eff*t_p)*np.cos(delta_omega/2 * T) - \
    delta_omega*np.sin(Omega_eff*t_p)*np.sin(delta_omega/2 * T)
        
    dB_domega = (delta_omega/(2*Omega_eff))*np.cos(Omega_eff*t_p)*\
    np.cos(delta_omega/2 * T) - np.sin(Omega_eff*t_p)*\
    np.sin(delta_omega/2 * T) - (delta_omega/2)*(T+t_p)*\
    np.sin(Omega_eff*t_p)*np.cos(delta_omega/2 * T) - \
    (Omega_eff*T + delta_omega**2/(4*Omega_eff)*t_p)*\
    np.cos(Omega_eff*t_p)*np.sin(delta_omega/2*T)

    A = alpha_0 * np.exp(-T/T1) * omega_rabi**2 / 2
        
    dP_domega = A*B*np.sin(Omega_eff*t_p) * ((B*delta_omega/Omega_eff**6)*\
                           np.sin(Omega_eff*t_p) - \
                           (B*delta_omega*t_p/(2*Omega_eff**5)) * \
                           np.cos(Omega_eff*t_p) -  2*dB_domega*\
                           (np.sin(Omega_eff*t_p)/Omega_eff**4))
        
    return dP_domega*scale

def residual(params, omega, polarization, sigma_x, sigma_y): 
    omega_larmor = params['omega_larmor'].value
    T = params['T'].value
    T1 = params['T1'].value
    offset = params['offset'].value
    scale = params['scale'].value
            
    model = analytical(omega, omega_larmor, T, T1, offset, scale)
        
    if sigma_x is not None and sigma_y is not None:
        derivative = analytical_derivative(omega, omega_larmor, \
                                                       T, T1)
        return (polarization - model) / np.sqrt(sigma_y**2 + \
               (derivative*sigma_x)**2)
    elif sigma_x is not None and sigma_y is None:
        derivative = analytical_derivative(omega, omega_larmor, \
                                                           T, T1)
        return (polarization - model) / (derivative*sigma_x)
    elif sigma_x is None and sigma_y is not None:
        return (polarization - model) / sigma_y
    elif sigma_x is None and sigma_y is None:
        return (polarization - model)