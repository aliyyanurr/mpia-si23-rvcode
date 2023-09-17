import numpy as np
import matplotlib.pyplot as plt
import kepler 

G = 6.674*1e-11 #m^3 kg^-1 s^-2
s_in_year = 3.154e+7 # s_in_year s = 1 year
Msun_in_kg = 1.98847e+30 # Msun_in_kg kg = 1 Msun
Mjup_in_kg = 1.898e+27 # Mjup_in_kg kg = 1 Mjup
au_in_m = 149597870700 # au_in_m m = 1 au


# THIS BLOCK OF CODE IS A PERIOD-TO-SEMIMAJOR AXIS CONVERTER (KEPLER 3RD LAW)

def atop_sunjup(a, M, m):
    """
    Converts semimajor axis to period, with the units in au, solar mass (larger), and jupiter mass (smaller).
    Returning period in the unit of second.
    Usage: atop_sunjup(semimajor_axis, larger_mass, smaller_mass)
    """
    P = (( (a*au_in_m)**3 * 4 * np.pi**2)/(G*((M*Msun_in_kg)+(m*Mjup_in_kg))))**0.5
    print(f'P = ', P / s_in_year, ' years')
    return P

def atop_sunsun(a, M, m):
    """
    Converts semimajor axis to period, with the units in au, solar mass (larger and smaller).
    Returning period in the unit of second.
    Usage: atop_sunsun(semimajor_axis, larger_mass, smaller_mass)
    """
    P = (( (a*au_in_m)**3 * 4 * np.pi**2)/(G*((M*Msun_in_kg)+(m*Msun_in_kg))))**0.5
    print(f'P = ', P / s_in_year, ' years')
    return P
    
def atop_si(a, M, m):
    """
    Converts semimajor axis to period, with the units in m and kg (all bodies).
    Returning period in the unit of second.
    Usage: atop_si(semimajor_axis, larger_mass, smaller_mass)
    """
    P = (( (a)**3 * 4 * np.pi**2)/(G*((M)+(m))))**0.5
    #print(f'P = ', P / s_in_year, ' years')
    return P #in second

def ptoa_sunjup(P,M,m):
    """
    Converts period to semimajor axis, with the units in year, solar mass (larger), and jupiter mass (smaller).
    Returning semimajor axis in the unit of m.
    Usage: ptoa_sunjup(period, larger_mass, smaller_mass)
    """
    a = (((P*s_in_year)**2 * G * ((M*Msun_in_kg)+(m*Mjup_in_kg)))/(4 * np.pi**2))**(1/3)
    print(f'a = ', a/au_in_m, 'au')
    return a

def ptoa_si(P,M,m):
    """
    Converts period to semimajor axis, with the units in s and kg (all bodies).
    Returning semimajor axis in the unit of m.
    Usage: ptoa_si(period, larger_mass, smaller_mass)
    """
    a = ((P**2 * G * (M+m))/(4 * np.pi**2))**(1/3)
    print(f'a = ', a/au_in_m, 'au')
    return a

# THIS BLOCK OF CODE IS FOR RADIAL VELOCITY SIMULATION

def arcocos(costheta):
    """
    Calculate the arccosine of an angle to give a full-circle angle (0-2pi).
    Input: cos_angle in float list
    Output: angle (rad)
    Usage: arcocos(cos_angle)
    """
    theta = np.arccos(costheta)
    tarr = []
    for i in range(0,len(theta)):
        if i == 0:
            tarr.append(theta[i])
        elif i == 1:
            selisih = theta[i] - theta[i-1]
            indi = np.sign(selisih)
            tarr.append(theta[i])
        else:
            selisih = theta[i] - theta[i-1]
            indi2 = np.sign(selisih)
            if indi2 != indi:
                theta = 2*np.pi - theta
                tarr.append(theta[i])
            else:
                tarr.append(theta[i])
    tarr = np.array(tarr)
    return(tarr)

def avg_ang_speed(per):
    """
    Calculate average angular speed of the smaller body.
    Input: period (s).
    Output: average angular speed, in rad/s.
    Usage: avg_ang_speed(period)
    """
    avg_ang_sp = 2*np.pi/per # rad / s
    return avg_ang_sp

def mean_anomaly(t, avg_ang_sp):
    """
    Calculate the mean anomaly of the smaller body, at each time in the orbit.
    Input: time (s), average angular speed (rad/s)
    Output: mean anomaly (rad)
    Usage: mean_anomaly(time, avg_ang_speed)
    """
    mean_anom = avg_ang_sp * t
    return mean_anom

def kepler_solve(mean_anom, ecc):
    """
    Solve the Kepler equation.
    Input: mean anomaly (rad), eccentricity
    Output: eccentric anomaly (rad), true anomaly (rad)
    Usage: kepler_solve(mean_anomaly, eccentricity)
    """
    ecc_anom, cos_TA, sin_TA = kepler.kepler(mean_anom, ecc)
    true_anom = arcocos(cos_TA)
    return ecc_anom, true_anom

def pos(sma, ecc, true_anomaly):
    """
    Calculate the position (r, theta) of a body in the orbit.
    Input: semimajor axis (any unit), eccentricity, true anomaly (rad)
    Output: r (related unit), true anomaly (rad) 
    Usage: pos(semimajor_axis, eccentricity, true_anomaly)
    """
    r = sma * (1 - ecc**2) / (1 + ecc * np.cos(true_anomaly))
    return r, true_anomaly

def radial_velocity(true_anom, omega, ecc, gamma, sma, incl, mass1, mass2):
    """
    Calculate the radial velocity of the smaller body around the larger body.
    Input: true anomaly (rad), omega (rad), eccentricity, gamma (m/s), period (s), inclination (rad), larger mass (kg), smaller mass (kg)
    Output: radial velocity (m/s)
    Usage: radial_velocity(true_anomaly, omega, eccentricity, gamma, period, inclination, larger_mass, smaller_mass)
    """
    G = 6.6743*1e-11 #m3 kg-1 s-2
    K = np.sqrt(G/(1-ecc**2))*mass2*np.sin(incl)*((mass1+mass2)**(-0.5))*(sma**(-0.5))
    rad_vel = K * (np.cos(true_anom + omega) + (ecc * np.cos(omega))) + gamma
    return rad_vel

def ind_radial_velocity(true_anom, omega, ecc, gamma, sma, incl, mass_obj, mass_total):
    """
    Calculate the individual radial velocity of each body towards the system barycenter. Can be used to calculate the relative radial velocity.
    Input: true anomaly (rad), omega (rad), eccentricity, gamma (m/s), period (s), inclination (rad), larger mass (kg), smaller mass (kg)
    Output: individual radial velocity (m/s)
    Usage: radial_velocity(true_anomaly, omega, eccentricity, gamma, semimajor axis of the object, inclination, larger_mass, smaller_mass)
    """
    G = 6.6743*1e-11 #m3 kg-1 s-2
    K = np.sqrt(G/(1-ecc**2))*mass_obj*np.sin(incl)*((mass_total)**(-0.5))*(sma**(-0.5))
    rad_vel = K * (np.cos(true_anom + omega) + (ecc * np.cos(omega))) + gamma
    return rad_vel