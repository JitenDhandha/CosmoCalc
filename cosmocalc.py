import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import Planck18, Planck15, Planck13, WMAP9, WMAP7, z_at_value
import astropy.cosmology.units as cu
import astropy.units as u
import argparse

def set_cosmo(cosmo_name):
    if cosmo_name == 'Planck18':
        cosmo = Planck18
    elif cosmo_name == 'Planck15':
        cosmo = Planck15
    elif cosmo_name == 'Planck13':
        cosmo = Planck13
    elif cosmo_name == 'WMAP9':
        cosmo = WMAP9
    elif cosmo_name == 'WMAP7':
        cosmo = WMAP7
    else:
        raise ValueError("Invalid cosmology name. Please choose one of the following: 'Planck18', 'Planck15', 'Planck13', 'WMAP9', 'WMAP7'")
    return cosmo

def print_cosmo(cosmo):
    
    print(f"Name: {cosmo.name}", end=', ')
    print(f"H0 = {cosmo.H0:.4f}", end=', ')
    print(f"Om0 = {cosmo.Om0:.4f}", end=', ')
    print(f"Ode0 = {cosmo.Ode0:.4f}", end=', ')
    print(f"Ok0 = {cosmo.Ok0:.4f}", end=', ')
    print(f"Ob0 = {cosmo.Ob0:.4f}", end=', ')
    print(f"Neff = {cosmo.Neff:.4f}")
    
def z_to_tage(cosmo, z):
    
    if z == 'imp':
        z_imp = [6000,3400,2000,1100,200,80,65,30,20,12,6,0.3] * cu.redshift
        tage_imp = cosmo.age(z_imp).to(u.Myr)
        print("Important redshifts:")
        for i in range(len(z_imp)):
            print(f"z = {z_imp[i]}, tage = {tage_imp[i]:.4f}")
            
    elif z == 'plot':
        zmin = 1e-1
        zmax = 2e4
        z_linspace = np.logspace(np.log10(zmax),np.log10(zmin),1000) * cu.redshift
        tage_linspace = cosmo.age(z_linspace).to(u.Myr)
        fig = plt.figure(figsize=(12,6))
        ax = fig.add_subplot(111)
        ax.plot(z_linspace, tage_linspace, color='black', linewidth=2.5)
        ax.axvspan(3400,zmax,color='red', alpha=0.20)
        ax.text(0.075,1.02,'Radiation\ndominated',fontsize=15,transform=ax.transAxes,ha='center')
        ax.axvspan(0.3,3400,color='green', alpha=0.20)
        ax.text(0.5,1.02,'Matter\ndominated',fontsize=15,transform=ax.transAxes,ha='center')
        ax.axvspan(zmin,0.3,color='blue', alpha=0.20)
        ax.text(0.96,1.02,'Lambda\ndominated',fontsize=15,transform=ax.transAxes,ha='center')
        ax.grid(visible=True, which='both', axis='both', color='grey', alpha=0.25, linestyle='-')
        ax.tick_params(axis='both', which='both', direction='in', top=True, right=True, labelsize=15)
        ax.set_xlabel(r'Redshift $z$', fontsize=17)
        ax.set_ylabel(r'Age of the Universe $t_{\rm age}$ [Myr]', fontsize=17)
        ax.invert_xaxis()
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim([zmax,zmin])
        fig.show()
    else:
        try: 
            z = float(eval(z)) * cu.redshift
            tage = cosmo.age(z).to(u.Myr)
            print(f"z = {z}, tage = {tage:.4f}")
        except (ValueError, TypeError):
            raise ValueError("Invalid redshift. Please provide a valid redshift, 'imp' or 'plot'")
        
def tage_to_z(cosmo, tage):
    try:
        tage = float(eval(tage)) * u.Myr
        z = z_at_value(cosmo.age, tage)
        print(f"tage = {tage}, z = {z:.4f}")
    except (ValueError, TypeError):
        raise ValueError("Invalid age of the Universe. Please provide a valid age of the Universe in Myrs.")

def H_at_z(cosmo, z):
    try: 
        z = float(eval(z))
        Hz = cosmo.H(z)
        print(f"z = {z}, H(z) = {Hz:.4f}")
    except (ValueError, TypeError):
        raise ValueError("Invalid redshift. Please provide a valid redshift.")
    
def comoving_to_deg(cosmo, z, x):
    try:
        z = float(eval(z)) * cu.redshift
        x = float(eval(x)) * u.Mpc
        theta = cosmo.arcsec_per_kpc_comoving(z) * x.to(u.kpc)
        theta = theta.to(u.deg)
        print(f"z = {z}, x = {x}, theta = {theta:.4f}")
    except (ValueError, TypeError):
        raise ValueError("Invalid redshift or comoving distance. Please provide a valid redshift and comoving distance.")

def deg_to_comoving(cosmo, z, theta):
    try:
        z = float(eval(z)) * cu.redshift
        theta = float(eval(theta)) * u.deg
        x = theta.to(u.arcsec) / cosmo.arcsec_per_kpc_comoving(z)
        x = x.to(u.Mpc)
        print(f"z = {z}, theta = {theta:.4f}, x = {x:.4f}")
    except (ValueError, TypeError):
        raise ValueError("Invalid redshift or angular size. Please provide a valid redshift and angular size.")
    
def deg2_to_comoving_area(cosmo, z, omega):
    try:
        z = float(eval(z)) * cu.redshift
        omega = float(eval(omega)) * u.deg**2
        A = omega.to(u.arcsec**2) / cosmo.arcsec_per_kpc_comoving(z)**2 
        A = A.to(u.Mpc**2)
        print(f"z = {z}, omega = {omega:.4f}, A = {A:.4f}")
    except (ValueError, TypeError):
        raise ValueError("Invalid redshift or angular area. Please provide a valid redshift and angular area.")
    
def deg2_to_comoving_volume(cosmo, z1, z2, omega):
    try:
        z1 = float(eval(z1)) * cu.redshift
        z2 = float(eval(z2)) * cu.redshift
        omega = float(eval(omega)) * u.deg**2
        V = omega.to(u.steradian)/(4*np.pi*u.steradian) * abs(cosmo.comoving_volume(z1) - cosmo.comoving_volume(z2))
        print(f"z1 = {z1}, z2 = {z2}, omega = {omega:.4f}, V = {V:.4f}")
    except (ValueError, TypeError):
        raise ValueError("Invalid redshifts or angular area. Please provide valid redshifts and angular area.")
        
def main():
    
    functions = {
        'z_to_tage': ['Calculate the age of the Universe at a given redshift', z_to_tage],
        'tage_to_z': ['Calculate the redshift at a given age of the Universe', tage_to_z],
        'H_at_z': ['Calculate the Hubble parameter at a given redshift', H_at_z],
        'print_cosmo': ['Print the cosmological parameters', print_cosmo],
        'comoving_to_deg': ['Convert comoving distance to angular size', comoving_to_deg],
        'deg_to_comoving': ['Convert angular size to comoving distance', deg_to_comoving],
        'deg2_to_comoving_area': ['Convert angular area to comoving area', deg2_to_comoving_area],
        'deg2_to_comoving_volume': ['Convert angular area to comoving volume', deg2_to_comoving_volume],
    }

    parser = argparse.ArgumentParser(description='Cosmology Calculator Tool, by Jiten Dhandha')
    parser.add_argument('-f','--function', type=str, help=f'Function to calculate: {", ".join(functions.keys())}')
    parser.add_argument('-c','--cosmology', type=str, default='Planck18', help='Name of the cosmology to use')
    parser.add_argument('-z','--redshift', type=str, help="Input redshift")
    parser.add_argument('-z2','--redshift2', type=str, help="Input redshift 2")
    parser.add_argument('-t','--tage', type=str, help='Age of the Universe in Myr')
    parser.add_argument('-x','--comoving', type=str, help='Comoving distance in Mpc')
    parser.add_argument('-theta','--angular_sep', type=str, help='Angular size in deg')
    parser.add_argument('-omega','--angular_area', type=str, help='Angular area in deg^2')
    args = parser.parse_args()
    cosmo = set_cosmo(args.cosmology)
    
    if args.function not in functions:
        raise ValueError(f"Invalid function. Please choose one of the following: {', '.join(functions.keys())}")
    else:
        f = functions[args.function][1]
        if args.function == 'z_to_tage':
            f(cosmo, args.redshift)
        elif args.function == 'tage_to_z':
            f(cosmo, args.tage)
        elif args.function == 'comoving_to_deg':
            f(cosmo, args.redshift, args.comoving)
        elif args.function == 'deg_to_comoving':
            f(cosmo, args.redshift, args.angular_sep)
        elif args.function == 'deg2_to_comoving_area':
            f(cosmo, args.redshift, args.angular_area)
        elif args.function == 'deg2_to_comoving_volume':
            f(cosmo, args.redshift, args.redshift2, args.angular_area)
        elif args.function == 'H_at_z':
            f(cosmo, args.redshift)
        elif args.function == 'print_cosmo':
            f(cosmo)
        
    plt.show()
    
if __name__ == '__main__':
    main()