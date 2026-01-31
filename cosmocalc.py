# Core packages for cosmology calculations
import numpy as np
import matplotlib.pyplot as plt
from astropy.cosmology import Planck18, Planck15, Planck13, WMAP9, WMAP7, FlatLambdaCDM, z_at_value
import astropy.cosmology.units as cu
import astropy.units as u

# CLI packages
from typing import Annotated
from cyclopts import App, Parameter
from rich.console import Console
from rich.panel import Panel
from rich.text import Text

app = App(help="Cosmology Calculator Tool, by Jiten Dhandha", help_format="rich", print_error=True)
console = Console()

class CLIError(Exception):
    def __init__(self, message: str, *, title: str = "Error"):
        super().__init__(message)
        self.message = message
        self.title = title

    def render(self):
        return Panel(
            Text(self.message),
            title=self.title,
            title_align="left",
            border_style="red",
        )

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
        try:
            cosmo = FlatLambdaCDM(eval(cosmo_name))
        except:
            raise CLIError("Invalid cosmology name.")
    return cosmo
    
@app.command
def print_cosmo(
    cosmology: Annotated[str, Parameter(name=['-c', '--cosmology'])] = "Planck18"
):
    """
    Print the cosmological parameters.
    
    Parameters
    ----------
    cosmology : str
        Cosmology model to use
    """
    cosmo = set_cosmo(cosmology)
    print(f"Name: {cosmo.name}", end=', ')
    print(f"H0 = {cosmo.H0:.4e}", end=', ')
    print(f"Om0 = {cosmo.Om0:.4e}", end=', ')
    print(f"Ode0 = {cosmo.Ode0:.4e}", end=', ')
    print(f"Ok0 = {cosmo.Ok0:.4e}", end=', ')
    print(f"Ob0 = {cosmo.Ob0:.4e}", end=', ')
    print(f"Neff = {cosmo.Neff:.4e}")
        
@app.command
def z_to_tage(
    z: Annotated[str, Parameter(name=['-z', '--redshift'])],
    cosmology: Annotated[str, Parameter(name=['-c', '--cosmology'])] = "Planck18"
):
    """
    Calculate the age of the Universe (in Myr) at a given redshift.
    
    Parameters
    ----------
    z : str
        Redshift value ('imp' for important redshifts, 'plot' for plotting, or numeric value)
    cosmology : str
        Cosmology model to use
    """
    cosmo = set_cosmo(cosmology)
    if z == 'imp':
        z_imp = [6000,3400,2000,1100,200,80,65,30,20,12,6,0.3] * cu.redshift
        tage_imp = cosmo.age(z_imp).to(u.Myr)
        print("Important redshifts:")
        for i in range(len(z_imp)):
            print(f"z = {z_imp[i]}, tage = {tage_imp[i]:.4e}")
            
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
            print(f"z = {z}, tage = {tage:.4e}")
        except (ValueError, TypeError, NameError):
            raise CLIError("Invalid redshift. Please provide a valid redshift, 'imp' or 'plot'")
    
@app.command  
def tage_to_z(
    tage: Annotated[str, Parameter(name=['-t', '--tage'])],
    cosmology: Annotated[str, Parameter(name=['-c', '--cosmology'])] = "Planck18"
):
    """
    Calculate the redshift at a given age of the Universe (in Myr).

    Parameters
    ----------
    tage : str
        Age of the Universe in Myr
    cosmology : str
        Cosmology model to use
    """
    
    cosmo = set_cosmo(cosmology)
    try:
        tage = float(eval(tage)) * u.Myr
        z = z_at_value(cosmo.age, tage)
        print(f"tage = {tage}, z = {z:.4e}")
    except (ValueError, TypeError, NameError):
        raise CLIError("Invalid age of the Universe. Please provide a valid age of the Universe in Myr.")    

@app.command
def time_between_z(
    z1: Annotated[str, Parameter(name=['-z1', '--redshift1'])],
    z2: Annotated[str, Parameter(name=['-z2', '--redshift2'])],
    cosmology: Annotated[str, Parameter(name=['-c', '--cosmology'])] = "Planck18"
):
    """
    Calculate the time (in Myr) between two redshifts.
    
    Parameters
    ----------
    z1 : str
        First redshift value
    z2 : str
        Second redshift value
    cosmology : str
        Cosmology model to use
    
    """
    cosmo = set_cosmo(cosmology)
    try:
        z1 = float(eval(z1)) * cu.redshift
        z2 = float(eval(z2)) * cu.redshift
        dt = np.abs(cosmo.age(z1).to(u.Myr) - cosmo.age(z2).to(u.Myr))
        print(f"z1 = {z1}, z2 = {z2}, dt = {dt:.4e}")
    except (ValueError, TypeError, NameError):
        raise CLIError("Invalid redshifts. Please provide valid redshifts.")

@app.command
def H_at_z(
    z: Annotated[str, Parameter(name=['-z', '--redshift'])],
    cosmology: Annotated[str, Parameter(name=['-c', '--cosmology'])] = "Planck18"
):
    """
    Calculate the Hubble parameter (in km/s/Mpc) at a given redshift.
    
    Parameters
    ----------
    z : str
        Redshift value
    cosmology : str
        Cosmology model to use
    """
    cosmo = set_cosmo(cosmology)
    try: 
        z = float(eval(z)) * cu.redshift
        Hz = cosmo.H(z)
        print(f"z = {z}, H(z) = {Hz:.4e}")
    except (ValueError, TypeError, NameError):
        raise CLIError("Invalid redshift. Please provide a valid redshift.")

@app.command
def comoving_length_to_deg(
    z: Annotated[str, Parameter(name=['-z', '--redshift'])],
    x: Annotated[str, Parameter(name=['-x', '--comoving_distance'])],
    cosmology: Annotated[str, Parameter(name=['-c', '--cosmology'])] = "Planck18"
):
    """
    Convert comoving length (in Mpc) at a given redshift to angular size (in deg).
    
    Parameters
    ----------
    z : str
        Redshift value
    x : str
        Comoving distance in Mpc
    cosmology : str
        Cosmology model to use
    """
    cosmo = set_cosmo(cosmology)
    try:
        z = float(eval(z)) * cu.redshift
        x = float(eval(x)) * u.Mpc
        theta = cosmo.arcsec_per_kpc_comoving(z) * x.to(u.kpc)
        theta = theta.to(u.deg)
        print(f"z = {z}, x = {x}, theta = {theta:.4e}")
    except (ValueError, TypeError, NameError):
        raise CLIError("Invalid redshift or comoving distance. Please provide a valid redshift and comoving distance.")
 
@app.command
def deg_to_comoving_length(
    z: Annotated[str, Parameter(name=['-z', '--redshift'])],
    theta: Annotated[str, Parameter(name=['-t', '--angular_size'])],
    cosmology: Annotated[str, Parameter(name=['-c', '--cosmology'])] = "Planck18"
):
    """
    Convert angular size (in deg) at a given redshift to comoving length (in Mpc).
    
    Parameters
    ----------
    z : str
        Redshift value
    theta : str
        Angular size in deg
    cosmology : str
        Cosmology model to use
    """
    cosmo = set_cosmo(cosmology)
    try:
        z = float(eval(z)) * cu.redshift
        theta = float(eval(theta)) * u.deg
        x = theta.to(u.arcsec) / cosmo.arcsec_per_kpc_comoving(z)
        x = x.to(u.Mpc)
        print(f"z = {z}, theta = {theta:.4e}, x = {x:.4e}")
    except (ValueError, TypeError, NameError):
        raise CLIError("Invalid redshift or angular size. Please provide a valid redshift and angular size.")

@app.command
def deg2_to_comoving_area(
    z: Annotated[str, Parameter(name=['-z', '--redshift'])],
    omega: Annotated[str, Parameter(name=['-o', '--angular_area'])],
    cosmology: Annotated[str, Parameter(name=['-c', '--cosmology'])] = "Planck18"
):
    """
    Convert angular area (in deg^2) to comoving area (in Mpc^2) at a given redshift.
    
    Parameters
    ----------
    z : str
        Redshift value
    omega : str
        Angular area in deg^2
    cosmology : str
        Cosmology model to use
    
    """
    cosmo = set_cosmo(cosmology)
    try:
        z = float(eval(z)) * cu.redshift
        omega = float(eval(omega)) * u.deg**2
        A = omega.to(u.arcsec**2) / cosmo.arcsec_per_kpc_comoving(z)**2 
        A = A.to(u.Mpc**2)
        print(f"z = {z}, omega = {omega:.4e}, A = {A:.4e}")
    except (ValueError, TypeError, NameError):
        raise CLIError("Invalid redshift or angular area. Please provide a valid redshift and angular area.")

@app.command
def deg2_to_comoving_volume(
    z1: Annotated[str, Parameter(name=['-z1', '--redshift1'])],
    z2: Annotated[str, Parameter(name=['-z2', '--redshift2'])],
    omega: Annotated[str, Parameter(name=['-o', '--angular_area'])],
    cosmology: Annotated[str, Parameter(name=['-c', '--cosmology'])] = "Planck18"
):
    """
    Convert angular area (in deg^2) between two redshifts to comoving volume (in Mpc^3).
    
    Parameters
    ----------
    z1 : str
        First redshift value
    z2 : str
        Second redshift value
    omega : str
        Angular area in deg^2
    cosmology : str
        Cosmology model to use
    """
    cosmo = set_cosmo(cosmology)
    try:
        z1 = float(eval(z1)) * cu.redshift
        z2 = float(eval(z2)) * cu.redshift
        omega = float(eval(omega)) * u.deg**2
        V = omega.to(u.steradian)/(4*np.pi*u.steradian) * abs(cosmo.comoving_volume(z1) - cosmo.comoving_volume(z2))
        print(f"z1 = {z1}, z2 = {z2}, omega = {omega:.4e}, V = {V:.4e}")
    except (ValueError, TypeError, NameError):
        raise CLIError("Invalid redshifts or angular area. Please provide valid redshifts and angular area.")

@app.command
def comoving_to_proper(
    z: Annotated[str, Parameter(name=['-z', '--redshift'])],
    x: Annotated[str, Parameter(name=['-x', '--comoving_distance'])],
    cosmology: Annotated[str, Parameter(name=['-c', '--cosmology'])] = "Planck18"
):
    """
    Convert comoving distance (in cMpc) at a given redshift to proper distance (in pMpc).
    
    Parameters
    ----------
    z : str
        Redshift value
    x : str
        Comoving distance in cMpc
    cosmology : str
        Cosmology model to use
    """
    cosmo = set_cosmo(cosmology)
    try:
        z = float(eval(z)) * cu.redshift
        x = float(eval(x))
        a = cosmo.scale_factor(z)
        y = x * a
        print(f"z = {z}, x = {x:.4e} cMpc --> y = {y:.4e} pMpc")
    except (ValueError, TypeError, NameError):
        raise CLIError("Invalid redshift or comoving distance. Please provide a valid redshift and comoving distance.")

@app.command
def proper_to_comoving(
    z: Annotated[str, Parameter(name=['-z', '--redshift'])],
    y: Annotated[str, Parameter(name=['-y', '--proper_distance'])],
    cosmology: Annotated[str, Parameter(name=['-c', '--cosmology'])] = "Planck18"
):
    """
    Convert proper distance (in pMpc) at a given redshift to comoving distance (in cMpc).
    
    Parameters
    ----------
    z : str
        Redshift value
    y : str
        Proper distance in pMpc
    cosmology : str
        Cosmology model to use
    """
    cosmo = set_cosmo(cosmology)
    try:
        z = float(eval(z)) * cu.redshift
        y = float(eval(y))
        a = cosmo.scale_factor(z)
        x = y / a
        print(f"z = {z}, y = {y:.4e} pMpc --> x = {x:.4e} cMpc")
    except (ValueError, TypeError, NameError):
        raise CLIError("Invalid redshift or proper distance. Please provide a valid redshift and proper distance.")

@app.command
def photon_unit_conversion(
    photon_unit_str: Annotated[str, Parameter(name=['-p', '--photon_unit_conversion'])]
):
    """
    Convert between different photon units.
    
    Parameters
    ----------
    photon_unit_str : str
        Photon unit conversion string (e.g. "21 cm to MHz", spaces required)
    
    """
    inputs = photon_unit_str.split(' ')
    try:
        value = float(eval(inputs[0]))
        from_unit = u.Unit(inputs[1])
        to_unit = u.Unit(inputs[3])
        quantity = value * from_unit
        converted_quantity = quantity.to(to_unit, equivalencies=u.spectral())
        print(f"{quantity} = {converted_quantity:.4e}")
    except (ValueError, TypeError, IndexError):
        raise CLIError('Invalid photon unit conversion string. Please provide a valid conversion string, e.g. "21 cm to MHz" (spaces required).')
    
if __name__ == '__main__':
    try:
        app()
    except CLIError as e:
        console.print(e.render())