# CosmoCalc

This is a simple command line tool for printing the cosmological parameters or calculating the following cosmological quantities:

* Hubble parameter at a given redshift
* Age of the universe at a given redshift
* Redshift at a given age of the universe
* Time between two redshifts
* Comoving length to angular degree at a given redshift
* Angular degree to comoving length at a given redshift
* Angular area to comoving area at a given redshift
* Angular area to comoving volume for a given redshift range
* Comoving length to proper length at a given redshift
* Proper length to comoving length at a given redshift
* Comoving distance to redshift
* Redshift to comoving distance
* Photon unit conversions

The code is written in Python and uses the `astropy` library for cosmological calculations. The cosmologies supported currently are: WMAP7, WMAP9, Planck13, Planck15, Planck18. The command line interface is built using the `cyclopts` and `rich` libraries.

## Installation and usage

1. Navigate to the directory where you want to install CosmoCalc and clone the repository:
```
git clone https://github.com/JitenDhandha/CosmoCalc.git
```

2. Make sure you have Python 3 installed along with the required packages. You can install the required packages using pip:
```
pip install numpy astropy matplotlib cyclopts
```

3. Add an alias to your terminal startup script (e.g. `.bashrc`, `.bash_profile`, or `.zshrc`) to easily run the tool from anywhere:
```
alias cosmocalc="python {directory}/CosmoCalc/cosmocalc.py"
```

4. You're all set! Now, you should be able to use CosmoCalc from your terminal. Test the installation by running:
```
cosmocalc --help
```

5. Here are some example commands you can run, with math operations supported via strings:

```
cosmocalc z_to_tage --redshift 1080 --cosmology Planck18
cosmocalc z_to_tage -z 1080 -c Planck18
cosmocalc time_between_z -z1 "1080-80/2" -z2 "1080+80/2"
```