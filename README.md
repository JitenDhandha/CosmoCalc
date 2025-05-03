# CosmoCalc

This is a simple command line tool for printing the cosmological parameters or calculating the following cosmological quantities:

* Hubble parameter at a given redshift
* Age of the universe at a given redshift
* Redshift at a given age of the universe
* Time between two redshifts
* Comoving length to angular degree
* Angular degree to comoving length
* Angular area to comoving area
* Angular area to comoving volume
* Comoving length to proper length at a given redshift
* Proper length to comoving length at a given redshift

It is written in Python and uses the `astropy` library for cosmological calculations. The cosmologies supported currently are: WMAP7, WMAP9, Planck13, Planck15, Planck18.

## Installation and usage

To install CosmoCalc, clone this repository in a directory of your choice and add the following line to your startup bash script (e.g. `.bashrc` or `.bash_profile` or `.zshrc`):

```
alias cosmocalc="python {directory}/CosmoCalc/cosmocalc.py"
```

Now, from anywhere in your terminal, you can run the command `cosmocalc` to calculate cosmological quantities, such as:

```
cosmocalc --function "z_to_tage" --redshift 20 --cosmology Planck18
cosmocalc -f "z_to_tage" -z 20 -c Planck18
```

Both are equivalent and will calculate the age of the universe at redshift 20 using the Planck 2018 cosmology. You can also use the `--help` flag to see all available options. Instead of supplying numbers, you can also supply strings with math operations, such as:

```
cosmocalc -f "time_between_z" -z "1080-80/2" -z2 "1080+80/2"
```

which is roughly the duration of the Recombination event!