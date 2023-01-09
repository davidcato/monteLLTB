# VoidDistances2020

### Installation/Compilation
Compile `vd2020` following:

```
cd EllipticIntegrals/src
make

cd ../../LLTBBackground/src
make lltbpackage

cd ../../wLLTBBackground/source
make all

cd ../../src
make

cd ../wLLTBBackground/bin

```

Test using:

```
./vd2020 b 0.7 0.0245 0.1225 0.0 0.7 -1.0 -0.25 0.75 1065.

```


### Usage
Binary `vd2020`, placed at `wLLTBBackground/bin/`, has two actions:

   - b  (basic/fastest computation, does not compute dipole needed to compute Y Compton & kSZ effects)
   - c  (complete/slowest computation, include dipole and should be used only if Y Compton & kSZ are needed) you should also pass YHe z_reio in the input line
   
The 'c' computation is significant slower than the 'b' action since it also offers the dipole as function of redshift including then computation of several non-trivial integrals at differents redshift. 
   
'b' action can be used by:

```
vd2020 b H0 om_b om_dm Om_k Om_DE w_DE delta0 z_Boundary z_drag

```

while 'c' action can be used through:

```
vd2020 b H0 om_b om_dm Om_k Om_DE w_DE delta0 z_Boundary z_drag YHe z_reio

```
For pedagogical purposes we have include a Jupyter notebook available in the folder `notebooks`. 
