# Initial Condition generation

To reproduce the initial conditions used in this work follow these steps.

## Toycluster

Initial conditions were constructed with [`ToyCluster`](https://github.com/LudwigBoess/Toycluster) with version `sha1` key: `fdfc98c8338541a3b68c2b43079ac95883fe3d2a`.
The `Makefile` in this folder can be used for compilation. Adjust library paths as needed.

## Parameter file

Please use the parameter file `CIZA_J2242_Pink_Xe_0-5.par` in this folder.

## Running

You can construct the initial conditions by running

```
export OMP_NUM_THREADS=16
./Toycluster CIZA_J2242_Pink_Xe_0-5.par
```