# monteLLTB

### Description
We embeded the `vd2020` code (available **[here](https://github.com/valkenburg/vd2020)**) into `montepython` to create `monteLLTB`: a cosmological solver and sampler for the $\Lambda$LTB model. Taking advantage of the likelihood and sampler structure of `montepython` we include the $\Lambda$LTB cosmology by adapting the likelihood computation scheme. We started defining the method `ini_LLTB` in `sampler.py`,
which executes the solver `vd2020` considering the current sampled point. Then, a call for `ini_LLTB` is included into the method `compute_lkl` to pass the $\Lambda$LTB solution to the corresponding likelihood. Note that this is possible since the method of the likehood `loglkl` now receives a new argument `LLTBin`, which contains the $\Lambda$LTB solution. We also modified the likelihoods in order to compute the observables according the $\Lambda$LTB predictions. Note that the output of `vd2020` is managed by the file `LLTB_functions.py`, which contains definitions of distances and metric functions. Finally, it is important to mention that we modified `vd2020` in order to customize the management of error, output precision and outputted functions. However, the core of the $\Lambda$LTB solver, the implementation to compute $R(t,r)$ through Carlson's elliptic integrals (**[Valkenburg 2011](https://arxiv.org/abs/1104.1082)**), remained unchanged.

### Prerequisites

   - CLASS: Cosmic Linear Anisotropy Solving System
   - Monte Python (version >= v3.3.0) with the Planck 20218 likelihoods
   - scikit-learn

### Installation
To install monteLLTB you should fisrt compile the `vd2020` following the instructions on `vd2020/README.PDF`. Once `vd2020` is installed you should modify `montepython` by doing:
```
cp -r montepython_files/* /path-to-your-montepython/

```
then modify the file `default.conf.template` to include the path to your `vd2020` installation (besides the path to `class` and `clik` likelihoods). Set `default.conf.template` as your default configuration.

Repeat the procedure for `class` by doing:
```
cp -r class_files/* /path-to-your-class/

```
Compile `class` again and have fun with inhomogeneous cosmology! 

### Usage 
`monteLLTB` was first introduced in (**[Camarena et al. 2021](https://arxiv.org/abs/2107.02296)**), please cite this paper if you make use of the code. 