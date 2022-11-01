<script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=default"></script>



#<p style="text-align: center;">Manual</p>


An Stokes inversion code based on the Bayes' Theorem.
The original version is uesd in the paper.
[MCMC Inversion of Stokes Profiles](https://ui.adsabs.harvard.edu/abs/2019ApJ...875..127L/citations). It is updated with a different parallel method. The IO of fits file is not support for the current version.

## Dependencies

* Openmpi
* cfitsio (No need for the current version)


## Compilation

* Compile the MCMC_INVERSION 

```
cd /src
make hard=y
```
 



## Run the MCMC_INVERSION 

```
cd /Outputs/run 
mpirun -np 4 ../MCMC_INVERSION 
```
By default, there should be an input file in the running folder.

* input.dat. configuration for the MCMC_INVERSION.

OR 

The path to the input file could also be specified like

```
mpirun -np 4 ../MCMC_INVERSION ./input.dat
```

## input file

The input.inv file is used to configure the run of the code.
The keywords are as follows (only)


>* path2profile  = "path\_to\_file"
>>* "path\_to\_file": Path to the Stokes profile file.   
> 
>* lambda = 
>>* float: specify the wavelength of the line; 
>
>* geffect = 
>>* float: the effective Landu factor.
>
>* gemc_simulation = 
>>* YES: do the GEMC simulation; 
>>* NO or Others: .
>
>* gemc_nchains = 
>>* integer: Chain number for the GEMC simualtion. 
>
>* dream_nchains = 
>>* integer: Chain number for the DREAM simualtion. 
>
>* dream_ngener = 
>>* integer: the number of the generation for the DREAM simualtion. 
>
>* max_pair =
>>* integer: max number of the chain pairs for sampling. 
>
>* ncr = 
>>* integer: number of the CR values. 
>
>* cr_update = 
>>* YES: Update the crossover probabilities; 
>>* NO or Others: 
>
>* distribution_update =
>>* YES: Update the targe distribution; 
>>* NO or Others: .
>
>* 1sigma =
>>* YES: Compute the 1 sigma region; 
>>* NO or Others: .
>
>* output_samples =
>>* YES: output the samples to a binay file; 
>>* NO or Others: .

## Profile file

At the moment, we only read the profile from a .dat. 
see the examples.

###Layout of the .dat file

```
num 
Lambda I Q U V
Lambda I Q U V
.
.
Lambda I Q U V
Lambda I Q U V
```
>* num is the number of wavelength grids.
>

## Output file 

###The Samples.bin
A binary file which saves the samples.
A python subroutine (in python_tool) is provide to read the file.

####keys
>* 'Num_Par', number of the paramter.
>* 'Nsamples', number of the samples for each parameter.
>* 'Samples', an array saves the samples.
