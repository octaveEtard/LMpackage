# LMpackage
Octave Etard (octave.etard11@imperial.ac.uk)

Matlab routines for fitting linear forward / backward models / CCA.

**Requires Matlab R2019b or newer** (tested on R2019b & R2020a/b).


## Introduction
This package consists of Matlab routines for fitting linear forward / backward models (deconvolution). It implements ridge-regularised linear models, with options to use different L2 penalties (e.g. curvature). This implementation strives for computational and memory efficiency to allow for more flexibility e.g. when fitting large models with high sampling rate, or to form generic cross-validated models faster.

It incorporates high-level "wrapper" functions that should make fitting models straightforward, with users needing only to specify functions to read-in the input data, and contains examples illustrating how to use the code.


## Quick start

### Installation
Add the `functions` folder to your path. The code is structured as a [Matlab package](https://uk.mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html) contained in `functions/+LM` to reduce the risk of shadowing any of the user's own functions. Functions can be called by using the `LM.` prefix (e.g. `out = LM.someFunction(x,y,z)`). Alternatively the required functions can also be [imported](https://uk.mathworks.com/help/matlab/matlab_oop/importing-classes.html) (e.g. `import LM.someFunction; out = someFunction(x,y,z)`).
Note a large number of calls of the type `LM.someFunction` in a **script** can lead to decreased performance. This can be fixed by encapsulating the calls in a function or by importing the relevant functions.

### Examples
The `example` folder contains detailed examples illustrating how to use this package based on synthetic data. These include:

  - `example_spikes_deconvolution.m` : deconvolution example that showcase the use of the low-level routines.
  - `example_forward.m` & `example_backward.m` : simple fitting of forward/backward models, using the wrapper functions so that one does not have to deal with the lower level functions mentioned above.

For more examples, have a look at [LMpackage_examples](https://github.com/octaveEtard/LMpackage_examples), which demonstrates the use of this package to fit forward / backward models with cross-validation (subject-specific and generic (population) models) on real EEG datasets.
