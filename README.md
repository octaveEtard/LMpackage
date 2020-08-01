# LMpackage
Octave Etard ( octave.etard11@imperial.ac.uk )

Matlab routines for fitting linear forward / backward models.

**Requires Matlab R2019b or newer** (tested on R2019b & R2020a).


## Introduction
This package consists of Matlab routines for fitting linear forward / backward models (deconvolution). It implements ridge-regularised linear models, with options to use different L2 penalties (e.g. curvature). This implementation strives for computational and memory efficiency to allow for more flexibility e.g. when fitting large models with high sampling rate, or to form generic cross-validated models faster.

It incorporates high-level "wrapper" functions that should make fitting models straightforward, with users needing only to specify functions to read-in the input data.

Please note that the documentation is incomplete at the moment. There are however examples illustrating how to use the code.


## Quick start

### Installation
Add the `functions` (and `functions/demo` if you wish to run the examples) folder(s) to your path. All function names start with `LM_` so as to reduce the risk of shadowing any of your own functions.

### Examples
The `example` folder contains detailed examples illustrating how to use this package based on synthetic data. These include:

  - `example_spikes_deconvolution.m` : deconvolution example that showcase the use of the low-level routines.
  - `example_forward.m` & `example_backward.m` : simple fitting of forward/backward models, using the wrapper functions so that one does not have to deal with the lower level functions mentioned above.

For more examples, have a look at [LMpackage_examples](https://github.com/octaveEtard/LMpackage_examples), which demonstrates the use of this package to fit forward / backward models with cross-validation (subject-specific and generic (population) models) on real EEG datasets.
