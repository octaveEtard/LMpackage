# LMpackage
Octave Etard ( octave.etard11@imperial.ac.uk )

Matlab routines for fitting linear forward / backward models. **Requires Matlab R2019b or newer** (tested on R2019b & R2020a).


## Introduction
This package consists of Matlab routines for fitting linear forward / backward models (deconvolution). This implementation strives for computational and memory efficiency to allow for more flexibility e.g. when fitting large models with high sampling rate, or to form generic cross-validated models faster.

It incorporates high-level "wrapper" functions that should make fitting models straightforward, with users needing only to specify functions to read-in the input data.

Please note that the documentation is incomplete at the moment. There are however examples illustrating how to use the code.


## Quick start

### Installation
Add the `functions` folder and subfolders to your path. All function names start with `LM_` so as to reduce the risk of shadowing any of your own functions. The examples require some [EEGLAB](https://sccn.ucsd.edu/eeglab/index.php) functions to load the data.

### Examples
The `example` folder contains detailed examples illustrating how to use this package. These include:

1. In `examples/synthetic`:

  - `example_spikes_deconvolution.m` : deconvolution example that showcase the use of the low-level routines.
  - `example_forward.m` & `example_backward.m` : simple fitting of forward/backward models, using the wrapper functions so that one does not have to deal with the lower level functions mentioned above.


2. In `example/dataHugo` & `example/dataOctave`:

  - examples for fitting forward/backward models with cross-validation (subject-specific or generic models, depending on the examples) on two different real EEG datasets.
