---------------------------------------------------
Simulator for Quantized MU-MIMO-OFDM Uplink
---------------------------------------------------
(c) 2016 Christoph Studer (studer@cornell.edu)
---------------------------------------------------

# Important information:

If you are using the simulator (or parts of it) for a publication, then you MUST cite our paper:

@article{StuderDurisi:2016,
  author={C. Studer and G. Durisi},
  journal={IEEE Transactions on Communications},
  title={Quantized Massive {MU-MIMO-OFDM} Uplink},
  year={2016},
  volume={64},
  number={6},
  pages={2387-2399},
  month={June}
}

and clearly mention this in your paper. If your code also makes use of the FASTA toolbox, then you MUST cite

@article{GoldsteinStuderBaraniuk:2014,
  Author = {Goldstein, Tom and Studer, Christoph and Baraniuk, Richard},
  Title = {A Field Guide to Forward-Backward Splitting with a {FASTA} Implementation},
  year = {2014},
  journal = {arXiv eprint},
  volume = {abs/1411.3406},
  url = {http://arxiv.org/abs/1411.3406},
  ee = {http://arxiv.org/abs/1411.3406}
}

@misc{FASTA:2014,
  Author = {Goldstein, Tom and Studer, Christoph and Baraniuk, Richard},
  title = {{FASTA}:  A Generalized Implementation of Forward-Backward Splitting},
  note = {http://arxiv.org/abs/1501.04979},
  month = {January},
  year = {2015}
}

See here for more details about FASTA: https://www.cs.umd.edu/~tomg/FASTA.html

If you are thinking of contacting us, please do not e-mail the author to ask for download instructions, installation guidelines, or the toolbox itself. Note that we will NOT help to debug user-generated code not included in the provided package. If, however, you notice a bug in our code, please be so kind to contact C. Studer (studer@cornell.edu). Note that some of the modes have been modified when working on the paper and not all parameter scripts will run right out-of-the-box. 

The package is supplied "as is", without any accompanying support services, maintenance, or future updates. We make no warranties, explicit or implicit, that the software contained in this package is free of error or that it will meet your requirements for any particular application. It should not be relied on for any purpose where incorrect results could result in loss of property, personal injury, liability or whatsoever. If you do use our software for any such purpose, it is at your own risk. The authors disclaim all liability of any kind, either direct or consequential, resulting from your use of these programs.

# How to start a simulation:

Before you run any simulation or generate any plots, run the allpath.m script which will include the necessary subfolders. As a first step, I recommend to run the err_plot.m script in the plots folder. This script generates all trade-off plots from the paper. The simulation results used for that are in the results folder (which is quite large: 395MB).

Starting a simulation is also quite straightforward: The simulator uses so-called parameter files (found in the param/ folder), which define all necessary simulation settings and also start the MIMO-OFDM simulations. For example, type

>> ERR_16x8_16QAM_Tap4_Q6_MMSESoft(0)

which starts a simulation in a quantized 16 BS antenna 8 user massive MU-MIMO system with 16-QAM modulation, 3 quantization bit, and the quantized MMSE detector. The argument 0 determines the random seed. Note that we highly recommend you to execute the code step-by-step (using Matlab's debug mode) in order to gain understanding of the simulator. Since the code is quite big, it may take you a few hours to do so :-)

# Version 1.00 (Sep 3, 2016) - studer@cornell.edu, initial version for public access
