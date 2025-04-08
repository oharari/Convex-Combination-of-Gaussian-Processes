# Convex Combination of Gaussian Processes for Bayesian Analysis of Deterministic Computer Experiments

This repository contains R scripts and data sets to help the reader follow the ideas and examples appearing in [our 2014 Technometrics paper](https://www.tandfonline.com/doi/abs/10.1080/00401706.2013.861629)(see [preprint](https://drive.google.com/file/d/19rzOfv4Zwhyv-cIdWWo09OwXmszR-TNV/view?usp=drive_link)).

Hereinafter are the different subfolders and their contents:

## 1D Codes and Designs 

- 1D Combined GP Simulation Designs.txt : the designs used for the 1D simulation.

  an 8X100 matrix, each row of which is a 1D LHD.


- 1D Combined GP Public.R : the code for the 1D Combined GP, using the Matern kernel.


- 1D Combined GP Two Families Public.R : a Bayesian convex combination of one Gaussian 
  process based on the Matern kernel and another one based on the Cubic Spline kernel.

## 2D Codes and Designs 

- Training Designs: a folder containing the 100 size 14 LHDs used for the 2D 
  simulations whose results appear in the paper.

- 2D Combined GP Isotropic Public.R : the code for the isotropic 2D Combined GP model.

- 2D Combined GP Isotropic Advanced.R : the same code, this time with the option to 
  choose te hyperparameters for the inverse-gamma priors by maximizing the marginalized 
  likelihood over a finite set of possibilities.

- 2D Combined GP Anisotropic Public : the code for the 2D anisotropic Combined GP model.

- hyperpars.matrix.txt : a 60X4 matrix, each row of which is a quadraplet of possible 
  hyperparameters for the two inverse-gamma priors. Read from 
  2D Combined GP Isotropic Advanced.R

- maximin 14 pts.txt : a size 14 (JMP generated) maximin LHD, used in Section 5 of the paper.

- maximin 100 pts.txt : a size 100 (JMP generated) maximin LHD, appears in the discussion.

## Batch Sequential ME Designs

- Batch Sequential ME Design.R : the code for Section 7 in the paper. takes a first batch size
  14 ME design, draws a random sample from the posterior, re-estimates the paramters (be it
  the posterior mean or MAP), plugs in the estimates and chooses a second batch of inputs. the
  user may then proceed to fit the Combined GP model.

- Initial ME Design.txt: a size 14 Maximum Entropy design corresponding to the prior means,
  p = 0.5, theta_1 = 1 and theta_2 = 4. Can be recreated should the user enable the
  relevant line in the code.	

- maximin 21 pts.txt : The size 21 maximin LHD used for the results of Table 5 in the paper.

- All_Subdesigns.txt : 1000 size 7 second batch designs, each corresponding to a single triplet 
  of parameters drawn from the posterior after the first batch was sampled, serialized into a 
  7000X2 matrix.
  
- k-medoids ME Design.txt : a size 21 design, the last 7 inputs of which are the result of a
  7-medoids clustering on the data in All_Subdesigns.txt. Appears in Figure 8 in the paper 
  (left panel).

- Plug-in ME 14 plus 7 Design.txt : the design appearing in the right panel of Figure 8.
  
## Ground Vibrations Emulator

- Combined GP Ground Vibrations.R : The code for fitting an isotropic Combined GP model to the
  9 dimensional data from the ground vibrations computer experiment described in

	"Development of Algorithm for the Evaluation of Ground Vibrations in Israel due to 
         Earthquakes" by David M. Steinberg, Sigal Levi, Zeev Somer and Gideon Leonard (2007).

- Training Sets: a folder containing 9 size 50 random samples from the original data and 9 
  size 90 samples.

- Test Sets: a folder containing the remaining data (to a total of 200) for each of the 
  training sets.

## Heat Exchanger Emulator 

- Combined GP Heat Exchanger.R : the code for fitting an isotropic Combined GP model to the 4
  dimensional data from the heat exchanger computer experiment described in

  "Building surrogate models with detailed and approximate simulations" by P. Z. G. Qian, 
   C. C. Seepersad, V. R. Joseph J. K. Allen and C. F. J. Wu (2006).

- hyperpars.matrix.txt : A 624X4 matrix, each row of which is an optional quadraplet of 
  hyperparameters for the inverse-gamma priors. The brave user can enable the relevant code
  in charge of selecting the best quadraplet based on the marginalized likelihood.

 - Qian Training Set.txt : the training set for the original experiment, standardized to lay
   within the unit 4 dimensional cube.

- Qian Test Set.txt : the test set for the original experiment.
