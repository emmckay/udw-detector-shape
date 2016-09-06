# udw-detector-shape
all the mathematica stuff that I need to run on torch to tell me what the shape of a SC qubit looks like using the UdW model!


%%%%%%%%%%%%%%%% GENERAL %%%%%%%%%%%%%%%%%%%%%%%

code organized by switching function used, then by smearing
goal of code to calculate transition probability <g|\rho|g>


%%%%%%%%%%%%%%%%%% COSINE TRAP SWITCHING %%%%%%%%%%

cosine trapezoidal switching function

C - general smearing.nb
	sets up the smearing independent portion of the expression for <g|\rho|g> as functions using NIntegrate, to be used for each specific smearing

C - gaussian smear.nb
	sets up gaussian smearing function & does numerical integration to find <g|\rho|g> for a set of widths of the gaussian & for three UV cutoffs.

C - lorentzian smear.nb
	sets up lorentzian smearing function & does numerical integration to find <g|\rho|g> for a set of widths of the lorentzian & for three UV cutoffs.



%%%%%%%%%%%%%%%%%% DELTA SWITCHING %%%%%%%%%%

delta switching function

D - general smearing.nb
	sets up the k integral

D - gaussian smearing.nb
	sets up gaussian smearing function & does integration to find <g|\rho|g> for a set of widths of the gaussian & for three UV cutoffs.



%%%%%%%%%%%%%%%%%% GAUSSIAN SWITCHING %%%%%%%%%%

gaussian switching function

G - general smearing.nb
	sets up the k integral

G - gaussian smear.nb
	sets up gaussian smearing function & integrates to find <g|\rho|g> for a set of widths of the gaussian & for three UV cutoffs.

G - lorentzian smear.nb
	sets up lorentzian smearing function & integrates to find <g|\rho|g> for a set of widths of the lorentzian & for three UV cutoffs.

G - delta smear.nb
	sets up delta smearing function & integrates to find <g|\rho|g> & for three UV cutoffs.	


%%%%%%%%%%%%%%%%%% SHARP TRAP SWITCHING %%%%%%%%%%

sharp trapezoidal switching function

S - general smearing.nb
	sets up the smearing independent portion of the expression for <g|\rho|g> as functions using NIntegrate, to be used for each specific smearing

S - gaussian smear.nb
	sets up gaussian smearing function & does numerical integration to find <g|\rho|g> for a set of widths of the gaussian & for three UV cutoffs.

S - lorentzian smear.nb
	sets up lorentzian smearing function & does numerical integration to find <g|\rho|g> for a set of widths of the lorentzian & for three UV cutoffs.

S - delta smear.nb
	sets up delta smearing function & does numerical integration to find <g|\rho|g> & for three UV cutoffs.