<<<<<<< HEAD:Sgsmear2.m
(* ::Package:: *)

(*this is the sharp trap switch gaussian smear business, using numerical integration for the first time*)


ftilGs[k_]:=E^(-(1/4) k^2 \[Sigma]^2)
ggintegrandG[k_]:=1/(2 \[Pi] \[Epsilon]^2) A^4 E^(-k \[CurlyEpsilon]) k ftilGs[k]^2 ((E^(-I (c+\[Epsilon]) (k+\[CapitalOmega])) (1-\[Alpha]) (-1+E^(I \[Epsilon] (k+\[CapitalOmega]))-I k \[Epsilon]-I \[Epsilon] \[CapitalOmega]) (E^(I c (k+\[CapitalOmega]))+I E^(I (c+\[Epsilon]) (k+\[CapitalOmega])) (I+\[Epsilon] (k+\[CapitalOmega]))))/(\[Epsilon]^2 (k+\[CapitalOmega])^4)+1/(\[Epsilon]^2 (k^2-\[CapitalOmega]^2)^4)(1-\[Alpha]) (-k^6 \[Epsilon]^2-2 \[CapitalOmega]^4-\[Epsilon]^2 \[CapitalOmega]^6-k^4 (-2+\[Epsilon]^2 \[CapitalOmega]^2)+k^2 (-12 \[CapitalOmega]^2+\[Epsilon]^2 \[CapitalOmega]^4)+2 (k^4+6 k^2 \[CapitalOmega]^2+\[CapitalOmega]^4) Cos[k \[Epsilon]] Cos[\[Epsilon] \[CapitalOmega]]-(6 k^4 \[Epsilon] \[CapitalOmega]+4 k^2 \[Epsilon] \[CapitalOmega]^3+2 \[Epsilon] \[CapitalOmega]^5) Cos[k \[Epsilon]] Sin[\[Epsilon] \[CapitalOmega]]+Sin[k \[Epsilon]] ((2 k^5 \[Epsilon]+4 k^3 \[Epsilon] \[CapitalOmega]^2-6 k \[Epsilon] \[CapitalOmega]^4) Cos[\[Epsilon] \[CapitalOmega]]+(8 k^3 \[CapitalOmega]+8 \[CapitalOmega]^3) Sin[\[Epsilon] \[CapitalOmega]]))-(2 \[Alpha] Sin[c (-k+\[CapitalOmega])]^2)/(-k+\[CapitalOmega])^2+((4-6 \[Alpha]) Sin[c (k+\[CapitalOmega])]^2)/(k+\[CapitalOmega])^2)


NIntegrate[ggintegrandG[k],{k,0,\[Infinity]}]
=======
(* ::Package:: *)

(*this is the sharp trap switch gaussian smear business, using numerical integration for the first time*)
(*do a numerical integration as a function of c and \[Sigma]*)
(*this package heavily inspired by Alison's XvT code on github - pabetism/QuadraticCoupling*)
(*we have to set parameters bc we are numerically integrating*)

\[Epsilon]=0.5;  (*could be a function of c?*)
\[CapitalOmega]=1; \[Alpha]=0; \[Lambda]=.1; \[CurlyEpsilon]=.25;

numSteps=200;

\[Sigma]Min=0;
\[Sigma]Max=2;
\[Sigma]StepSize=Abs[\[Sigma]Min-\[Sigma]Max]/numSteps;

cMin=0;
cMax=2;
cStepSize=Abs[cMin-cMax]/numSteps;

ftilGs[k_,\[Sigma]_] := E^(-(1/4) k^2 \[Sigma]^2)

ggintegrandG[k_,\[Sigma]_,c_] := 1/(2 \[Pi] \[Epsilon]^2) A^4 E^(-k \[CurlyEpsilon]) k ftilGs[k,\[Sigma]]^2 ((E^(-I (c+\[Epsilon]) (k+\[CapitalOmega])) (1-\[Alpha]) (-1+E^(I \[Epsilon] (k+\[CapitalOmega]))-I k \[Epsilon]-I \[Epsilon] \[CapitalOmega]) 
	(E^(I c (k+\[CapitalOmega]))+I E^(I (c+\[Epsilon]) (k+\[CapitalOmega])) (I+\[Epsilon] (k+\[CapitalOmega]))))/(\[Epsilon]^2 (k+\[CapitalOmega])^4)+1/(\[Epsilon]^2 (k^2-\[CapitalOmega]^2)^4)(1-\[Alpha]) 
	(-k^6 \[Epsilon]^2-2 \[CapitalOmega]^4-\[Epsilon]^2 \[CapitalOmega]^6-k^4 (-2+\[Epsilon]^2 \[CapitalOmega]^2) + k^2 (-12 \[CapitalOmega]^2+\[Epsilon]^2 \[CapitalOmega]^4)+2 (k^4+6 k^2 \[CapitalOmega]^2+\[CapitalOmega]^4) Cos[k \[Epsilon]] Cos[\[Epsilon] \[CapitalOmega]]
	-(6 k^4 \[Epsilon] \[CapitalOmega]+4 k^2 \[Epsilon] \[CapitalOmega]^3+2 \[Epsilon] \[CapitalOmega]^5) Cos[k \[Epsilon]] Sin[\[Epsilon] \[CapitalOmega]]+Sin[k \[Epsilon]] ((2 k^5 \[Epsilon]+4 k^3 \[Epsilon] \[CapitalOmega]^2-6 k \[Epsilon] \[CapitalOmega]^4) Cos[\[Epsilon] \[CapitalOmega]]
	+(8 k^3 \[CapitalOmega]+8 \[CapitalOmega]^3) Sin[\[Epsilon] \[CapitalOmega]]))-(2 \[Alpha] Sin[c (-k+\[CapitalOmega])]^2)/(-k+\[CapitalOmega])^2 + ((4 - 6 \[Alpha]) Sin[c (k+\[CapitalOmega])]^2)/(k+\[CapitalOmega])^2)

table=ParallelTable[
	NIntegrate[ggintegrandG[k,\[Sigma],c],{k,0,\[Infinity]}],
	{\[Sigma],\[Sigma]Min,\[Sigma]Max,\[Sigma]StepSize},
	{c,cMin,cMax,cStepSize}]

date=StringJoin["SharpTrapGSmear",DateString["ISODateTime"]]
identifier=FileNameJoin[{"../Plotting/TorchOutput",date}];

datafile=StringJoin[{identifier,".csv"}]
metafile=StringJoin[{identifier,"_metadata.csv"}]

Print["The run identifier is ",date]

metadata={{"epsilon","Omega","alpha","lambda","curly epsilon","numSteps","\[Sigma]Min","\[Sigma]Max","\[Sigma]StepSize","cMin","cMax","cStepSize"},
	{\[Epsilon],\[CapitalOmega],\[Alpha],\[Lambda],\[CurlyEpsilon],numSteps,\[Sigma]Min,\[Sigma]Max,\[Sigma]StepSize,cMin,cMax,cStepSize}}

Export[datafile,table]
Export[metafile,metadata]
>>>>>>> 64142c2d88f546070b878f4c5088f725da3997f9:Sgsmear.m
