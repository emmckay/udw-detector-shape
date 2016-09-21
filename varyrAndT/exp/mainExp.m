(* ::Package:: *)

(*Import["/home/emmckay/udw-detector-shape/Packages/Initialize.m"]
Import["/home/emmckay/udw-detector-shape/Packages/CutoffFunctions.m"]*)

(* configuration for starting remote kernels *)
Needs["SubKernels`RemoteKernels`"]

(* initialize the kernels on all machines defined in the host file *)
hosts=Import["machines_tmp","List"]
math = "/opt/apps/software/Mathematica/10.2.0/Executables/WolframKernel" <>
   " -wstp -linkmode Connect `4` -linkname `2` -subkernel -noinit >&  /dev/null &";

(* on the master node initialize only one kernel less, since one is already running *)
imin=2;
imax=Length[hosts];
idelta=1;

Do[  Print["starting Kernel: ",i," on ",hosts[[i]]];
     LaunchKernels[RemoteMachine[hosts[[i]], "/usr/bin/ssh " <> $UserName <> "@" <> hosts[[i]] <> " \"" <> math <> "\"", i]];,
     {i,imin,imax,idelta}
]


gaussianCutoff[k_,\[CurlyEpsilon]_]:=E^(-Abs[k]^2/(2 \[CurlyEpsilon]^2));
lorentzianCutoff[k_,\[CurlyEpsilon]_]:=\[CurlyEpsilon]^2/(Abs[k]^2+\[CurlyEpsilon]^2);
expCutoff[k_,\[CurlyEpsilon]_]:=E^(-Abs[k]/(2 \[CurlyEpsilon]));
sharpCutoff[k_,\[CurlyEpsilon]_]:=HeavisideTheta[k+\[CurlyEpsilon],-k+\[CurlyEpsilon]];


(*smearing function and fourier transform*)
f[x_]:=1/(\[Sigma]*\[Sqrt]\[Pi]) E^(-x^2/\[Sigma]^2)
ftil[k_]:=FourierTransform[f[x],x,k,FourierParameters->{1,-1},Assumptions->\[Sigma]\[Element]Reals&&\[Sigma]>0] 

(*switching function*)
chiT1[t_,r_,T_]:=1/(r+T) 1/2*(1+Cos[\[Pi]/r (t+T/2)])
chiT2[t_,r_,T_]:=1/(r+T)
chiT3[t_,r_,T_]:=1/(r+T) 1/2*(1+Cos[\[Pi]/r (t-T/2)])

(*all those integrands*)
(*tpInt1[t_?NumericQ,k_?NumericQ,r_?NumericQ,T_?NumericQ]:= tpInt1[t,k,r,T]=   NIntegrate[chiT1[tp,r,T]* Cos[k*(t-tp)]*Cos[\[CapitalOmega]*(t-tp)],{tp,-T/2-r,t}, MaxRecursion->maxRec, PrecisionGoal->precGoal, WorkingPrecision->workPrec];
tInt1[ k_?NumericQ,r_?NumericQ,T_?NumericQ]:=             tInt1[k,r,T]=    \[Alpha]*NIntegrate[chiT1[t,r,T]*  tpInt1[t,k,r,T],{t,-T/2-r,-T/2},           MaxRecursion->maxRec, PrecisionGoal->precGoal, WorkingPrecision->workPrec];
tpInt2[t_?NumericQ,k_?NumericQ,r_?NumericQ,T_?NumericQ]:= tpInt2[t,k,r,T]=   NIntegrate[chiT2[tp,r,T]* Cos[k*(t-tp)]*Cos[\[CapitalOmega]*(t-tp)],{tp,-T/2,t},   MaxRecursion->maxRec, PrecisionGoal->precGoal, WorkingPrecision->workPrec];
tInt2[ k_?NumericQ,r_?NumericQ,T_?NumericQ]:=             tInt2[k,r,T]=    \[Alpha]*NIntegrate[chiT2[t,r,T]*  tpInt2[t,k,r,T],{t,-T/2,T/2},              MaxRecursion->maxRec, PrecisionGoal->precGoal, WorkingPrecision->workPrec];
tpInt3[t_?NumericQ,k_?NumericQ,r_?NumericQ,T_?NumericQ]:= tpInt3[t,k,r,T]=   NIntegrate[chiT3[tp,r,T]* Cos[k*(t-tp)]*Cos[\[CapitalOmega]*(t-tp)],{tp,T,t},      MaxRecursion->maxRec, PrecisionGoal->precGoal, WorkingPrecision->workPrec];
tInt3[ k_?NumericQ,r_?NumericQ,T_?NumericQ]:=             tInt3[k,r,T]=    \[Alpha]*NIntegrate[chiT3[t,r,T]*  tpInt3[t,k,r,T],{t,T/2,T/2+r},             MaxRecursion->maxRec, PrecisionGoal->precGoal, WorkingPrecision->workPrec];*)

tpInt1[t_?NumericQ,k_?NumericQ,r_?NumericQ,T_?NumericQ]:= tpInt1[t,k,r,T]=   NIntegrate[chiT1[tp,r,T]* Cos[k*(t-tp)]*Cos[\[CapitalOmega]*(t-tp)],{tp,-T/2-r,t}];
tInt1[ k_?NumericQ,r_?NumericQ,T_?NumericQ]:=             tInt1[k,r,T]=    \[Alpha]*NIntegrate[chiT1[t,r,T]*  tpInt1[t,k,r,T],{t,-T/2-r,-T/2}];
tpInt2[t_?NumericQ,k_?NumericQ,r_?NumericQ,T_?NumericQ]:= tpInt2[t,k,r,T]=   NIntegrate[chiT2[tp,r,T]* Cos[k*(t-tp)]*Cos[\[CapitalOmega]*(t-tp)],{tp,-T/2,t}];
tInt2[ k_?NumericQ,r_?NumericQ,T_?NumericQ]:=             tInt2[k,r,T]=    \[Alpha]*NIntegrate[chiT2[t,r,T]*  tpInt2[t,k,r,T],{t,-T/2,T/2}];
tpInt3[t_?NumericQ,k_?NumericQ,r_?NumericQ,T_?NumericQ]:= tpInt3[t,k,r,T]=   NIntegrate[chiT3[tp,r,T]* Cos[k*(t-tp)]*Cos[\[CapitalOmega]*(t-tp)],{tp,T,t}];
tInt3[ k_?NumericQ,r_?NumericQ,T_?NumericQ]:=             tInt3[k,r,T]=    \[Alpha]*NIntegrate[chiT3[t,r,T]*  tpInt3[t,k,r,T],{t,T/2,T/2+r}];

(*the integrand we really want*)
kIntegrand[k_,r_,T_]:=-(1/\[Pi])*k*ftil[k]^2*(tInt1[k,r,T]+tInt2[k,r,T]+tInt3[k,r,T])


(*set parameters*)
\[Alpha] = 1;
\[CapitalOmega] = 1;
\[CurlyEpsilon] = 10 \[CapitalOmega];
\[Lambda] = 0.1;
\[Sigma] = 10^(-4)/\[CapitalOmega];
maxRec = 15;
precGoal = 20;
workPrec = 50;

numSteps = 9;
rMin = 0.2/\[CapitalOmega];  rMax = 2/\[CapitalOmega]; rStepSize = Abs[rMin-rMax]/numSteps;
rValues = Table[r,{r,rMin,rMax,rStepSize}];

TMin = 0.4/\[CapitalOmega];  TMax = 4/\[CapitalOmega]; TStepSize = Abs[TMin-TMax]/numSteps;
TValues = Table[T,{T,TMin,TMax,TStepSize}];

(*table = ParallelTable[
	\[Alpha] + \[Lambda]^2 * NIntegrate[ expCutoff[k,\[CurlyEpsilon]]*kIntegrand[k,r,T], {k,0,\[Infinity]},
		MaxRecursion->maxRec,
		PrecisionGoal->precGoal,
		WorkingPrecision->workPrec],
	{T,TMin,TMax,TStepSize},{r,rMin,rMax,rStepSize}]*)

table = Monitor[ParallelTable[
	\[Alpha] + \[Lambda]^2 * NIntegrate[ expCutoff[k,\[CurlyEpsilon]]*kIntegrand[k,r,T], {k,0,\[Infinity]}],
	{T,TMin,TMax,TStepSize}],T];

(*save data*)

cutoffModel = "Exp";
metaData = {{"alpha","Omega","cutoff scale","coupling","cutoff model"},
		{\[Alpha], \[CapitalOmega], \[CurlyEpsilon], \[Lambda], cutoffModel}};

date = DateString["ISODateTime"];
identifier = FileNameJoin[{"/home/emmckay/udw-detector-shape/Plotting/TorchOutput",cutoffModel,date}];

dataFile = StringJoin[{identifier,".csv"}];
metaFile = StringJoin[{identifier,"_metadata.csv"}];

Export[metaFile,metaData];
Export[dataFile,table];
