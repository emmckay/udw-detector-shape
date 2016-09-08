(* ::Package:: *)

Import["/home/emmckay/udw-detector-shape/Packages/Initialize.m"]
Import["/home/emmckay/udw-detector-shape/Packages/CutoffFunctions.m"]

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

table = ParallelTable[
	\[Alpha] + \[Lambda]^2 * NIntegrate[ CutoffFunctions`lorentzianCutoff[k,\[CurlyEpsilon]]*kIntegrand[k,r,T], {k,0,\[Infinity]},
		MaxRecursion->maxRec,
		PrecisionGoal->precGoal,
		WorkingPrecision->workPrec],
	{T,TMin,TMax,TStepSize},{r,rMin,rMax,rStepSize}]

(*save data*)

cutoffModel = "Lorentzian";
metaData = {{"alpha","Omega","cutoff scale","coupling","cutoff model"},
		{\[Alpha], \[CapitalOmega], \[CurlyEpsilon], \[Lambda], cutoffModel}};

date = DateString["ISODateTime"];
identifier = FileNameJoin[{"/home/emmckay/udw-detector-shape/Plotting/TorchOutput",cutoffModel,date}];

dataFile = StringJoin[{identifier,".csv"}];
metaFile = StringJoin[{identifier,"_metadata.csv"}];

Export[metaFile,metaData];
Export[dataFile,table];
