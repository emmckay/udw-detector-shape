(* ::Package:: *)

Import["/home/emmckay/udw-detector-shape/Packages/Initialize.m"]
Import["/home/emmckay/udw-detector-shape/Packages/CutoffFunctions.m"]

(*set parameters*)
\[Alpha] = 1;
\[CapitalOmega] = 1;
\[CurlyEpsilon] = 10 \[CapitalOmega];
\[Lambda] = 0.1;

numSteps = 9;
rMin = 0.2/\[CapitalOmega];  rMax = 2/\[CapitalOmega]; rStepSize = Abs[rMin-rMax]/numSteps;
rValues = Table[r,{r,rMin,rMax,rStepSize}];

TMin = 0.4/\[CapitalOmega];  TMax = 4/\[CapitalOmega]; TStepSize = Abs[TMin-TMax]/numSteps;
TValues = Table[T,{T,TMin,TMax,TStepSize}];

table = ParallelTable[
	\[Alpha] + \[Lambda]^2 * NIntegrate[ lorentzianCutoff[k,\[CurlyEpsilon]]*kIntegrand[k,T], {k,0,\[Infinity]},
		MaxRecursion->maxRec,
		PrecisionGoal->precGoal,
		WorkingPrecision->workPrec],
	{T,TMin,TMax,TStepSize},{r,rMin,rMax,rStepSize}]

(*save data*)

cutoffModel = "Lorentzian";
metaData = {{"alpha","Omega","cutoff scale","coupling","cutoff model"},
		{\[Alpha], \[CapitalOmega], \[CurlyEpsilon], \[Lambda], cutoffModel}}

date = DateString["ISODateTime"];
identifier = FileNameJoin[{"/home/emmckay/udw-detector-shape/Plotting/TorchOutput",cutoffModel,date}];

dataFile = StringJoin[{identifier,".csv"}];
metaFile = StringJoin[{identifier,"_metadata.csv"}];

Export[metaFile,metaData];
Export[dataFile,table];
