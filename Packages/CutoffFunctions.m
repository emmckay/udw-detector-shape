(* ::Package:: *)

BeginPackage["CutoffFunctions`"]

Begin["`Private`"]

gaussianCutoff[k_,\[CurlyEpsilon]_]:=E^(-Abs[k]^2/(2 \[CurlyEpsilon]^2));
lorentzianCutoff[k_,\[CurlyEpsilon]_]:=\[CurlyEpsilon]^2/(Abs[k]^2+\[CurlyEpsilon]^2);
expCutoff[k_,\[CurlyEpsilon]_]:=E^(-Abs[k]/(2 \[CurlyEpsilon]));
sharpCutoff[k_,\[CurlyEpsilon]_]:=HeavisideTheta[k+\[CurlyEpsilon],-k+\[CurlyEpsilon]];

End[]
EndPackage[]
