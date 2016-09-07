(* ::Package:: *)

BeginPackage["Initialize`"]

kIntegrand::usage="the integrand"

Begin["`Private`"]

(*smearing scale*)
\[Sigma]0=10^-4/\[CapitalOmega];    \[Sigma]=\[Sigma]0;    \[Sigma]Norm=\[Sigma]/\[Sigma]0;

(*smearing function and fourier transform*)
f[x_]:=1/(\[Sigma]*\[Sqrt]\[Pi]) E^(-x^2/\[Sigma]^2)
ftil[k_]:=FourierTransform[f[x],x,k,FourierParameters->{1,-1},Assumptions->\[Sigma]\[Element]Reals&&\[Sigma]>0] 

(*switching function*)
chiT1[t_,r_,T_]:=1/(r+T) 1/2*(1+Cos[\[Pi]/r (t+T/2)])
chiT2[t_,r_,T_]:=1/(r+T)
chiT3[t_,r_,T_]:=1/(r+T) 1/2*(1+Cos[\[Pi]/r (t-T/2)])

(*all those integrands*)
tpInt1[t_?NumericQ,k_?NumericQ,r_?NumericQ,T_?NumericQ]:=tpInt1[t,k,r,T]=NIntegrate[chiT1[tp,r,T]*Cos[k*(t-tp)]*Cos[\[CapitalOmega]*(t-tp)],{tp,-T/2-r,t}]
tInt1[k_?NumericQ,r_?NumericQ,T_?NumericQ]:=tInt1[k,r,T]=\[Alpha]*NIntegrate[chiT1[t,r,T] tpInt1[t,k,r,T],
{t,-T/2-r,-T/2}]
tpInt2[t_,k_,r_,T_]:=Integrate[chiT2[tp,r,T]*Cos[k*(t-tp)]*Cos[\[CapitalOmega]*(t-tp)],{tp,-T/2,t}]
tInt2[k_,r_,T_]:=\[Alpha]*Integrate[chiT2[t,r,T] tpInt2[t,k,r,T],{t,-T/2,T/2}]
tpInt3[t_?NumericQ,k_?NumericQ,r_?NumericQ,T_?NumericQ]:=tpInt3[t,k,r,T]=NIntegrate[chiT3[tp,r,T]*Cos[k*(t-tp)]*Cos[\[CapitalOmega]*(t-tp)],{tp,T,t}]
tInt3[k_?NumericQ,r_?NumericQ,T_?NumericQ]:=tInt3[k,r,T]=\[Alpha]*NIntegrate[chiT3[t,r,T] tpInt3[t,k,r,T],{t,T/2,T/2+r}]

(*the integrand we really want*)
kIntegrand[k_,r_,T_]:=-(1/\[Pi])*k*ftil[k]^2*(tInt1[k,r,T]+tInt2[k,r,T]+tInt3[k,r,T])

End[]
EndPackage[]
