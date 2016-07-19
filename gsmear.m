(* ::Package:: *)

(* ::Package:: *)
(*integrate this gaussian smearing & gaussian switching with polyn cutoff thing -- fingers crossed*)

ftils[k_]:=E^(-(1/4) k^2 \[Sigma]^2);
ggintegrand[k_]:=ftils[k]^2*(\[Eta]/(4\[Sqrt]\[Pi])*(1-\[Alpha])*E^(-\[Eta]^2/2 (k^2+\[CapitalOmega]^2))
-\[Eta]^2/(8\[Sqrt]2) \[Alpha]*E^(-\[Eta]^2/2 (k+\[CapitalOmega])^2)*(1+E^(2*\[Eta]*\[CapitalOmega]*k)));
cutoff[k_]:=((k+2)/(5*\[CapitalOmega]))^-2;
gg=2*Integrate[ggintegrand[k]*cutoff[k],{k,0,\[Infinity]}]
