(* ::Package:: *)

(*this is the sharp trap switch gaussian smear business, using numerical integration for the first time*)


ftilGs[k_]:=E^(-(1/4) k^2 \[Sigma]^2)
ggintegrandG[k_]:=1/(2 \[Pi] \[Epsilon]^2) A^4 E^(-k \[CurlyEpsilon]) k ftilGs[k]^2 ((E^(-I (c+\[Epsilon]) (k+\[CapitalOmega])) (1-\[Alpha]) (-1+E^(I \[Epsilon] (k+\[CapitalOmega]))-I k \[Epsilon]-I \[Epsilon] \[CapitalOmega]) (E^(I c (k+\[CapitalOmega]))+I E^(I (c+\[Epsilon]) (k+\[CapitalOmega])) (I+\[Epsilon] (k+\[CapitalOmega]))))/(\[Epsilon]^2 (k+\[CapitalOmega])^4)+1/(\[Epsilon]^2 (k^2-\[CapitalOmega]^2)^4)(1-\[Alpha]) (-k^6 \[Epsilon]^2-2 \[CapitalOmega]^4-\[Epsilon]^2 \[CapitalOmega]^6-k^4 (-2+\[Epsilon]^2 \[CapitalOmega]^2)+k^2 (-12 \[CapitalOmega]^2+\[Epsilon]^2 \[CapitalOmega]^4)+2 (k^4+6 k^2 \[CapitalOmega]^2+\[CapitalOmega]^4) Cos[k \[Epsilon]] Cos[\[Epsilon] \[CapitalOmega]]-(6 k^4 \[Epsilon] \[CapitalOmega]+4 k^2 \[Epsilon] \[CapitalOmega]^3+2 \[Epsilon] \[CapitalOmega]^5) Cos[k \[Epsilon]] Sin[\[Epsilon] \[CapitalOmega]]+Sin[k \[Epsilon]] ((2 k^5 \[Epsilon]+4 k^3 \[Epsilon] \[CapitalOmega]^2-6 k \[Epsilon] \[CapitalOmega]^4) Cos[\[Epsilon] \[CapitalOmega]]+(8 k^3 \[CapitalOmega]+8 \[CapitalOmega]^3) Sin[\[Epsilon] \[CapitalOmega]]))-(2 \[Alpha] Sin[c (-k+\[CapitalOmega])]^2)/(-k+\[CapitalOmega])^2+((4-6 \[Alpha]) Sin[c (k+\[CapitalOmega])]^2)/(k+\[CapitalOmega])^2)


NIntegrate[ggintegrandG[k],{k,0,\[Infinity]}]
