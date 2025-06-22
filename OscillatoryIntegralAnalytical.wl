(* ::Package:: *)

(*Author: Dr.Maseim Bassis Kenmoe, 21.06.2025 Regensburg*)

BeginPackage["OscillatoryIntegralAnalytical`", {"DerivativeOfModulusOfPCF`"}]; Null

OscillatoryIntegralAnalytical::usage = "OscillatoryIntegralAnalytical[Order_] computes the symbolic k-th order of the Oscillatory multiple Integral..";

Begin["`Private`"]; Null


OscillatoryIntegralAnalytical[Order_Integer]:=Module[{Result}, 

factor[k_]:=(-1)^(k+1)/(2^(3k-1) k!);

OscillatoryIntegral[k_] := factor[k]DevLandauZener[k];

DevLandauZener[k_]:=Evaluate[Sum[Binomial[k,n](-(\[Pi]/2))^(k-n) (-((2(k-n))/\[Pi]))DerivativeOfModulusOfPCF[n],{n,0,k}
                              ]
                       ]/.{im->Im, re->Re};

Result=OscillatoryIntegral[Order][[1]]; 
Result

]

End[]

EndPackage[]
