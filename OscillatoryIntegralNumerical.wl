(* ::Package:: *)

(*Author: Dr.Maseim Bassis Kenmoe, 21.06.2025 Regensburg*)

BeginPackage["OscillatoryIntegralNumerical`"]; Null

OscillatoryIntegralNumerical::usage = "OscillatoryIntegralNumerical[Order_] computes the symbolic k-th order of the Oscillatory multiple Integral..";

Begin["`Private`"]; Null

formatD[expr_] :=
	RowBox[{SubsuperscriptBox[
		StyleBox["D",FontSlant->"Plain"],
		RowBox[{"-","1"}],
		RowBox[{"(",ToBoxes[expr],")"}]
	]}];

formatDConj[expr_] :=
	RowBox[{SubsuperscriptBox[StyleBox["D",FontSlant->"Plain"],
	RowBox[{"-","1"}],
	RowBox[{"(",ToBoxes[expr],")","*"}]
	]}];

list1[val_] := Flatten[
   Table[
	 If[a!=b,
	   {
	     modAnaPCF[a+r,0] modAnaPCFconj[b+r,0]->re[modAnaPCF[a+r,0] modAnaPCFconj[b+r,0]]+I im[modAnaPCF[a+r,0] modAnaPCFconj[b+r,0]],
	     modAnaPCF[b+r,0] modAnaPCFconj[a+r,0]->re[modAnaPCF[a+r,0] modAnaPCFconj[b+r,0]]-I im[modAnaPCF[a+r,0] modAnaPCFconj[b+r,0]]
        },
        Unevaluated[Sequence[]]
		  ],
          {a,0,val},{b,0,val}
	]
];

list2[val_] := Join[
   Table[
     modAnaPCF[a+r,0] -> ModifiedParabolicCylinderD[a,\[Nu],-I(1-I)t], 
     {a,0,val}
   ],
   Table[
     modAnaPCFconj[a+r,0] -> ModifiedParabolicCylinderDconj[a,\[Nu],I(1+I)t], 
     {a,0,val}
   ]
];

DevModParabolicCylinderD[r_, n_, \[Nu]_] := 
  Sum[
    Sum[
      Sum[
        Binomial[n, k] * P[n - k, m][\[Nu]] * Q[k, l][\[Nu]] *
        modAnaPCF[r + m, \[Nu]] * modAnaPCFconj[r + l, \[Nu]],
        {l, 0, k}
      ],
      {m, 0, n - k}
    ],
    {k, 0, n}
  ];

P[n_,m_][\[Nu]_] := Which[
        n<m, 0,
        m<0, 0,
        n==m, (-I)^n,
        n!=m,Binomial[n,m](-I)^m (TnuPower[n-m-1,g1[nu]]/.{nu->\[Nu]})
 ];

Q[n_,m_][\[Nu]_] := Which[
        n<m, 0,
        m<0, 0,
        n==m, (I)^n,
        n!=m, Binomial[n,m](I)^m (TnuPowerConj[n-m-1,g2[nu]]/.{nu->\[Nu]})
];

f1[nu_] := PolyGamma[0,-(-I nu-1)]-Log[2]/2; f2[nu_]:=PolyGamma[0,-(I nu-1)]-Log[2]/2;

g1[nu_] := -I f1[nu];  g2[nu_]:=I f2[nu];

T1[expr_] := D[expr,nu]-I f1[nu]*expr; T2[expr_]:=D[expr,nu]+I f2[nu]*expr;

TnuPower[n_Integer,expr_] := Nest[T1,expr,n];

TnuPowerConj[n_Integer,expr_] := Nest[T2,expr,n];

DevOrder[n_]:=FullSimplify[
                    Collect[
                       Expand[
                          FullSimplify[
                              DevModParabolicCylinderD[r,n,0]
                          ]
                       ],modAnaPCF[r,0] modAnaPCFconj[r,0]]
                    /.list1[n]
                  ]/.list2[n];

OscillatoryIntegralNumerical[Order_Integer]:=Module[{Result}, 

factor[k_]:=(-1)^(k+1)/(2^(3k-1) k!);

OscillatoryIntegral[k_] := factor[k]DevLandauZener[k];

DevLandauZener[k_]:=Evaluate[Sum[Binomial[k,n](-(\[Pi]/2))^(k-n) (-((2(k-n))/\[Pi]))DevOrder[n],{n,0,k}
                              ]
                       ]/.{im->Im, re->Re};

Result=OscillatoryIntegral[Order]; 
Result

]

End[]

EndPackage[]
