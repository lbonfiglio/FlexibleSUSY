BeginPackage["OneLoopDecaysUtils`"];

GetProcessType::usage="";
ExtractFormFactors::usage="";
ToFSConventions::usage="";

GetGraphCombinatorialFactor::usage="";
GetGraphNumber::usage="";
GetGraphInsertions::usage="";

CollectLorentzStructures::missterms="Missing terms from form factor decomposition: ``";
ExtractFormFactors::missterms="Could not associate the following terms with a matrix element: ``";

Begin["`Private`"];

GetGenericFieldSymbol[{field_[indices__], properties__}] := field;
GetGenericFieldSymbol[{field_, properties__}] := field;

GetGenericProcess[Rule[inFields_List, outFields_List]] :=
    Rule[GetGenericFieldSymbol /@ inFields, GetGenericFieldSymbol /@ outFields];

GetProcessType[FormCalc`Amp[process_][expr_]] := GetProcessType[process];

GetGraphCombinatorialFactor[FeynArts`FeynmanGraph[s_, level_][insertions__]] := s;

GetGraphNumber[FeynArts`FeynmanGraph[s_, level_ == n_][insertions__]] := n;

GetGraphInsertions[FeynArts`FeynmanGraph[s_, level_][insertions__]] := List[insertions];

IsSSSDecay[process_] := GetProcessType[process] === Rule[{S}, {S, S}];
IsSFFDecay[process_] := GetProcessType[process] === Rule[{S}, {F, F}];
IsSVVDecay[process_] := GetProcessType[process] === Rule[{S}, {V, V}];
IsSSVDecay[process_] := Or[GetProcessType[process] === Rule[{S}, {S, V}],
                           GetProcessType[process] === Rule[{S}, {V, S}]];

IsLoopFunction[LoopTools`A0[__]] := True;
IsLoopFunction[LoopTools`A00[__]] := True;
IsLoopFunction[LoopTools`A0i[__]] := True;
IsLoopFunction[LoopTools`B0[__]] := True;
IsLoopFunction[LoopTools`DB0[__]] := True;
IsLoopFunction[LoopTools`B00[__]] := True;
IsLoopFunction[LoopTools`DB00[__]] := True;
IsLoopFunction[LoopTools`B001[__]] := True;
IsLoopFunction[LoopTools`B0i[__]] := True;
IsLoopFunction[LoopTools`B1[__]] := True;
IsLoopFunction[LoopTools`DB1[__]] := True;
IsLoopFunction[LoopTools`B11[__]] := True;
IsLoopFunction[LoopTools`DB11[__]] := True;
IsLoopFunction[LoopTools`B111[__]] := True;
IsLoopFunction[LoopTools`C0[__]] := True;
IsLoopFunction[LoopTools`C0i[__]] := True;
IsLoopFunction[LoopTools`D0[__]] := True;
IsLoopFunction[LoopTools`D0i[__]] := True;
IsLoopFunction[LoopTools`E0[__]] := True;
IsLoopFunction[LoopTools`E0i[__]] := True;
IsLoopFunction[LoopTools`F0[__]] := True;
IsLoopFunction[LoopTools`F0i[__]] := True;
IsLoopFunction[f_] := False;

IsNonTrivialLorentzStructure[FormCalc`Pair[FormCalc`k[i_Integer], FormCalc`k[j_Integer]]] := False;
IsNonTrivialLorentzStructure[FormCalc`Pair[__]] := True;
IsNonTrivialLorentzStructure[FormCalc`Eps[__]] := True;
IsNonTrivialLorentzStructure[FormCalc`DiracChain[__]] := True;
IsNonTrivialLorentzStructure[FormCalc`WeylChain[__]] := True;
IsNonTrivialLorentzStructure[_] := False;

JoinStructureLists[args__StructureList] := Join[args];

CollectCoefficients[StructureList[pairs__]] :=
    Module[{likeTerms, sumTerms},
           likeTerms = Gather[List[pairs], (First[#1] === First[#2])&];
           sumTerms[terms_List] :=
               Module[{lorentzStructure},
                      lorentzStructure = First[First[terms]];
                      {lorentzStructure, Simplify[Plus @@ (Last /@ terms)]}
                     ];
           StructureList @@ (sumTerms /@ likeTerms)
          ];

MultiplyStructureLists[factors__StructureList] := Distribute[Times[factors], StructureList];

FactorOutLorentzStructure[expr_Times] :=
    Module[{factors, lorentz, scalar},
           factors = List @@ expr;
           factored = FactorOutLorentzStructure /@ factors;
           MultiplyStructureLists @@ factored
          ]

IsPositiveInteger[n_] := IntegerQ[n] && n > 0;

FactorOutLorentzStructure[Power[expr_, exponent_?IsPositiveInteger]] :=
    Module[{i},
           MultiplyStructureLists @@ Table[FactorOutLorentzStructure[expr], {i, 1, exponent}]
          ];

FactorOutLorentzStructure[expr_Plus] :=
   CollectCoefficients[JoinStructureLists @@ (FactorOutLorentzStructure /@ (List @@ expr))];

FactorOutLorentzStructure[expr_?IsNonTrivialLorentzStructure] := StructureList[{expr, 1}];

FactorOutLorentzStructure[expr_] := StructureList[{1, expr}];

ToTermList[expr_Plus] := List @@ expr;
ToTermList[expr_Times] := List[expr];

CollectLorentzStructures[FormCalc`Amp[process_][expr_]] :=
    Module[{couplingRules, loopFuncRules, compactExpr, terms, result, difference},
           couplingRules = Rule[#, Unique["coup"]]& /@ DeleteDuplicates[Cases[expr, G[_][_][__][__]]];
           loopFuncRules = Rule[#, Unique["lpFn"]]& /@ DeleteDuplicates[Cases[expr, f_[args__] /; IsLoopFunction[f[args]]]];
           compactExpr = expr /. couplingRules /. loopFuncRules;
           terms = ToTermList[Expand[compactExpr]];
           result = CollectCoefficients[JoinStructureLists @@ (FactorOutLorentzStructure /@ terms)];
           result = List @@ (result /. (Reverse /@ Join[couplingRules, loopFuncRules]));
           (* Check result is consistent with original expression *)
           difference = Simplify[expr - Plus @@ ((Times @@ #)& /@ result)];
           If[difference =!= 0,
              Message[CollectLorentzStructures::missterms, difference];
             ];
           result
          ];

ExtractFormFactors[FormCalc`Amp[process_][expr_]] :=
    CollectLorentzStructures[FormCalc`Amp[process][expr]];

CouplingToSARAHCpRules[] :=
    {
     RuleDelayed[FeynArts`G[_][0][fields__][1], SARAH`Cp[fields][1]],
     RuleDelayed[FeynArts`G[_][0][fields__][FeynArts`NonCommutative[Global`ChiralityProjector[1]]], SARAH`Cp[fields][SARAH`PR]],
     RuleDelayed[FeynArts`G[_][0][fields__][FeynArts`NonCommutative[Global`ChiralityProjector[-1]]], SARAH`Cp[fields][SARAH`PL]],
     RuleDelayed[FeynArts`G[_][0][fields__][FormCalc`Private`ga[6]], SARAH`Cp[fields][SARAH`PR]],
     RuleDelayed[FeynArts`G[_][0][fields__][FormCalc`Private`ga[7]], SARAH`Cp[fields][SARAH`PL]],
     RuleDelayed[FeynArts`G[_][0][fields__][Global`MetricTensor[KI1[i1_], KI1[i2_]]], SARAH`Cp[fields][SARAH`g[i1, i2]]],
     RuleDelayed[FeynArts`G[_][0][fields__]["d_"[KI1[i1_], KI1[i2_]]], SARAH`Cp[fields][SARAH`g[i1, i2]]],
     RuleDelayed[FeynArts`G[_][0][fields__][FeynArts`Mom[i1_] - FeynArts`Mom[i2_]], SARAH`Cp[fields][SARAH`Mom[{fields}[[i1]]] - SARAH`Mom[{fields}[[i2]]]]],
     RuleDelayed[FeynArts`G[_][0][fields__][NonCommutativeMultiply[KI1[i1_], FormCalc`Private`ga[6]]], SARAH`Cp[fields][SARAH`LorentzProduct[SARAH`gamma[i1], SARAH`PR]]],
     RuleDelayed[FeynArts`G[_][0][fields__][NonCommutativeMultiply[KI1[i1_], FormCalc`Private`ga[7]]], SARAH`Cp[fields][SARAH`LorentzProduct[SARAH`gamma[i1], SARAH`PL]]]
    };

ToSARAHCouplings[expr_] := expr /. CouplingToSARAHCpRules[];

LoopFunctionToFSNotationRules[] :=
    {
     RuleDelayed[LoopTools`B0i[LoopTools`bb0, args__], SARAH`B0[args]],
     RuleDelayed[LoopTools`B0i[LoopTools`bb00, args__], SARAH`B00[args]],
     RuleDelayed[LoopTools`B0i[LoopTools`bb1, args__], SARAH`B1[args]],
     RuleDelayed[LoopTools`B0i[LoopTools`bb11, args__], SARAH`B11[args]],
     RuleDelayed[LoopTools`C0i[LoopTools`cc0, args__], SARAH`C0[args]],
     RuleDelayed[LoopTools`C0i[LoopTools`cc1, args__], SARAH`C1[args]],
     RuleDelayed[LoopTools`C0i[LoopTools`cc2, args__], SARAH`C2[args]]
    };

ToFSLoopFunctions[expr_] := expr /. LoopFunctionToFSNotationRules[];

ToFSConventions[expr_] := ToFSLoopFunctions[ToSARAHCouplings[expr]];

End[];

EndPackage[];
