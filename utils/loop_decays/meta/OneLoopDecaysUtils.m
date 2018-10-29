BeginPackage["OneLoopDecaysUtils`", {"CConversion`", "TextFormatting`"}];

GetProcessType::usage="";
ExtractFormFactors::usage="";
ToFSConventions::usage="";

CreateAmplitudeFunctionDecl::usage="";
CreateAmplitudeFunctionDef::usage="";

DefaultStandardMatrixElements::nodefault="No default basis defined for process ``";
ExtractFormFactors::missterms="Could not associate the following terms with a matrix element: ``";

Begin["`Private`"];

GetGenericFieldSymbol[{field_[indices__], properties__}] := field;
GetGenericFieldSymbol[{field_, properties__}] := field;

GetGenericProcess[Rule[inFields_List, outFields_List]] :=
    Rule[GetGenericFieldSymbol /@ inFields, GetGenericFieldSymbol /@ outFields];

GetProcessType[Amp[process_][expr_]] := GetProcessType[process];

IsSSSDecay[process_] := GetProcessType[process] === Rule[{S}, {S, S}];
IsSFFDecay[process_] := GetProcessType[process] === Rule[{S}, {F, F}];
IsSVVDecay[process_] := GetProcessType[process] === Rule[{S}, {V, V}];
IsSSVDecay[process_] := Or[GetProcessType[process] === Rule[{S}, {S, V}],
                           GetProcessType[process] === Rule[{S}, {V, S}]];
(*
DefaultStandardMatrixElements[process_?IsSSSDecay] := {1};
DefaultStandardMatrixElements[process_?IsSFFDecay] :=
DefaultStandardMatrixElements[process_?IsSSVDecay] :=
DefaultStandardMatrixElements[process_?IsSVVDecay] :=
*)
DefaultStandardMatrixElements[process_] :=
    (
     Message[DefaultStandardMatrixElements::nodefault, GetProcessType[process]];
     Abort[];
    )

IsLoopFunction[A0[__]] := True;
IsLoopFunction[A00[__]] := True;
IsLoopFunction[A0i[__]] := True;
IsLoopFunction[B0[__]] := True;
IsLoopFunction[DB0[__]] := True;
IsLoopFunction[B00[__]] := True;
IsLoopFunction[DB00[__]] := True;
IsLoopFunction[B001[__]] := True;
IsLoopFunction[B0i[__]] := True;
IsLoopFunction[B1[__]] := True;
IsLoopFunction[DB1[__]] := True;
IsLoopFunction[B11[__]] := True;
IsLoopFunction[DB11[__]] := True;
IsLoopFunction[B111[__]] := True;
IsLoopFunction[C0[__]] := True;
IsLoopFunction[C0i[__]] := True;
IsLoopFunction[D0[__]] := True;
IsLoopFunction[D0i[__]] := True;
IsLoopFunction[E0[__]] := True;
IsLoopFunction[E0i[__]] := True;
IsLoopFunction[F0[__]] := True;
IsLoopFunction[F0i[__]] := True;
IsLoopFunction[f_] := False;

IsNonTrivialLorentzStructure[Pair[__]] := True;
IsNonTrivialLorentzStructure[Eps[__]] := True;
IsNonTrivialLorentzStructure[DiracChain[__]] := True;
IsNonTrivialLorentzStructure[WeylChain[__]] := True;
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

CollectLorentzStructures[Amp[process_][expr_]] :=
    Module[{couplingRules, loopFuncRules, compactExpr, terms, result},
           couplingRules = Rule[#, Unique["coup"]]& /@ DeleteDuplicates[Cases[expr, G[_][_][__][__]]];
           loopFuncRules = Rule[#, Unique["lpFn"]]& /@ DeleteDuplicates[Cases[expr, f_[args__] /; IsLoopFunction[f[args]]]];
           compactExpr = expr /. couplingRules /. loopFuncRules;
           terms = ToTermList[Expand[compactExpr]];
           result = CollectCoefficients[JoinStructureLists @@ (FactorOutLorentzStructure /@ terms)];
           List @@ (result /. (Reverse /@ Join[couplingRules, loopFuncRules]))
          ];

ExtractFormFactors[Amp[process_][expr_]] :=
    CollectLorentzStructures[Amp[process][expr]];

CouplingToSARAHCpRules[] :=
    {
     RuleDelayed[G[_][0][fields__][1], SARAH`Cp[fields][1]],
     RuleDelayed[G[_][0][fields__][NonCommutative[Global`ChiralityProjector[1]]], SARAH`Cp[fields][SARAH`PR]],
     RuleDelayed[G[_][0][fields__][NonCommutative[Global`ChiralityProjector[-1]]], SARAH`Cp[fields][SARAH`PL]],
     RuleDelayed[G[_][0][fields__][FormCalc`Private`ga[6]], SARAH`Cp[fields][SARAH`PR]],
     RuleDelayed[G[_][0][fields__][FormCalc`Private`ga[7]], SARAH`Cp[fields][SARAH`PL]],
     RuleDelayed[G[_][0][fields__][Global`MetricTensor[KI1[i1_], KI1[i2_]]], SARAH`Cp[fields][SARAH`g[i1, i2]]],
     RuleDelayed[G[_][0][fields__]["d_"[KI1[i1_], KI1[i2_]]], SARAH`Cp[fields][SARAH`g[i1, i2]]],
     RuleDelayed[G[_][0][fields__][Mom[i1_] - Mom[i2_]], SARAH`Cp[fields][SARAH`Mom[{fields}[[i1]]] - SARAH`Mom[{fields}[[i2]]]]],
     RuleDelayed[G[_][0][fields__][NonCommutativeMultiply[KI1[i1_], FormCalc`Private`ga[6]]], SARAH`Cp[fields][SARAH`LorentzProduct[SARAH`gamma[i1], SARAH`PR]]],
     RuleDelayed[G[_][0][fields__][NonCommutativeMultiply[KI1[i1_], FormCalc`Private`ga[7]]], SARAH`Cp[fields][SARAH`LorentzProduct[SARAH`gamma[i1], SARAH`PL]]]
    };

ToSARAHCouplings[expr_] := expr /. CouplingToSARAHCpRules[];

LoopFunctionToFSNotationRules[] :=
    {
     RuleDelayed[B0i[bb0, args__], SARAH`B0[args]],
     RuleDelayed[B0i[bb00, args__], SARAH`B00[args]],
     RuleDelayed[B0i[bb1, args__], SARAH`B1[args]],
     RuleDelayed[B0i[bb11, args__], SARAH`B11[args]],
     RuleDelayed[C0i[cc0, args__], C0[args]],
     RuleDelayed[C0i[cc1, args__], SARAH`C1[args]],
     RuleDelayed[C0i[cc2, args__], SARAH`C2[args]]
    };

ToFSLoopFunctions[expr_] := expr /. LoopFunctionToFSNotationRules[];

ToFSConventions[expr_] := ToFSLoopFunctions[ToSARAHCouplings[expr]];

CreateSavedLoopFunction[loopFunction_, name_:"", loopFunctionsNamespace_:"passarino_veltman"] :=
    Module[{saveName},
           If[name != "",
              saveName = name;,
              saveName = CreateLocalLoopFunctionName[loopFunction];
             ];
           "const auto " <> saveName <> " = " loopFunctionsNamespace <> "::" <>
           CConversion`RValueToCFormString[loopFunction] <> ";"
          ];

CreateSavedLoopFunctions[loopFunctions_List] :=
    Utils`StringJoinWithSeparator[CreateSavedLoopFunction /@ loopFunctions, "\n"];

CreateAmplitudeFunctionDecl[] :=
    Module[{returnType, name}
           returnType = CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarType]];
           name = CreateDiagramName[insertions];
           returnType <> " " <> name <> "(" <> args <> ");"
          ];

CreateAmplitudeFunctionDef[formFactors_List] :=
    Module[{loopFunctions, returnType, saveLoopFunctions = "", name = "", body = ""},
           name = CreateDiagramName[insertions];
           loopFunctions = DeleteDuplicates[Cases[formFactors, f_[args__] /; IsLoopFunction[f[args]]]];
           saveLoopFunctions = CreateSavedLoopFunctions[loopFunctions];
           body = body <> saveLoopFunctions;
           returnType = CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarType]];
           returnType <> " " <> name <> "(" <> args <> ")\n{\n" <>
           TextFormatting`IndentText[body] <> "}\n"
          ];

End[];

EndPackage[];
