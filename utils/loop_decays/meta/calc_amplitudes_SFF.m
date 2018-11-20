status = 0;

baseDir = DirectoryName[$InputFileName];
AppendTo[$Path, baseDir];

CheckDirectoryContains[dir_String, files_List] :=
    Module[{workDir = Directory[], result = False},
           If[DirectoryQ[dir],
              result = FileNames[files, {dir}] =!= {};
             ];
           result
          ];

ToFieldString[Times[-1, f_]] := ToFieldString[f];
ToFieldString[field_[indices__Integer]] := ToString[field] <> StringJoin[ToString /@ List[indices]];
ToFieldString[field_] := ToString[field];

ToGenericFieldString[Times[-1, f_]] := ToGenericFieldString[f];
ToGenericFieldString[field_[indices__Integer]] := ToFieldString[field];
ToGenericFieldString[field_] := ToFieldString[field];

CreateDiagramsOutputFileName[Rule[{initial_}, finalState_List]] :=
    "diagrams_" <> ToGenericFieldString[initial] <> StringJoin[Sort[ToGenericFieldString /@ finalState]] <> ".m";

CreateAmplitudesOutputFileName[Rule[{initial_}, finalState_List]] :=
    "amplitudes_" <> ToGenericFieldString[initial] <> StringJoin[Sort[ToGenericFieldString /@ finalState]] <> ".m";

CreateAmplitudesExprsOutputFileName[Rule[{initial_}, finalState_List]] :=
    "amplitudes_exprs_" <> ToGenericFieldString[initial] <> StringJoin[Sort[ToGenericFieldString /@ finalState]] <> ".m";

CreateFormFactorsOutputFileName[Rule[{initial_}, finalState_List]] :=
    "form_factors_" <> ToGenericFieldString[initial] <> StringJoin[Sort[ToGenericFieldString /@ finalState]] <> ".m";

loadFeynArts = Needs["FeynArts`"];
loadFormCalc = Needs["FormCalc`"];

If[loadFeynArts === $Failed || loadFormCalc === $Failed,
   Quit[1];
  ];

WriteFeynArtsOutputFile[fileName_, expr_] :=
    Module[{result, commentStr},
           commentStr = "(* Generated with " <> $FeynArtsVersion <>
                        " at " <> DateString[] <> " *)";
           result = Put[OutputForm[commentStr], fileName];
           If[result === $Failed,
              Return[result];
             ];
           PutAppend[expr, fileName]
          ];

WriteFormCalcOutputFile[fileName_, expr_] :=
    Module[{result, commentStr},
           commentStr = "(* Generated with " <> $FormCalcVersion <>
                        " at " <> DateString[] <> " *)";
           result = Put[OutputForm[commentStr], fileName];
           If[result === $Failed,
              Return[result];
             ];
           PutAppend[expr, fileName]
          ];

WriteFormFactorsOutputFile[fileName_, expr_] :=
    Module[{result, commentStr},
           commentStr = "(* Generated at " <> DateString[] <> " *)";
           result = Put[OutputForm[commentStr], fileName];
           If[result === $Failed,
              Return[result];
             ];
           PutAppend[expr, fileName]
          ];

CalculateAmplitudes[amplitudeHead_, amplitudeExpr_] :=
    (
     FormCalc`ClearProcess[];
     FormCalc`CalcFeynAmp[amplitudeHead[amplitudeExpr], OnShell -> True, FermionChains -> Chiral] //. Subexpr[] //. Abbr[] //. GenericList[]
    );

topologies = FeynArts`CreateTopologies[1, 1 -> 2, ExcludeTopologies -> Internal];

process = {S[1]} -> {F[4], -F[4]};

diagramsOutputFile = CreateDiagramsOutputFileName[process];
amplitudesOutputFile = CreateAmplitudesOutputFileName[process];
amplitudesExprsOutputFile = CreateAmplitudesExprsOutputFileName[process];
formFactorsOutputFile = CreateFormFactorsOutputFileName[process];

diags = FeynArts`InsertFields[topologies, process, InsertionLevel -> Generic, Model -> "SM"];
diagsOutputStatus = WriteFeynArtsOutputFile[FileNameJoin[{resultsDir, diagramsOutputFile}], diags];
If[diagsOutputStatus === $Failed,
   status = 2;
  ];

amplitudes = FeynArts`CreateFeynAmp[diags];
ampsOutputStatus = WriteFeynArtsOutputFile[FileNameJoin[{resultsDir, amplitudesOutputFile}], amplitudes];
If[ampsOutputStatus === $Failed,
   status = 2;
  ];

GetGraphID[FeynAmp[GraphID[details__], rest__]] := GraphID[details];
graphIDs = List @@ (GetGraphID /@ amplitudes);

genericAmplitudes = FeynArts`PickLevel[Generic][amplitudes];

amplitudesExprs = CalculateAmplitudes[Head[genericAmplitudes], #]& /@ genericAmplitudes;
exprsOutputStatus = WriteFormCalcOutputFile[FileNameJoin[{resultsDir, amplitudesExprsOutputFile}], amplitudesExprs];
If[exprsOutputStatus === $Failed,
   status = 2;
  ];

loadUtils = Needs["OneLoopDecaysUtils`"];
If[loadUtils === $Failed,
   Quit[1];
  ];

GetExternalLegInsertions[diagrams_] :=
    Module[{genericDiags, topologies, externalFields,
            genericInsertions, externalInsertions},
           genericDiags = FeynArts`PickLevel[Generic][diagrams];
           topologies = List @@ (First /@ genericDiags);
           externalFields = Cases[#, Propagator[Incoming|Outgoing][v1_, v2_, f_] :> f, {0, Infinity}]& /@ topologies;
           genericInsertions = List @@ ((List @@ #)& /@ Last /@ genericDiags);
           externalInsertions = MapIndexed[With[{insertionsList = #1, topIdx = First[#2]},
                                                Select[#, MemberQ[externalFields[[topIdx]], First[#]]&]& /@ insertionsList]&,
                                           genericInsertions];
           externalInsertions = DeleteDuplicates[#, (#1[[1]] === #2[[1]])&]& /@ externalInsertions;
           If[Or @@ (Length[#] != 1)& /@ externalInsertions,
              Print["Error: multiple different external field insertions for a single topology are not handled."];
              Quit[1];
             ];
           (List @@ #[[1]])& /@ externalInsertions
          ];

ExternalMomentumSymbol[i_] := k[i];

StripSign[-f_] := f;
StripSign[f_] := f;

GenericFieldType[-field_] := GenericFieldType[field];
GenericFieldType[field_[indices__]] := GenericFieldType[field];
GenericFieldType[field_] := StripSign[field];

GetSquaredMassSymbol[mass_] := mass^2;

CreateLoopFunctionArgReplacementRules[fieldIdx_, mass_] :=
    Module[{loopFnHeads, squaredMass},
           squaredMass = GetSquaredMassSymbol[mass];
           squaredMass = squaredMass /. Index[Generation, i_] :> Symbol["Gen" <> ToString[i]];
           loopFnHeads = { SARAH`B0, SARAH`C1, SARAH`C2, SARAH`C0 };
           Rule[#[x___, squaredMass, y___], #[x, Pair[ExternalMomentumSymbol[fieldIdx],
                                                      ExternalMomentumSymbol[fieldIdx]], y]]& /@ loopFnHeads
          ];

CreateLoopFunctionArgReplacementRule[Field[fieldIdx_], mass_] :=
    CreateLoopFunctionArgReplacementRule[fieldIdx, mass];

CreateSquaredMassReplacementRules[fieldIdx_, mass_] :=
    Module[{squaredMass},
           squaredMassSymbol = GetSquaredMassSymbol[mass];
           {Rule[squaredMassSymbol, Pair[ExternalMomentumSymbol[fieldIdx], ExternalMomentumSymbol[fieldIdx]]],
            Rule[mass^2, Pair[ExternalMomentumSymbol[fieldIdx], ExternalMomentumSymbol[fieldIdx]]]}
          ];

CreateSquaredMassReplacementRules[Field[fieldIdx_], mass_] :=
    CreateSquaredMassReplacementRules[fieldIdx, mass];

CreateGenericMassReplacementRules[fieldIdx_, field_, mass_] :=
    Module[{genericField},
           genericField = GenericFieldType[field][Index[Generic, fieldIdx]];
           {Rule[mass, Mass[genericField]]} /. Index[Generation, i_] :> Symbol["Gen" <> ToString[i]]
          ];

CreateGenericMassReplacementRules[Field[fieldIdx_], field_, mass_] :=
    CreateGenericMassReplacementRules[fieldIdx, field, mass];

CreateCouplingIndexReplacementRules[fieldIdx_, field_] :=
    Module[{genericField},
           genericField = GenericFieldType[field][Index[Generic, fieldIdx]];
           {Rule[SARAH`Cp[x___, field, y___], SARAH`Cp[x, genericField, y]],
            Rule[SARAH`Cp[x___, -field, y___], SARAH`Cp[x, -genericField, y]]}
          ];

CreateCouplingIndexReplacementRules[Field[fieldIdx_], field_] :=
    CreateCouplingIndexReplacementRules[fieldIdx, field];

ToGenericExpressions[externalLegInsertions_, statesInfo_, formFactors_List] :=
    Module[{incomingInfo, incomingField, incomingMass,
            outgoingInfo, outgoingFields, outgoingMasses,
            fieldMasses, loopFnArgReplacements, squaredMassReplacements, massReplacements,
            couplingReplacements, result},
           result = formFactors;
           incomingInfo = statesInfo[[1]];
           outgoingInfo = statesInfo[[2]];
           incomingField = Select[externalLegInsertions, (Last[#] === First[incomingInfo[[1]]])&];
           outgoingFields = Complement[externalLegInsertions, incomingField];
           incomingMass = (Rule[#[[1]], #[[3]]]& /@ incomingInfo) /. (Reverse /@ externalLegInsertions);
           outgoingMasses = (Rule[-#[[1]], #[[3]]]& /@ outgoingInfo) /. (Reverse /@ externalLegInsertions);
           fieldMasses = Join[incomingMass, outgoingMasses];
           loopFnArgReplacements = Flatten[CreateLoopFunctionArgReplacementRules[#[[1]], #[[2]]]& /@ fieldMasses];
           squaredMassReplacements = Flatten[CreateSquaredMassReplacementRules[#[[1]], #[[2]]]& /@ fieldMasses];
           massReplacements = Flatten[CreateGenericMassReplacementRules[#[[1]], #[[1]] /. externalLegInsertions, #[[2]]]& /@ fieldMasses];
           couplingReplacements = Flatten[CreateCouplingIndexReplacementRules[#[[1]], #[[2]]]& /@ externalLegInsertions];
           formFactors //. loopFnArgReplacements /. squaredMassReplacements /. massReplacements //. couplingReplacements
          ];

CanonicalizeCouplings[formFactors_List] :=
    Module[{countGhosts, countVectors, couplingsUUV, dupCouplingsUUV,
            couplingSubs, result},
           result = formFactors;
           countGhosts[fields_List] := Count[fields, U[__] | -U[__]];
           countVectors[fields_List] := Count[fields, V[__] | -V[__]];
           couplingsUUV = Cases[formFactors, SARAH`Cp[fields__][lor_] /; (countGhosts[List[fields]] == 2 &&
                                                                          countVectors[List[fields]] == 1), {0, Infinity}];
           dupCouplingsUUV = Select[GatherBy[couplingsUUV, #[[0]]&], (Length[#] > 1) &];
           (* Assume that only one component of the kinematic vector for a given U-U-V
              coupling has non-vanishing coefficient *)
           If[dupCouplingsUUV =!= {},
              couplingSubs = DeleteCases[#, SARAH`Cp[U[a___, Index[Generic, i_], b___],
                                                     U[c___, Index[Generic, j_], d___],
                                                     V[x___, Index[Generic, k_], y___]][SARAH`Mom[U[a___, Index[Generic, i_], b___]]]]& /@ dupCouplingsUUV;
              couplingSubs = Rule[#, 0]& /@ Flatten[couplingSubs];
              result = result /. couplingSubs;
             ];
           result
          ];

ConvertFeynmanGraphToList[graph_] :=
    {OneLoopDecaysUtils`GetGraphNumber[graph],
     OneLoopDecaysUtils`GetGraphCombinatorialFactor[graph],
     OneLoopDecaysUtils`GetGraphInsertions[graph]
    };

CollectDiagramInfo[ids_, diagrams_, formFactors_] :=
    Module[{i, j, nTopologies = Length[diagrams],
            topology, insertions, nGenericInsertions, genericInsertions,
            genericAmps, genericIDs, count = 1},
           First[Last[
               Reap[
                    For[i = 1, i <= nTopologies, i++,
                        topology = diagrams[[i, 1]];
                        insertions = diagrams[[i, 2]];
                        genericInsertions = List @@ insertions;
                        nGenericInsertions = Length[genericInsertions];
                        genericAmps = List @@ formFactors[[count ;; count + nGenericInsertions - 1]];
                        genericIDs = ids[[count ;; count + nGenericInsertions - 1]];
                        count += nGenericInsertions;
                        MapThread[Sow[List[#1, topology, #2, Simplify[#3]]]&, {genericIDs, genericInsertions, genericAmps}];
                       ];
                   ]]]
          ];

ReplaceExternalFieldInsertions[insertions_, diagram_] :=
    Module[{genericFieldReplacements},
           genericFieldReplacements = Rule[#[[2]], GenericFieldType[#[[2]]]]& /@ insertions;
           {diagram[[1]], diagram[[2]], diagram[[3]] /. genericFieldReplacements, diagram[[4]]}
          ];

Print["Extracting form factors ..."];

formFactors = OneLoopDecaysUtils`ExtractFormFactors /@ amplitudesExprs;

Print["Converting form factors ..."];

formFactors = OneLoopDecaysUtils`ToFSConventions /@ formFactors;

externalInsertions = GetExternalLegInsertions[diags];
If[Length[externalInsertions] != 1,
   Print["Error: expected only a single topology."];
   Quit[1];
  ];
externalInsertions = First[externalInsertions];
statesInfo = amplitudesExprs[[0,1,2]];
formFactors = CanonicalizeCouplings /@ ToGenericExpressions[externalInsertions, statesInfo, #]& /@ formFactors;

Print["Combining graph info ..."];

contributions = CollectDiagramInfo[graphIDs, diags, formFactors];
contributions = ReplaceExternalFieldInsertions[externalInsertions, #]& /@ contributions;

formFactorsOutputStatus = WriteFormFactorsOutputFile[FileNameJoin[{resultsDir, formFactorsOutputFile}], contributions];
If[formFactorsOutputStatus === $Failed,
   status = 3;
  ];
