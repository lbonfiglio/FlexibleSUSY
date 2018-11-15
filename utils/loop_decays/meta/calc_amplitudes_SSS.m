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

ToFieldString[field_[indices__Integer]] := ToString[field] <> StringJoin[ToString /@ List[indices]];
ToFieldString[field_] := ToString[field];

ToGenericFieldString[field_[indices__Integer]] := ToFieldString[field];
ToGenericFieldString[field_] := ToFieldString[field];

CreateDiagramsOutputFileName[Rule[{initial_}, finalState_List]] :=
    "diagrams_" <> ToGenericFieldString[initial] <> StringJoin[Sort[ToGenericFieldString /@ finalState]] <> ".m";

CreateAmplitudesOutputFileName[Rule[{initial_}, finalState_List]] :=
    "amplitudes_" <> ToGenericFieldString[initial] <> StringJoin[Sort[ToGenericFieldString /@ finalState]] <> ".m";

CreateAmplitudesExprsOutputFileName[Rule[{initial_}, finalState_List]] :=
    "amplitudes_exprs_" <> ToGenericFieldString[initial] <> StringJoin[Sort[ToGenericFieldString /@ finalState]] <> ".m";

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

CreateAmplitudes[diagrams_, generic_:True] :=
    Module[{amps},
           amps = FeynArts`CreateFeynAmp[diagrams];
           If[generic,
              FeynArts`PickLevel[Generic][amps],
              amps
             ]
          ];

CalculateAmplitudes[amplitudeHead_, amplitudeExpr_] :=
    (
     FormCalc`ClearProcess[];
     FormCalc`CalcFeynAmp[amplitudeHead[amplitudeExpr], OnShell -> False] //. Subexpr[] //. Abbr[] //. GenericList[]
    );

topologies = FeynArts`CreateTopologies[1, 1 -> 2, ExcludeTopologies -> Internal];

process = {S[1]} -> {S[1], S[1]};

diagramsOutputFile = CreateDiagramsOutputFileName[process];
amplitudesOutputFile = CreateAmplitudesOutputFileName[process];
amplitudesExprsOutputFile = CreateAmplitudesExprsOutputFileName[process];

classDiags = FeynArts`InsertFields[topologies, process, InsertionLevel -> Classes, Model -> "SM"];
diagsOutputStatus = WriteFeynArtsOutputFile[FileNameJoin[{resultsDir, diagramsOutputFile}], classDiags];
If[diagsOutputStatus === $Failed,
   status = 2;
  ];

amplitudes = FeynArts`CreateFeynAmp[classDiags];
ampsOutputStatus = WriteFeynArtsOutputFile[FileNameJoin[{resultsDir, amplitudesOutputFile}], amplitudes];
If[ampsOutputStatus === $Failed,
   status = 2;
  ];

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

ConvertFeynmanGraphToList[graph_] :=
    {OneLoopDecaysUtils`GetGraphNumber[graph],
     OneLoopDecaysUtils`GetGraphCombinatorialFactor[graph],
     OneLoopDecaysUtils`GetGraphInsertions[graph]
    };

CollectDiagramInfo[diagrams_, formFactors_] :=
    Module[{i, j, nTopologies = Length[diagrams],
            topology, insertions, nGenericInsertions, genericInsertions,
            genericAmps, count = 1},
           First[Last[
               Reap[
                    For[i = 1, i <= nTopologies, i++,
                        topology = diagrams[[i, 1]];
                        insertions = diagrams[[i, 2]];
                        genericInsertions = List @@ (#[[1]]& /@ insertions);
                        nGenericInsertions = Length[genericInsertions];
                        genericAmps = List @@ formFactors[[count ;; count + nGenericInsertions - 1]];
                        count += nGenericInsertions;
                        MapThread[Sow[List[topology, #1, Simplify[#2]]]&, {genericInsertions, genericAmps}];
                       ];
                   ]]]
          ];

Print["Extracting form factors ..."];
Print["expr = ", amplitudesExprs];
formFactors = OneLoopDecaysUtils`ExtractFormFactors /@ amplitudesExprs;

Print["Converting form factors ..."];
(*formFactors = OneLoopDecaysUtils`ToFSConventions /@ formFactors;*)

Print["Combining graph info ..."];
contributions = CollectDiagramInfo[classDiags, formFactors];
Print["contributions = ", contributions];
