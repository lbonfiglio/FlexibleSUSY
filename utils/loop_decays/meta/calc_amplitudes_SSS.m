status = 0;

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

CalculateAmplitudes[amplitudes_] :=
    Module[{result},
           FormCalc`ClearProcess[];
           result = FormCalc`CalcFeynAmp[Head[amplitudes][#]]& /@ amplitudes;
           result //. Subexpr[] //. Abbr[] //. GenericList[]
          ];

topologies = FeynArts`CreateTopologies[1, 1 -> 2, ExcludeTopologies -> Internal];

process = {S[1]} -> {S[1], S[1]};

diagramsOutputFile = CreateDiagramsOutputFileName[process];
amplitudesOutputFile = CreateAmplitudesOutputFileName[process];
amplitudesExprsOutputFile = CreateAmplitudesExprsOutputFileName[process];

classDiags = FeynArts`InsertFields[topologies, process, InsertionLevel -> Classes];
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

amplitudesExprs = CalculateAmplitudes[genericAmplitudes];
exprsOutputStatus = WriteFormCalcOutputFile[FileNameJoin[{resultsDir, amplitudesExprsOutputFile}], amplitudesExprs];
If[exprsOutputStatus === $Failed,
   status = 2;
  ];