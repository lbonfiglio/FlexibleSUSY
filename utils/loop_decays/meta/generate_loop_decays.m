status = 0;

baseDir = DirectoryName[$InputFileName];
fsDir = ParentDirectory[baseDir, 3];

AppendTo[$Path, FileNameJoin[{fsDir, "meta"}]];
AppendTo[$Path, baseDir];

loadFSStatus = Needs["WriteOut`"];
If[loadFSStatus === $Failed,
   Quit[1];
  ];

loadUtilsStatus = Needs["OneLoopDecaysUtils`"];
If[loadUtilsStatus === $Failed,
   Quit[1];
  ];

ToFieldString[field_[indices__Integer]] := ToString[field] <> StringJoin[ToString /@ List[indices]];
ToFieldString[field_] := ToString[field];

ToGenericFieldString[field_[indices__Integer]] := ToFieldString[field];
ToGenericFieldString[field_] := ToFieldString[field];

CreateProcessName[Rule[{initial_}, finalState_List]] :=
    ToGenericFieldString[initial] <> StringJoin[Sort[ToGenericFieldString /@ finalState]];

CreateFormFactorsOutputFileName[process_] :=
    "form_factors_" <> CreateProcessName[process] <> ".m";

ReadFormFactorsExprs[process_, dataDir_] :=
    Module[{formFactorsFile},
           formFactorsFile = FileNameJoin[{dataDir, CreateFormFactorsOutputFileName[process]}];
           If[FileExistsQ[formFactorsFile],
              result = Get[formFactorsFile],
              Print["Error: input file ", formFactorsFile, " does not exist."];
              Quit[2];
             ]
          ];

GetInsertionsAsList[FeynmanGraph[props__][insertions__]] := List[insertions];
GetInsertionsAsList[insertions_List] := insertions;

GetExternalMomentumCType[psq_] := CConversion`CreateCType[CConversion`ScalarType[realScalarCType]];
GetLoopMassCType[msq_] := CConversion`CreateCType[CConversion`ScalarType[realScalarCType]];
GetCouplingCType[cp_] := "const " <> CConversion`CreateCType[CConversion`ScalarType[complexScalarCType]] <> "&";

ExternalMomentumSymbol[idx_] := k[idx];

GetExternalMomenta[process_, diagram_] :=
    Module[{nExt},
           nExt = Length[Join[process[[1]], process[[2]]]];
           Table[ExternalMomentumSymbol[i], {i, 1, nExt}]
          ];

GetLoopMasses[process_, diagram_] :=
    Module[{formFactors, masses},
           formFactors = Last[diagram];
           masses = Cases[formFactors, Mass[field_, Loop], {0, Infinity}];
           DeleteDuplicates[masses]
          ];

GetCouplings[process_, diagram_] :=
    Module[{formFactors},
           formFactors = Last[diagram];
           couplings = Join[Cases[formFactors, SARAH`Cp[fields__], {0, Infinity}],
                            Cases[formFactors, SARAH`Cp[fields__][lor__], {0, Infinity}]];
           DeleteDuplicates[couplings]
          ];

CreateGraphIDString[GraphID[Topology == t_, Generic == g_, Number == n_]] :=
    "t" <> ToString[t] <> "g" <> ToString[g] <> "n" <> ToString[n];

(* @todo fix an appropriate naming convention *)
CreateOneLoopDiagramName[process_, graphID_, insertions_] :=
    Module[{processName, idString, nExt, internalFields, internalFieldsOrder},
           processName = CreateProcessName[process];
           idString = CreateGraphIDString[graphID];
           nExt = Length[Join[process[[1]], process[[2]]]];
           internalFields = DeleteCases[GetInsertionsAsList[insertions], (Field[i_] -> field_) /; i <= nExt];
           internalFieldsOrder = Ordering[internalFields /. (Field[i_] -> f_) :> i];
           internalFieldsLabel = StringJoin[ToGenericFieldString /@ ((#[[2]]& /@ internalFields)[[internalFieldsOrder]])];
           "diagram_" <> processName <> "_" <> idString <> "_" <> internalFieldsLabel
          ];

CreateOneLoopDiagramDeclaration[process_, diagram_] :=
    Module[{returnType, diagramName, externalMomenta, loopMasses, couplings, args = ""},
           returnType = CConversion`CreateCType[CConversion`ScalarType[complexScalarCType]];
           diagramName = CreateOneLoopDiagramName[process, diagram[[1]], diagram[[3]]];
           externalMomenta = GetExternalMomenta[process, diagram];
           loopMasses = GetLoopMasses[process, diagram];
           couplings = GetCouplings[process, diagram];
           args = StringJoin[Riffle[Join[GetExternalMomentumCType /@ externalMomenta,
                                         GetLoopMassCType /@ loopMasses,
                                         GetCouplingCType /@ couplings], ", "]];
           returnType <> " calculate_" <> diagramName <> "(" <> args <> ");"
          ];

CreateOneLoopDiagramDeclarations[process_, diagrams_] :=
    StringJoin[Riffle[CreateOneLoopDiagramDeclaration[process, #]& /@ diagrams, "\n\n"]];

CreateOneLoopDiagramDefinitions[diagrams_TopologyList, amplitudesExprs_] :=
    Module[{},
           ""
          ];

CreateDiagramEvaluators[process_, diagrams_] :=
    Module[{decls, defs},
           decls = CreateOneLoopDiagramDeclarations[process, diagrams];
           defs = CreateOneLoopDiagramDefinitions[process, diagrams];
           {decls, defs}
          ];

genericProcesses = {
    {S} -> {S, S} (*,
    {S} -> {F, F},
    {S} -> {V, V},
    {S} -> {S, V} *)
};

genericOneLoopDiagramEvaluatorDecls = "";
genericOneLoopDiagramEvaluatorDefs = "";

For[i = 1, i <= Length[genericProcesses], i++,
    process = genericProcesses[[i]];
    Print["Creating C++ code for process: ", process, " ..."];
    Print["... reading input files"];
    diagramExprs = ReadFormFactorsExprs[process, resultsDir];
    Print["... generating evaluation functions"];
    {decls, defs} = CreateDiagramEvaluators[process, diagramExprs];
    genericOneLoopDiagramEvaluatorDecls = genericOneLoopDiagramEvaluatorDecls <> decls;
    genericOneLoopDiagramEvaluatorDefs = genericOneLoopDiagramEvaluatorDefs <> defs;
   ];

Print["Writing output files ..."];
decayAmplitudesFiles = {{FileNameJoin[{templatesDir, "one_loop_decay_diagrams.hpp.in"}],
                         FileNameJoin[{resultsDir, "one_loop_decay_diagrams.hpp"}]},
                        {FileNameJoin[{templatesDir, "one_loop_decay_diagrams.cpp.in"}],
                         FileNameJoin[{resultsDir, "one_loop_decay_diagrams.cpp"}]}};
WriteOut`ReplaceInFiles[decayAmplitudesFiles,
                        { "@genericOneLoopDiagramEvaluatorDecls@" -> TextFormatting`WrapLines[genericOneLoopDiagramEvaluatorDecls],
                          "@genericOneLoopDiagramEvaluatorDefs@"  -> genericOneLoopDiagramEvaluatorDefs
                        }];

Print["Generating C++ code finished"];
