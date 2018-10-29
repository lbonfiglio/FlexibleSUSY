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

CreateDiagramsOutputFileName[Rule[{initial_}, finalState_List]] :=
    "diagrams_" <> ToGenericFieldString[initial] <> StringJoin[Sort[ToGenericFieldString /@ finalState]] <> ".m";

CreateAmplitudesOutputFileName[Rule[{initial_}, finalState_List]] :=
    "amplitudes_" <> ToGenericFieldString[initial] <> StringJoin[Sort[ToGenericFieldString /@ finalState]] <> ".m";

CreateAmplitudesExprsOutputFileName[Rule[{initial_}, finalState_List]] :=
    "amplitudes_exprs_" <> ToGenericFieldString[initial] <> StringJoin[Sort[ToGenericFieldString /@ finalState]] <> ".m";

ReadAmplitudeExprs[process_, dataDir_] :=
    Module[{diagramsFile, amplitudesFile, amplitudesExprsFile,
            diagrams, amplitudes, amplitudeExprs},
           diagramsFile = FileNameJoin[{dataDir, CreateDiagramsOutputFileName[process]}];
           amplitudesFile = FileNameJoin[{dataDir, CreateAmplitudesOutputFileName[process]}];
           amplitudesExprsFile = FileNameJoin[{dataDir, CreateAmplitudesExprsOutputFileName[process]}];

           If[FileExistsQ[diagramsFile],
              diagrams = Get[diagramsFile];,
              Print["Error: input file ", diagramsFile, " does not exist."];
              Quit[2];
             ];
           If[FileExistsQ[amplitudesFile],
              amplitudes = Get[amplitudesFile];,
              Print["Error: input file ", amplitudesFile, " does not exist."];
              Quit[2];
             ];
           If[FileExistsQ[amplitudesExprsFile],
              amplitudesExprs = Get[amplitudesExprsFile];,
              Print["Error: input file ", amplitudesExprsFile, " does not exist."];
              Quit[2];
             ];

           {diagrams, amplitudes, amplitudesExprs}
          ];

GetFormFactors[amplitudesExprs_] :=
    OneLoopDecaysUtils`ExtractFormFactors /@ amplitudesExprs;

CreateOneLoopDiagramDeclaration[diagram_] := "";

CreateOneLoopDiagramDeclarations[diagrams_TopologyList] :=
    Module[{},
           ""
          ];

CreateOneLoopDiagramDefinitions[diagrams_TopologyList, amplitudesExprs_] :=
    Module[{},
           ""
          ];

CreateDiagramEvaluators[diagrams_, formFactors_] :=
    Module[{decls, defs},
           decls = CreateOneLoopDiagramDeclarations[diagrams];
           defs = CreateOneLoopDiagramDefinitions[diagrams, formFactors];
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
    {diagrams, amplitudes, amplitudesExprs} = ReadAmplitudeExprs[process, resultsDir];
    Print["... extracting form factors"];
    formFactors = GetFormFactors[amplitudesExprs];
    Print[formFactors];
    Print["... generating evaluation functions"];
    {decls, defs} = CreateDiagramEvaluators[diagrams, formFactors];
    genericOneLoopDiagramEvaluatorDecls = genericOneLoopDiagramEvaluatorDecls <> decls;
    genericOneLoopDiagramEvaluatorDefs = genericOneLoopDiagramEvaluatorDefs <> defs;
   ];

Print["Writing output files ..."];
decayAmplitudesFiles = {{FileNameJoin[{templatesDir, "decay_amplitudes.hpp.in"}],
                         FileNameJoin[{resultsDir, "decay_amplitudes.hpp"}]},
                        {FileNameJoin[{templatesDir, "decay_amplitudes.cpp.in"}],
                         FileNameJoin[{resultsDir, "decay_amplitudes.cpp"}]}};
WriteOut`ReplaceInFiles[decayAmplitudesFiles,
                        { "@genericOneLoopDiagramEvaluatorDecls@" -> genericOneLoopDiagramEvaluatorDecls,
                          "@genericOneLoopDiagramEvaluatorDefs@"  -> genericOneLoopDiagramEvaluatorDefs
                        }];

Print["Generating C++ code finished"];
