(* :Copyright:

   ====================================================================
   This file is part of FlexibleSUSY.

   FlexibleSUSY is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published
   by the Free Software Foundation, either version 3 of the License,
   or (at your option) any later version.

   FlexibleSUSY is distributed in the hope that it will be useful, but
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with FlexibleSUSY.  If not, see
   <http://www.gnu.org/licenses/>.
   ====================================================================

*)

BeginPackage["Decays`", {"SARAH`", "CConversion`", "CXXDiagrams`", "TreeMasses`", "TextFormatting`", "Utils`", "Vertices`"}];

FSParticleDecay::usage="head used for storing details of an particle decay,
in the format
   FSParticleDecay[particle, {final state particle}, {diagram 1, diagram 2, ...}]
";

CreateCompleteParticleList::usage="";
GetDecaysForParticle::usage="";
GetVerticesForDecays::usage="gets required vertices for a list of decays";

CallDecaysCalculationFunctions::usage="creates calls to functions calculating
decays of the given particles.";
CreateDecaysCalculationPrototypes::usage="creates prototypes for convenience
functions calculating all decays for each particle.";
CreateDecaysCalculationFunctions::usage="creates definitions for convenience
functions calculating all decays for each particle.";
CreatePartialWidthCalculationPrototypes::usage="create prototypes for
functions computing partial widths of all decays.";
CreatePartialWidthCalculationFunctions::usage="creates definitions for
functions computing partial widths of all decays.";
CreateDecaysGetterFunctions::usage="create getters for specific particle decays";
CreateDecayTableGetterPrototypes::usage="create getter prototypes for C++ decay table";
CreateDecayTableGetterFunctions::usage="create getter definitions for C++ decay table";
CreateDecayTableInitialization::usage="create C++ initializer for decay table."

Begin["`Private`"];

GetInitialState[FSParticleDecay[particle_, finalState_List, diagrams_List]] := particle;
GetFinalState[FSParticleDecay[particle_, finalState_List, diagrams_List]] := finalState;
GetDecayDiagrams[FSParticleDecay[particle_, finalState_List, diagrams_List]] := diagrams;

GetDecayTopologies[nProducts_, nLoops_] :=
    (
     Print["Error: decay topology with ", nProducts, " particles and ", nLoops, " loops not supported"];
     Quit[1]
    );

(* tree-level two-body decay, with
   vertex 1 = incoming state
   vertex 2 = outgoing state
   vertex 3 = outgoing state
   vertex 4 = internal vertex *)
GetDecayTopologies[2, 0] :=
    {
     {{0,0,0,1},
      {0,0,0,1},
      {0,0,0,1},
      {1,1,1,0}}
    };

GetDecayTopologies[2, 1] :=
    {
     {{0,0,0,1,0,0},
      {0,0,0,0,1,0},
      {0,0,0,0,0,1},
      {1,0,0,0,1,1},
      {0,1,0,1,0,1},
      {0,0,1,1,1,0}}
     ,
     {{0,0,0,1,0},
      {0,0,0,0,1},
      {0,0,0,0,1},
      {1,0,0,0,2},
      {0,1,1,2,0}}
     ,
     {{0,0,0,1,0},
      {0,0,0,1,0},
      {0,0,0,0,1},
      {1,1,0,0,2},
      {0,0,1,2,0}}
     ,
     {{0,0,0,0,1},
      {0,0,0,1,0},
      {0,0,0,0,1},
      {0,1,0,0,2},
      {1,0,1,2,0}}
    };

GetDecayTopologies[nProducts_] := Join @@ (GetDecayTopologies[nProducts, #]& /@ {0, 1});

CreateCompleteParticleList[particles_List] := DeleteDuplicates[Join[particles, SARAH`AntiField[#]& /@ particles]];

GenericScalarName[] := "scalar";
GenericVectorName[] := "vector";
GenericFermionName[] := "fermion";
GenericGhostName[] := "ghost";

GetGenericTypeName[p_?TreeMasses`IsScalar] := GenericScalarName[];
GetGenericTypeName[p_?TreeMasses`IsVector] := GenericVectorName[];
GetGenericTypeName[p_?TreeMasses`IsFermion] := GenericFermionName[];
GetGenericTypeName[p_?TreeMasses`IsGhost] := GenericGhostName[];

GetDecaysForParticle[particle_, Infinity, allowedFinalStateParticles_List] :=
    (
     Print["Error: number of final state particles must be finite."];
     Quit[1];
    )

GetDecaysForParticle[particle_, {Infinity}, allowedFinalStateParticles_List] :=
    (
     Print["Error: number of final state particles must be finite."];
     Quit[1];
    )

GetDecaysForParticle[particle_, {n_, Infinity}, allowedFinalStateParticles_List] :=
    (
     Print["Error: number of final state particles must be finite."];
     Quit[1];
    )

GetDecaysForParticle[particle_, {Infinity, n_}, allowedFinalStateParticles_List] :=
    (
     Print["Error: number of final state particles must be finite."];
     Quit[1];
    )

GetDecaysForParticle[particle_, {Infinity, Infinity}, allowedFinalStateParticles_List] :=
    (
     Print["Error: number of final state particles must be finite."];
     Quit[1];
    )

GetDecaysForParticle[particle_, maxNumberOfProducts_Integer /; maxNumberOfProducts >= 2,
                     allowedFinalStateParticles_List] :=
    GetDecaysForParticle[particle, {2, maxNumberOfProducts}, allowedFinalStateParticles];

GetDecaysForParticle[particle_, {minNumberOfProducts_Integer /; minNumberOfProducts >= 2,
                                 maxNumberOfProducts_Integer /; maxNumberOfProducts >= 2},
                                allowedFinalStateParticles_List] :=
    Module[{i, finalStateSizes},
           finalStateSizes = Table[{i}, {i, minNumberOfProducts, maxNumberOfProducts}];
           Flatten[GetDecaysForParticle[particle, #, allowedFinalStateParticles]& /@ finalStateSizes, 1]
          ];

(* returns False if state consists only of Goldstones or ghosts, for any set of generation indices *)
IsPhysicalFinalState[finalState_List] :=
    Module[{goldstones, onlyGoldstones},
           goldstones = Select[finalState, TreeMasses`IsGoldstone];
           isAlwaysGoldstone = (TreeMasses`GetDimensionWithoutGoldstones[#] == 0)& /@ goldstones;
           (Select[finalState, TreeMasses`IsGhost] === {}) && !(Or @@ isAlwaysGoldstone)
          ];

IsElectricChargeConservingDecay[initialParticle_, finalState_List] :=
    Module[{chargeSum},
           chargeSum = Simplify[Plus @@ (Join[{-TreeMasses`GetElectricCharge[initialParticle]}, TreeMasses`GetElectricCharge /@ finalState])];
           PossibleZeroQ[chargeSum]
          ];

(* @todo handle more than 2 particles in final state and non-SM color representations *)
IsColorInvariantDecay[initialParticle_, finalState_List] :=
    Module[{initialStateRep, finalStateReps, result = True},
           If[Length[finalState] == 2,
              initialStateRep = TreeMasses`GetColorRepresentation[initialParticle];
              finalStateReps = Sort[TreeMasses`GetColorRepresentation /@ finalState];
              Switch[initialStateRep,
                     S, result = ((finalStateReps === {S, S}) ||
                                  (finalStateReps === {T, T}) ||
                                  (finalStateReps === {-T, -T}) ||
                                  (finalStateReps === Sort[{-T, T}]) ||
                                  (finalStateReps === {O, O}));,
                     T|-T, result = ((finalStateReps === Sort[{T, S}]) ||
                                     (finalStateReps === Sort[{-T, S}]) ||
                                     (finalStateReps === Sort[{O, S}]));,
                     O, result = ((finalStateReps === Sort[{O, S}]) ||
                                  (finalStateReps === {T, T}) ||
                                  (finalStateReps === {-T, -T}) ||
                                  (finalStateReps === Sort[{-T, T}]));,
                     _, result = True; (* unhandled case *)
                    ];
             ];
           result
          ];

IsPossibleNonZeroDiagram[diagram_, useDependences_:False] :=
    Module[{vertices, vertexVals, isPossibleNonZeroVertex},
           vertices = CXXDiagrams`VerticesForDiagram[diagram];
           isPossibleNonZeroVertex[vertex_] := MemberQ[vertex[[2 ;;]][[All, 1]], Except[0]];
           vertexVals = SARAH`Vertex[#, UseDependences -> useDependences]& /@ vertices;
           And @@ (isPossibleNonZeroVertex /@ vertexVals)
          ];

ContainsOnlySupportedVertices[diagram_] :=
    Module[{vertices, vertexTypes, unsupportedVertices},
           vertices = CXXDiagrams`VerticesForDiagram[diagram];
           vertexTypes = Vertices`VertexTypeForFields /@ vertices;
           unsupportedVertices = Complement[vertexTypes, Vertices`VertexTypes[]];
           If[unsupportedVertices =!= {},
              MapIndexed[(If[!MemberQ[Vertices`VertexTypes[], vertexTypes[[First[#2]]]],
                             Print["Warning: vertex with fields ", #1, " is not currently supported."];
                             Print["    Diagrams involving this vertex will be discarded."];
                            ];)&, vertices];
             ];
           unsupportedVertices === {}
          ];

IsSupportedDiagram[diagram_] := ContainsOnlySupportedVertices[diagram];

GetContributingDiagramsForDecayGraph[initialField_, finalFields_List, graph_] :=
    Module[{externalFields, diagrams},
           externalFields = Join[{1 -> initialField}, MapIndexed[(First[#2] + 1 -> #1)&, finalFields]];
           diagrams = CXXDiagrams`FeynmanDiagramsOfType[graph, externalFields];
           Select[diagrams, IsPossibleNonZeroDiagram]
          ];

GetContributingGraphsForDecay[initialParticle_, finalParticles_List] :=
    Module[{nFinalParticles = Length[finalParticles], topologies, diagrams},
           topologies = GetDecayTopologies[nFinalParticles, 0];
           diagrams = Flatten[GetContributingDiagramsForDecayGraph[initialParticle, finalParticles, #]& /@ topologies, 1];
           Select[diagrams, IsSupportedDiagram]
          ];

(* defines a fixed ordering for final state particles  *)
(* @todo decide on what this ordering actually will be *)
OrderFinalState[initialParticle_, finalParticles_List] :=
    Module[{orderedFinalState},
           orderedFinalState = First[Vertices`SortCp[SARAH`Cp[Join[{initialParticle}, finalParticles]]]];
           Drop[orderedFinalState, First[Position[orderedFinalState, initialParticle]]]
          ];

GetDecaysForParticle[particle_, {exactNumberOfProducts_Integer}, allowedFinalStateParticles_List] :=
    Module[{genericFinalStates, finalStateParticlesClassified,
            isPossibleDecay, concreteFinalStates, decays},
           If[exactNumberOfProducts > 2,
              Print["Error: decays with ", exactNumberOfProducts,
                    " final particles are not currently supported."];
              Quit[1];
             ];
           genericFinalStates = GetAllowedGenericFinalStates[particle, exactNumberOfProducts];
           (* @todo checks on colour and Lorentz structure *)
           isPossibleDecay[finalState_] := (IsPhysicalFinalState[finalState] &&
                                            IsElectricChargeConservingDecay[particle, finalState] &&
                                            IsColorInvariantDecay[particle, finalState]);
           concreteFinalStates = Join @@ (GetParticleCombinationsOfType[#, allowedFinalStateParticles, isPossibleDecay]& /@ genericFinalStates);
           concreteFinalStates = OrderFinalState[particle, #] & /@ concreteFinalStates;
           decays = FSParticleDecay[particle, #, GetContributingGraphsForDecay[particle, #]]& /@ concreteFinalStates;
           Select[decays, (GetDecayDiagrams[#] =!= {})&]
          ];

GetDecaysForParticle[particle_, n_, allowedFinalStateParticles_List] :=
    (
     Print["Error: invalid number of final state particles: ", n];
     Quit[1];
    )

GatherParticlesByType[particles_List] :=
    Module[{areSameType},
           areSameType[p1_, p2_] := Or @@ ((#[p1] && #[p2])& /@ { TreeMasses`IsScalar,
                                                                  TreeMasses`IsVector,
                                                                  TreeMasses`IsFermion,
                                                                  TreeMasses`IsGhost });
           Gather[particles, areSameType]
          ];

(* returns a list of lists of the form {{particle type 1, {particles of type 1}}, {particle type 2, {particles of type 2}}, ...} *)
GetClassifiedParticleLists[particles_List] :=
    Module[{classified, foundTypes},
           classified = {GetGenericTypeName[First[#]], #}& /@ GatherParticlesByType[particles];
           foundTypes = First /@ classified;
           If[Length[Union[foundTypes]] != Length[foundTypes],
              Print["Error: particles incorrectly classified: ", classified];
              Quit[1];
             ];
           classified
          ];

BaseMulticombination[k_] := Module[{i}, Table[1, {i, 1, k}]];

(* see, e.g., http://www.martinbroadhurst.com/multicombinations.html *)
NextMulticombination[n_, combination_] :=
    Module[{k = Length[combination], i, incrementable, pos, val},
           incrementable = Position[combination, x_ /; x < n];
           If[Length[incrementable] == 0,
              {},
              pos = First[Last[incrementable]];
              val = combination[[pos]] + 1;
              Join[Take[combination, pos - 1], {val}, Table[val, {i, 1, k - pos}]]
             ]
          ];

NextMulticombinationsList[setSizes_List, combinations_List] :=
    Module[{numCombinations = Length[combinations], next},
           next = combinations;
           For[i = numCombinations, i > 0, i--,
               nextCombination = NextMulticombination[setSizes[[i]], combinations[[i]]];
               If[nextCombination =!= {},
                  next[[i]] = nextCombination;
                  Return[next];,
                  next[[i]] = BaseMulticombination[Length[next[[i]]]];
                 ];
              ];
           {}
          ];

GetParticleCombinationsOfType[genericState_List, particles_List, isValidTest_:Function[True]] :=
    Module[{genericTypeCounts, classifiedParticles, indexLists, candidate, combinations},
           genericTypeCounts = {#, Count[genericState, #]}& /@ DeleteDuplicates[genericState];
           classifiedParticles = Select[GetClassifiedParticleLists[particles], MemberQ[genericState, First[#]]&];
           genericTypeCounts = genericTypeCounts[[Flatten[Position[First /@ classifiedParticles, First[#]]& /@ genericTypeCounts]]];
           indexLists = BaseMulticombination[Last[#]]& /@ genericTypeCounts;
           combinations = Reap[
               While[indexLists =!= {},
                     candidate = Flatten[MapIndexed[With[{pos = First[#2], indices = #1},
                                                    Last[classifiedParticles[[pos]]][[#]]& /@ indices]&, indexLists]];
                     If[isValidTest[candidate],
                        Sow[candidate];
                       ];
                     indexLists = NextMulticombinationsList[Length[Last[#]]& /@ classifiedParticles, indexLists];
                    ];
               ];
           Flatten[Last[combinations], 1]
          ];

GetAllowedGenericFinalStates[particle_ /; (TreeMasses`IsScalar[particle] || TreeMasses`IsVector[particle]),
                             n_Integer] :=
    Switch[n,
           2, {{GenericScalarName[], GenericScalarName[]},
               {GenericScalarName[], GenericVectorName[]},
               {GenericVectorName[], GenericVectorName[]},
               {GenericFermionName[], GenericFermionName[]}},
           _, Print["Error: cannot determine allowed generic final states for n = ", n]; Quit[1];
          ];

GetAllowedGenericFinalStates[particle_?TreeMasses`IsFermion, n_Integer] :=
    Switch[n,
           2, {{GenericScalarName[], GenericFermionName[]},
               {GenericVectorName[], GenericFermionName[]}},
           _, Print["Error: cannot determine allowed generic final states for n = ", n]; Quit[1];
          ];

GetVerticesForDecay[decay_FSParticleDecay] :=
    Module[{diagrams = GetDecayDiagrams[decay]},
           DeleteDuplicates[Flatten[CXXDiagrams`VerticesForDiagram /@ diagrams, 1]]
          ];

GetVerticesForDecays[particleDecays_List] :=
    DeleteDuplicates[Flatten[GetVerticesForDecay /@ particleDecays, 1]]

LoopOverIndex[loopBody_String, index_, start_, stop_, type_:CConversion`ScalarType[CConversion`integerScalarCType]] :=
    Module[{idxStr, startStr, stopStr},
           idxStr = CConversion`ToValidCSymbolString[index];
           startStr = CConversion`ToValidCSymbolString[start];
           stopStr = CConversion`ToValidCSymbolString[stop];
           "for (" <> CConversion`CreateCType[type] <> " " <> idxStr <>
           " = " <> startStr <> "; " <> idxStr <> " < " <> stopStr <> "; ++" <>
           idxStr <> ") {\n" <> TextFormatting`IndentText[loopBody] <> "\n}"
          ];

(* generates a loop over the given indices, in the form
   {{idx1, start, stop}, {idx2, start, stop}, ...} with
   the first list entry being the innermost loop *)
LoopOverIndexCollection[loopBody_String, indices_List] :=
    Fold[LoopOverIndex[#1, Sequence @@ #2]&, loopBody, indices];

CreateGenericPartialWidthCalculationName[initialState_, finalState_List, fieldsNamespace_] :=
    "get_partial_width<" <> TreeMasses`CreateFieldClassName[initialState, prefixNamespace -> fieldsNamespace] <> "," <>
    Utils`StringJoinWithSeparator[TreeMasses`CreateFieldClassName[#, prefixNamespace -> fieldsNamespace]& /@ finalState, ","] <> " >";

CreatePartialWidthCalculationName[decay_FSParticleDecay, scope_:""] :=
    Module[{initialState, initialStateName,
            finalState, finalStateName},
           initialState = GetInitialState[decay];
           initialStateName = CConversion`ToValidCSymbolString[initialState];
           finalState = GetFinalState[decay];
           finalStateName = StringJoin[CConversion`ToValidCSymbolString /@ finalState];
           scope <> If[scope != "", "::", ""] <> "partial_width_" <> initialStateName <> "_to_" <> finalStateName
          ];

CreatePartialWidthCalculationPrototype[decay_FSParticleDecay] :=
    Module[{returnType = "", functionName = "", functionArgs = "",
            initialStateDim, finalStateDims},
           initialStateDim = TreeMasses`GetDimension[GetInitialState[decay]];
           finalStateDims = TreeMasses`GetDimension /@ GetFinalState[decay];
           functionArgs = "const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates&" <>
                          If[initialStateDim > 1, ", int", ""] <>
                          StringJoin[If[# > 1, ", int", ""]& /@ finalStateDims];
           functionName = CreatePartialWidthCalculationName[decay];
           returnType = CConversion`CreateCType[CConversion`ScalarType[CConversion`realScalarCType]];
           returnType <> " " <> functionName <> "(" <> functionArgs <> ") const;"
          ];

CreatePartialWidthCalculationFunction[decay_FSParticleDecay, fieldsNamespace_] :=
    Module[{i, returnType = "", functionName = "", functionArgs = "",
            initialState = GetInitialState[decay], initialStateDim,
            finalState = GetFinalState[decay], finalStateDims, setFieldIndices, body = ""},
           initialStateDim = TreeMasses`GetDimension[GetInitialState[decay]];
           finalStateDims = TreeMasses`GetDimension /@ GetFinalState[decay];
           functionArgs = "const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model" <>
                          If[initialStateDim > 1, ", int gI1", ""] <>
                          StringJoin[MapIndexed[If[#1 > 1, ", int gO" <> ToString[First[#2]], ""]&, finalStateDims]];
           functionName = CreatePartialWidthCalculationName[decay, "CLASSNAME"];
           returnType = CConversion`CreateCType[CConversion`ScalarType[CConversion`realScalarCType]];
           setFieldIndices[field_, indicesName_, indexVal_] :=
               Module[{i, dim, numIndices, result = ""},
                      dim = TreeMasses`GetDimension[field];
                      numIndices = CXXDiagrams`NumberOfFieldIndices[field];
                      result = "const field_indices<" <> TreeMasses`CreateFieldClassName[field, prefixNamespace -> fieldsNamespace] <> " >::type " <> indicesName;
                      If[numIndices == 0 || dim <= 1,
                         result = result <> "{};\n";,
                         result = result <> "{{" <> ToString[indexVal] <>
                                  StringJoin[Table[", 0", {i, 1, numIndices - 1}]] <> "}};\n";
                        ];
                      result
                     ];
           body = FlexibleSUSY`FSModelName <> "_evaluation_context context{model};\n" <>
                  StringJoin[setFieldIndices[#[[1]], #[[2]], #[[3]]]& /@
                                 Join[{{initialState, "in_indices", If[initialStateDim > 1, "gI1", ""]}},
                                      MapIndexed[{#1, "out_" <> ToString[First[#2]] <> "_indices",
                                                  If[finalStateDims[[First[#2]]] > 1, "gO" <> ToString[First[#2]], ""]}&, finalState]]];
           body = body <> "\nreturn " <> CreateGenericPartialWidthCalculationName[initialState, finalState, fieldsNamespace] <>
                  "(context, in_indices" <> StringJoin[Table[", out_" <> ToString[i] <> "_indices", {i, 1, Length[finalState]}]] <> ");\n";
           returnType <> " " <> functionName <> "(" <> functionArgs <> ") const\n{\n" <>
           TextFormatting`IndentText[body] <> "}\n"
          ];

CreatePartialWidthCalculationPrototypes[particleDecays_List] :=
    Module[{allDecays},
           allDecays = Flatten[Last /@ particleDecays];
           Utils`StringJoinWithSeparator[CreatePartialWidthCalculationPrototype /@ allDecays, "\n"]
          ];

CreatePartialWidthCalculationFunctions[particleDecays_List, fieldsNamespace_] :=
    Module[{allDecays},
           allDecays = Flatten[Last /@ particleDecays];
           Utils`StringJoinWithSeparator[CreatePartialWidthCalculationFunction[#, fieldsNamespace]& /@ allDecays, "\n"]
          ];

CallPDGCodeGetter[SARAH`bar[particle_], args__] :=
    "-" <> CallPDGCodeGetter[particle, args];

CallPDGCodeGetter[Susyno`LieGroups`conj[particle_], args__] :=
   "-" <> CallPDGCodeGetter[particle, args];

CallPDGCodeGetter[particle_, idx_String, namespace_] :=
    Module[{dim = TreeMasses`GetDimension[particle], particleStr, result = ""},
           particleStr = namespace <> If[namespace != "", "::", ""] <> CConversion`ToValidCSymbolString[particle];
           result = namespace <> If[namespace != "", "::", ""] <> "get_pdg_code_for_particle(" <>
                    particleStr;
           If[dim > 1,
              result = result <> ", " <> idx;
             ];
           result <> ")"
          ];

CallPartialWidthCalculation[decay_FSParticleDecay] :=
    Module[{i, initialState = GetInitialState[decay], initialStateDim,
            finalState = GetFinalState[decay], finalStateDims, functionArgs = "",
            pdgsList = "", loopIndices, body = ""},
           initialStateDim = TreeMasses`GetDimension[initialState];
           finalStateDims = TreeMasses`GetDimension /@ finalState;
           finalStateStarts = TreeMasses`GetDimensionStartSkippingGoldstones /@ finalState;
           functionArgs = "model_" <> If[initialStateDim > 1, ", gI1", ""] <>
                          MapIndexed[If[#1 > 1, ", gO" <> ToString[First[#2]], ""]&, finalStateDims];
           pdgsList = MapIndexed[With[{idx = First[#2]},
                                      CallPDGCodeGetter[#1, If[finalStateDims[[idx]] > 1, "gO" <> ToString[idx], ""], FlexibleSUSY`FSModelName <> "_info"]]&,
                                 finalState];
           pdgsList = "{" <> Utils`StringJoinWithSeparator[pdgsList, ", "] <> "}";
           body = "decays.set_decay(" <> CreatePartialWidthCalculationName[decay] <> "(" <> functionArgs <> "), " <> pdgsList <> ");";
           loopIndices = Reverse[Select[MapIndexed[With[{idx = First[#2]},
                                                        If[#1 > 1,
                                                           {"gO" <> ToString[idx], finalStateStarts[[idx]] - 1, #1},
                                                           {}
                                                          ]
                                                       ]&, finalStateDims], (# =!= {})&]];
           If[loopIndices =!= {},
              body = LoopOverIndexCollection[body, loopIndices];
             ];
           body <> "\n"
          ];

CreateDecaysCalculationFunctionName[particle_, scope_:""] :=
    scope <> If[scope =!= "", "::", ""] <> "calculate_" <> CConversion`ToValidCSymbolString[particle] <> "_decays";

CreateDecaysCalculationPrototype[particle_] :=
    "void " <> CreateDecaysCalculationFunctionName[particle] <> "();";

CreateDecaysCalculationPrototypes[decayLists_List] :=
    Utils`StringJoinWithSeparator[CreateDecaysCalculationPrototype[First[#]]& /@ decayLists, "\n"];

CreateLocalScope[body_] := "{\n" <> TextFormatting`IndentText[body] <> "}\n";

CreateDecaysCalculationFunction[decaysList_] :=
    Module[{particle = First[decaysList], particleDim, particleStart,
            decayChannels = Last[decaysList],
            runToScale = "", body = ""},
           particleDim = TreeMasses`GetDimension[particle];
           particleStart = TreeMasses`GetDimensionStartSkippingGoldstones[particle];
           runToScale = "const auto& decay_mass = PHYSICAL(" <>
                        CConversion`ToValidCSymbolString[TreeMasses`GetMass[particle]] <>
                        ");\nmodel_.run_to(decay_mass" <> If[particleDim > 1, "(gI1)", ""] <> ");\n";
           body = StringJoin[CallPartialWidthCalculation /@ decayChannels];
           body = "auto& decays = decay_table.get_" <> CConversion`ToValidCSymbolString[particle] <>
                  "_decays(" <> If[particleDim > 1, "gI1", ""] <> ");\n\n" <> body;
           body = "auto model_ = model;\n\nif (run_to_decay_particle_scale) {\n" <>
                  TextFormatting`IndentText[runToScale] <> "}\n\n" <> body;
           If[particleDim > 1,
              body = LoopOverIndexCollection[body, {{"gI1", particleStart - 1, particleDim}}] <> "\n";
             ];
           "void " <> CreateDecaysCalculationFunctionName[particle, "CLASSNAME"] <>
           "()\n{\n"
           <> TextFormatting`IndentText[body] <> "}\n"
          ];

CreateDecaysCalculationFunctions[particleDecays_List] :=
    Utils`StringJoinWithSeparator[CreateDecaysCalculationFunction /@ particleDecays, "\n"];

CallDecaysFunction[particle_, arg_:"model", obj_:""] :=
    obj <> CreateDecaysCalculationFunctionName[particle] <> "();\n"

CallThreadedDecaysFunction[particle_, ptr_:"this", pool_:"tp"] :=
    pool <> ".run_task([" <> ptr <> "] () { " <>
    If[ptr === "this", "", ptr <> "->"] <>
    CreateDecaysCalculationFunctionName[particle] <> "(); });\n";

CallDecaysCalculationFunctions[particles_List, enableDecaysThreads_] :=
    Module[{result = ""},
           If[enableDecaysThreads,
              result = "Thread_pool tp(std::min(std::thread::hardware_concurrency(), " <>
                       ToString[Length[particles]] <> "u));\n\n" <>
                       StringJoin[CallThreadedDecaysFunction /@ particles];
              ,
              result = StringJoin[CallDecaysFunction /@ particles];
             ];
           result
          ];

CreateDecaysGetterFunctionName[particle_] :=
    "get_" <> CConversion`ToValidCSymbolString[particle] <> "_decays";

CreateDecaysGetterFunction[particle_] :=
    Module[{dim, body = ""},
           dim = TreeMasses`GetDimension[particle];
           body = "return decay_table." <> CreateDecayTableEntryGetterName[particle] <>
                  "(" <> If[dim > 1, "i", ""] <> ");";
           "const Decays_list& " <> CreateDecaysGetterFunctionName[particle] <> "(" <>
           If[dim > 1, "int i", ""] <> ") const { " <> body <> " }"
          ];

CreateDecaysGetterFunctions[particles_List] :=
    Utils`StringJoinWithSeparator[CreateDecaysGetterFunction /@ particles, "\n"];

CreateDecayTableEntryGetterName[particle_] :=
    "get_" <> CConversion`ToValidCSymbolString[particle] <> "_decays";

CreateDecayTableEntryGetterPrototype[particle_] :=
    Module[{dim},
           dim = TreeMasses`GetDimension[particle];
           "Decays_list& " <> CreateDecayTableEntryGetterName[particle] <> "(" <>
           If[dim > 1, "int", ""] <> ");"
          ];

CreateDecayTableEntryConstGetterPrototype[particle_] :=
    Module[{dim},
           dim = TreeMasses`GetDimension[particle];
           "const Decays_list& " <> CreateDecayTableEntryGetterName[particle] <> "(" <>
           If[dim > 1, "int", ""] <> ") const;"
          ];

CreateDecayTableEntryGetterFunctionBody[particle_, rows_List] :=
    Module[{i, dim, idxName = "gI1", errMsg = "", body = ""},
           dim = TreeMasses`GetDimension[particle];
           If[dim != Length[rows],
              Print["Error: number of rows (", Length[rows], ") does not match size of"];
              Print["    ", particle, " multiplet."];
              Quit[1];
             ];
           If[dim == 1,
              body = "return decay_table[" <> ToString[rows[[1]]] <> "];\n";
              ,
              body = "switch (" <> idxName <> ") {\n";
              For[i = 0, i < dim, i++,
                  body = body <> "case " <> ToString[i] <> ": return decay_table[" <> ToString[rows[[i + 1]]] <> "]; break;\n";
                 ];
              body = body <> "}\n\n";
              errMsg = "std::ostringstream sstr;\n" <>
                       "sstr << \"invalid particle index \" << std::to_string(" <> idxName <> ") << '\\n';\n\n" <>
                       "throw OutOfBoundsError(sstr.str());\n";
              body = body <> errMsg;
             ];
           body
          ];

CreateDecayTableEntryGetterFunction[particle_, rows_List, scope_:"CLASSNAME"] :=
    Module[{dim, body},
           dim = TreeMasses`GetDimension[particle];
           body = CreateDecayTableEntryGetterFunctionBody[particle, rows];
           "Decays_list& " <> scope <> If[scope != "", "::", ""] <>
           CreateDecayTableEntryGetterName[particle] <> "(" <>
           If[dim > 1, "int gI1", ""] <> ")\n{\n" <>
           TextFormatting`IndentText[body] <> "}"
          ];

CreateDecayTableEntryConstGetterFunction[particle_, rows_List, scope_:"CLASSNAME"] :=
    Module[{dim, body},
           dim = TreeMasses`GetDimension[particle];
           body = CreateDecayTableEntryGetterFunctionBody[particle, rows];
           "const Decays_list& " <> scope <> If[scope != "", "::", ""] <>
           CreateDecayTableEntryGetterName[particle] <> "(" <>
           If[dim > 1, "int gI1", ""] <> ") const\n{\n" <>
           TextFormatting`IndentText[body] <> "}"
          ];

CreateDecayTableEntryGetterFunction[particle_, row_Integer, scope_:"CLASSNAME"] :=
    CreateDecayTableEntryGetterFunction[particle, {row}, scope];

CreateDecayTableEntryConstGetterFunction[particle_, row_Integer, scope_:"CLASSNAME"] :=
    CreateDecayTableEntryConstGetterFunction[particle, {row}, scope];

CreateDecayTableGetterPrototypes[decayParticles_List] :=
    Utils`StringJoinWithSeparator[(CreateDecayTableEntryGetterPrototype[#] <> "\n" <>
                                   CreateDecayTableEntryConstGetterPrototype[#])& /@ decayParticles, "\n"];

CreateDecayTableGetterFunctions[decayParticles_List, scope_:"CLASSNAME"] :=
    Module[{i, dims, offsets, rowAssignments, defs = ""},
           dims = TreeMasses`GetDimension /@ decayParticles;
           offsets = If[Length[dims] == 1, {0}, Join[{0}, Accumulate[dims[[2;;]]]]];
           rowAssignments = MapIndexed[{decayParticles[[First[#2]]], Table[offsets[[First[#2]]] + i, {i, 0, #1 - 1}]}&, dims];
           defs = (CreateDecayTableEntryGetterFunction[#[[1]], #[[2]], scope] <> "\n\n" <>
                   CreateDecayTableEntryConstGetterFunction[#[[1]], #[[2]], scope])& /@ rowAssignments;
           Utils`StringJoinWithSeparator[defs, "\n"]
          ];

CreateDecayTableInitialization[decayParticles_List] :=
    Module[{i, dims, dimsWithoutGoldstones, starts, pdgCodes, initializerList = ""},
           dims = TreeMasses`GetDimension /@ decayParticles;
           dimsWithoutGoldstones = TreeMasses`GetDimensionWithoutGoldstones /@ decayParticles;
           starts = TreeMasses`GetDimensionStartSkippingGoldstones /@ decayParticles;
           pdgCodes = Parameters`GetPDGCodesForParticle /@ decayParticles;
           For[i = 1, i <= Length[decayParticles], i++,
               If[dims[[i]] != Length[pdgCodes[[i]]],
                  Print["Error: number of PDG codes does not match size of ", decayParticles[[i]], " multiplet."];
                  Quit[1];
                 ];
               If[dimsWithoutGoldstones[[i]] > 0,
                  initializerList = initializerList <> If[initializerList == "", "", ", "] <>
                                    Utils`StringJoinWithSeparator[("Decays_list(" <> ToString[#] <> ")")& /@ pdgCodes[[i, starts[[i]] ;;]], ", "];
                 ];
              ];
           ": decay_table({" <> initializerList <> "})"
          ];

End[]

EndPackage[];
