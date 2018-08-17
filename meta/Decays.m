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

CreateCompleteParticleList::usage="";
GetDecaysForParticle::usage="";

CallDecaysCalculationFunctions::usage="creates calls to functions calculating
decays of the given particles.";
CreateDecaysCalculationPrototypes::usage="creates prototypes for convenience
functions calculating all decays for each particle.";
CreateDecaysCalculationFunctions::usage="creates definitions for convenience
functions calculating all decays for each particle.";

FSParticleDecay::usage="head used for storing details of an particle decay,
in the format
   FSParticleDecay[particle, {final state particle}, {diagram 1, diagram 2, ...}]
";

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

IsElectricChargeConservingDecay[initialParticle_, finalState_List] :=
    Module[{chargeSum},
           chargeSum = Simplify[Plus @@ (Join[{-TreeMasses`GetElectricCharge[initialParticle]}, TreeMasses`GetElectricCharge /@ finalState])];
           PossibleZeroQ[chargeSum]
          ];

GetContributingDiagramsForDecayGraph[initialField_, finalFields_List, graph_] :=
    Module[{externalFields, diagrams},
           externalFields = Join[{1 -> initialField}, MapIndexed[(First[#2] + 1 -> #1)&, finalFields]];
           CXXDiagrams`FeynmanDiagramsOfType[graph, externalFields]
          ];

GetContributingGraphsForDecay[initialParticle_, finalParticles_List] :=
    Module[{nFinalParticles = Length[finalParticles], topologies},
           topologies = GetDecayTopologies[nFinalParticles];
           Flatten[GetContributingDiagramsForDecayGraph[initialParticle, finalParticles, #]& /@ topologies, 1]
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
           isPossibleDecay[finalState_] := (IsElectricChargeConservingDecay[particle, finalState]);
           concreteFinalStates = Join @@ (GetParticleCombinationsOfType[#, allowedFinalStateParticles, isPossibleDecay]& /@ genericFinalStates);
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

(*
CreatePartialWidthCalculationName[decayParticle_, productParticles_List] :=
    Module[{inParticleNoIndices, outParticlesNoIndices},
           inParticleNoIndices = Vertices`StripFieldIndices[decayParticle];
           outParticlesNoIndices = Vertices`StripFieldIndices /@ productParticles;

           "get_partial_width_" <> CConversion`ToValidCSymbolString[inParticleNoIndices] <>
           "_to_" <> StringJoin[CConversion`ToValidCSymbolString /@ outParticlesNoIndices]
          ];

CreatePartialWidthCalculationProtoype[decayParticle_, productParticles_List] :=
    "double " <> 

CreatePartialWidthCalculationFunction[decayParticle_, productParticles_List] :=
*)

CreatePartialWidthCalculationName[initialState_, finalState_List] :=
    "get_partial_width<" <> CConversion`ToValidCSymbolString[initialState] <> "," <>
    Utils`StringJoinWithSeparator[CXXDiagrams`CXXNameOfField /@ finalState, ","] <> " >";

CallPartialWidthCalculation[decay_FSParticleDecay] :=
    Module[{initialState = GetInitialState[decay],
            finalState = GetFinalState[decay], result = ""},
           result = CreatePartialWidthCalculationName[initialState, finalState] <> "();\n";
           result
          ];

CreateDecaysCalculationFunctionName[particle_, scope_:""] :=
    scope <> If[scope =!= "", "::", ""] <> "calculate_" <> CConversion`ToValidCSymbolString[particle] <> "_decays";

CreateDecaysCalculationPrototype[particle_] :=
    "void " <> CreateDecaysCalculationFunctionName[particle] <> "(const " <>
    FlexibleSUSY`FSModelName <> "_mass_eigenstates&);";

CreateDecaysCalculationPrototypes[decayLists_List] :=
    Utils`StringJoinWithSeparator[CreateDecaysCalculationPrototype[First[#]]& /@ decayLists, "\n"];

LoopOverIndex[loopBody_String, index_, start_, stop_, type_:CConversion`ScalarType[CConversion`integerScalarCType]] :=
    Module[{idxStr, startStr, stopStr},
           idxStr = CConversion`ToValidCSymbolString[index];
           startStr = CConversion`ToValidCSymbolString[start];
           stopStr = CConversion`ToValidCSymbolString[stop];
           "for (" <> CConversion`CreateCType[type] <> " " <> idxStr <>
           " = " <> startStr <> "; " <> idxStr <> " < " <> stopStr <> "; ++" <>
           idxStr <> ") {\n" <> TextFormatting`IndentText[loopBody] <> "\n}\n"
          ];

(* generates a loop over the given indices, in the form
   {{idx1, start, stop}, {idx2, start, stop}, ...} with
   the first list entry being the innermost loop *)
LoopOverIndexCollection[loopBody_String, indices_List] :=
    Fold[LoopOverIndex[#1, Sequence @@ #2]&, loopBody, indices];

CreateDecaysCalculationFunction[decaysList_] :=
    Module[{particle = First[decaysList], decayChannels = Last[decaysList],
            dim, dimStart, loopIndices, result = "", body = ""},
           dimStart = TreeMasses`GetDimensionStartSkippingGoldstones[particle] - 1;
           dim = TreeMasses`GetDimension[particle];
           loopIndices = If[dim > 1 || TreeMasses`GetDimensionWithoutGoldstones[particle] > 1,
                            {{"gI1", dimStart, dim}},
                            {}
                           ];
           body = "auto model = model_;\n\n" <> StringJoin[CallPartialWidthCalculation /@ decayChannels];
           body = LoopOverIndexCollection[body, loopIndices];
           "void " <> CreateDecaysCalculationFunctionName[particle, "CLASSNAME"] <>
           "(const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model_)\n{\n"
           <> TextFormatting`IndentText[body] <> "}\n"
          ];

CreateDecaysCalculationFunctions[particleDecays_List] :=
    Utils`StringJoinWithSeparator[CreateDecaysCalculationFunction /@ particleDecays, "\n"];

CallDecaysFunction[particle_, arg_:"model", obj_:""] :=
    obj <> CreateDecaysCalculationFunctionName[particle] <> "(" <> arg <> ");\n"

CallThreadedDecaysFunction[particle_, ptr_:"this", pool_:"tp", arg_:"model"] :=
    pool <> ".run_task([" <> ptr <> ", &" <> arg <> "] () { " <>
    If[ptr === "this", "", ptr <> "->"] <>
    CreateDecaysCalculationFunctionName[particle] <> "(" <> arg <> "); });\n";

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

End[]

EndPackage[];
