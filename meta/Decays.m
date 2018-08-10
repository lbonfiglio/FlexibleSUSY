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

CallDecaysCalculationFunctions::usage="creates calls to functions calculating
decays of the given particles.";
CreateDecaysCalculationPrototypes::usage="creates prototypes for convenience
functions calculating all decays for each particle.";
CreateDecaysCalculationFunctions::usage="creates definitions for convenience
functions calculating all decays for each particle.";

Begin["`Private`"];

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

CreateDecaysCalculationFunctionName[particle_, scope_:""] :=
    scope <> If[scope =!= "", "::", ""] <> "calculate_" <> CConversion`ToValidCSymbolString[particle] <> "_decays";

CreateDecaysCalculationPrototype[particle_] :=
    "void " <> CreateDecaysCalculationFunctionName[particle] <> "(const " <>
    FlexibleSUSY`FSModelName <> "_mass_eigenstates&);";

CreateDecaysCalculationPrototypes[particles_List] :=
    Utils`StringJoinWithSeparator[CreateDecaysCalculationPrototype /@ particles, "\n"];

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

CreateDecaysCalculationFunction[particle_] :=
    Module[{generationIndex, dim, dimStart, loopIndices, result = "", body = ""},
           dimStart = TreeMasses`GetDimensionStartSkippingGoldstones[particle];
           dim = TreeMasses`GetDimension[particle];
           generationIndex = Select[TreeMasses`GetParticleIndices[particle], CXXDiagrams`IsGenerationIndex];
           If[Length[generationIndex] > 1,
              generationIndex = First[generationIndex];
             ];
           loopIndices = If[generationIndex === {} || dim == 0, {}, {generationIndex, dimStart, dim}];
           Print["dimStart = ", dimStart];
           Print["dim = ", dim];
           Print["indices = ", TreeMasses`GetParticleIndices[particle]];
           Print["generationIndex = ", generationIndex];
           Print["loopIndices = ", loopIndices];
           body = "auto model = model_;\n";
           body = LoopOverIndexCollection[body, loopIndices];
           "void " <> CreateDecaysCalculationFunctionName[particle, "CLASSNAME"] <>
           "(const " <> FlexibleSUSY`FSModelName <> "_mass_eigenstates& model_)\n{\n"
           <> TextFormatting`IndentText[body] <> "}\n"
          ];

CreateDecaysCalculationFunctions[particles_List] :=
    Utils`StringJoinWithSeparator[CreateDecaysCalculationFunction /@ particles, "\n"];

CallDecaysFunction[particle_, arg_:"model", obj_:""] :=
    obj <> CreateDecaysCalculationFunctionName[particle] <> "(" <> arg <> ");\n"

CallThreadedDecaysFunction[particle_, ptr_:"this", pool_:"tp", arg_:"&model"] :=
    pool <> ".run_task([" <> ptr <> ", " <> arg <> "] () { " <>
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
