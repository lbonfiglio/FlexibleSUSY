(* ::Package:: *)

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

BeginPackage["CXXDiagrams`", {"SARAH`", "TextFormatting`", "TreeMasses`", "Vertices`", "Parameters`","CConversion`", "SelfEnergies`", "Utils`"}];

(* This module generates c++ code intended to be used similarly to SARAH's fields and Vertex[] function *)

FeynmanDiagramsOfType::usage="";
VerticesForDiagram::usage="";

AtomHead::usage="";
LorentzConjugateOperation::usage="";
LorentzConjugate::usage="";
RemoveLorentzConjugation::usage="";
CreateFieldStructs::usage="";
CreateNamedFieldAliases::usage="";
CreateSelfConjugateFieldsDefinitions::usage="";
CreateFieldTypeLists::usage="";
CreateFieldTraitsDefinitions::usage="";
CreateMassFunctions::usage="";
CreatePhysicalMassFunctions::usage="";
CreateUnitCharge::usage="";
CreateStrongCoupling::usage="";
NumberOfFieldIndices::usage="";
FieldInfo::usage="";
includeLorentzIndices::usage="";

Begin["`Private`"];

(* TODO: Better name than LorentzConjugate *)
LorentzConjugateOperation[field_] := If[FermionQ[field] || GhostQ[field],
                                        "bar",
                                        "conj"]
LorentzConjugate[field_] := SARAH`AntiField[field]

RemoveLorentzConjugation[p_] := p
RemoveLorentzConjugation[SARAH`bar[p_]] := p
RemoveLorentzConjugation[Susyno`LieGroups`conj[p_]] := p

AtomHead[x_] := If[AtomQ[x], x, AtomHead[Head[x]]];

ParticleTypeAsString::argx = "Unknown type of particle `1`. Supported types are scalar, fermion, vector and ghost.";
ParticleTypeAsString[part_] := Module[
   {},
   If[TreeMasses`IsScalar[part],    Return["scalar"]];
   If[TreeMasses`IsVector[part],    Return["vector"]];
   If[TreeMasses`IsFermion[part],   Return["fermion"]];
   If[TreeMasses`IsGhost[part],   Return["ghost"]];

   Message[ParticleTypeAsString::argx, part]; Abort[];
];

(* todo: shouldn't there be an anti-triplet? *)
ParticleColorRepAsString::argx = "Unknown color representation `1` for particle `2`. Supported representations are singlet, (anti-)triplet and octet.";
ParticleColorRepAsString[part_] :=
   Module[{rep = TreeMasses`GetColorRepresentation[part]},
      Switch[rep,
         S, "singlet",
         T, "triplet",
         -T, "anti_triplet",
         O, "octet",
         _, Message[ParticleColorRepAsString::argx, rep, part]; Abort[];
      ]
   ];

CreateFieldIndexBoundsDefinition[field_] :=
    Module[{indexBounds},
           indexBounds = IndexBoundsForField[field];
           "using index_bounds = boost::mpl::pair<\n" <>
              "  boost::mpl::vector_c<int" <>
                   StringJoin[", " <> ToString[#] & /@
                     (indexBounds[[1]] - 1)] <> ">,\n" <>
              "  boost::mpl::vector_c<int" <>
                   StringJoin[", " <> ToString[#] & /@
                     indexBounds[[2]]] <> ">\n" <>
              ">;\n"
          ];

CreateFieldSMFlags[field_] :=
    Module[{dim = TreeMasses`GetDimension[field], flags = ""},
           flags = If[dim == 1,
                      CConversion`CreateCBoolValue @ TreeMasses`IsSMParticle[field],
                      StringJoin @ Riffle[CConversion`CreateCBoolValue /@
                                          (TreeMasses`IsSMParticle[#] & /@ Table[field[{k}], {k, dim}]),
                                          ", "]
                     ];
           "using sm_flags = boost::mpl::vector_c<bool, " <> flags <> ">;"
          ];

CreateFieldStruct[field_] :=
    Module[{fieldName = TreeMasses`CreateFieldClassName[field],
            particleType = "", particleColorRep = "", particleMassless = "",
            indexBounds = "", numGenerations = "", smFlags = "", numIndices = "",
            electricCharge = "", conjugateField = "", body = ""},
           particleType = "static constexpr auto particle_type = field_traits::ParticleType::" <> ParticleTypeAsString[field] <> ";";
           particleColorRep = "static constexpr auto color_rep = field_traits::ParticleColorRep::" <> ParticleColorRepAsString[field] <> ";";
           particleMassless = "static constexpr auto massless = " <> If[TreeMasses`IsMassless[field], "true", "false"] <> ";";
           indexBounds = CreateFieldIndexBoundsDefinition[field];
           numGenerations = "static constexpr int numberOfGenerations = " <> ToString @ TreeMasses`GetDimension[field] <> ";";
           smFlags = CreateFieldSMFlags[field];
           numIndices = "static constexpr int numberOfFieldIndices = " <> ToString @ NumberOfFieldIndices[field] <> ";";
           electricCharge = "static constexpr double electric_charge = " <>
                            CConversion`RValueToCFormString[TreeMasses`GetElectricCharge[field]] <> ";";
           conjugateField = "using lorentz_conjugate = " <> TreeMasses`CreateFieldClassName[LorentzConjugate[field]] <> ";";

           body = particleType <> "\n" <>
                  particleColorRep <> "\n" <>
                  particleMassless <> "\n" <>
                  indexBounds <> "\n" <>
                  numGenerations <> "\n" <>
                  smFlags <> "\n" <>
                  numIndices <> "\n" <>
                  electricCharge <> "\n" <>
                  conjugateField <> "\n";

           "struct " <> fieldName <> " {\n" <> TextFormatting`IndentText[body] <> "};\n"
          ];

CreateFieldStructs[fields_List] :=
    Utils`StringJoinWithSeparator[CreateFieldStruct /@ fields, "\n"];

CreateNamedFieldAliases[] :=
       "// Named fields\n" <>
       "using Photon = " <> TreeMasses`CreateFieldClassName[SARAH`Photon] <> ";\n" <>
       "using Electron = " <> TreeMasses`CreateFieldClassName[AtomHead @ TreeMasses`GetSMElectronLepton[]] <> ";\n";

IsLorentzSelfConjugate[field_] := field === LorentzConjugate[field];

CreateSelfConjugateFieldDefinition[field_, namespacePrefix_] := "";
CreateSelfConjugateFieldDefinition[field_?IsLorentzSelfConjugate, namespacePrefix_] :=
    Module[{fieldName},
           fieldName = TreeMasses`CreateFieldClassName[field, prefixNamespace -> namespacePrefix];
           "template<> struct " <> LorentzConjugateOperation[field] <> "<" <> fieldName <> ">" <>
           " { using type = " <> fieldName <> "; };"
          ];

CreateSelfConjugateFieldsDefinitions[fields_List, namespacePrefix_:""] :=
    Utils`StringJoinWithSeparator[CreateSelfConjugateFieldDefinition[#, namespacePrefix]& /@ Select[fields, IsLorentzSelfConjugate   ], "\n"] <> "\n";

CreateFieldTypeLists[fields_] :=
    Module[{scalars, fermions, vectors, ghosts},
           scalars = Select[fields, TreeMasses`IsScalar];
           fermions = Select[fields, TreeMasses`IsFermion];
           vectors = Select[fields, TreeMasses`IsVector];
           ghosts = Select[fields, TreeMasses`IsGhost];

           "using scalars = boost::mpl::vector<" <>
           Utils`StringJoinWithSeparator[TreeMasses`CreateFieldClassName /@ scalars, ", "] <> ">;\n" <>
           "using fermions = boost::mpl::vector<" <>
           Utils`StringJoinWithSeparator[TreeMasses`CreateFieldClassName /@ fermions, ", "] <> ">;\n" <>
           "using vectors = boost::mpl::vector<" <>
           Utils`StringJoinWithSeparator[TreeMasses`CreateFieldClassName /@ vectors, ", "] <> ">;\n" <>
           "using ghosts = boost::mpl::vector<" <>
           Utils`StringJoinWithSeparator[TreeMasses`CreateFieldClassName /@ ghosts, ", "] <> ">;\n"
          ];

CreateFieldTypeTraitDefinition[field_?TreeMasses`IsScalar, namespacePrefix_] :=
    "template<>\nstruct is_scalar<" <> TreeMasses`CreateFieldClassName[field, prefixNamespace -> namespacePrefix] <> " > : public std::true_type {};";
CreateFieldTypeTraitDefinition[field_?TreeMasses`IsFermion, namespacePrefix_] :=
    "template<>\nstruct is_fermion<" <> TreeMasses`CreateFieldClassName[field, prefixNamespace -> namespacePrefix] <> " > : public std::true_type {};";
CreateFieldTypeTraitDefinition[field_?TreeMasses`IsVector, namespacePrefix_] :=
    "template<>\nstruct is_vector<" <> TreeMasses`CreateFieldClassName[field, prefixNamespace -> namespacePrefix] <> " > : public std::true_type {};";
CreateFieldTypeTraitDefinition[field_?TreeMasses`IsGhost, namespacePrefix_] :=
    "template<>\nstruct is_ghost<" <> TreeMasses`CreateFieldClassName[field, prefixNamespace -> namespacePrefix] <> " > : public std::true_type {};";

CreateFieldTraitsDefinition[field_, namespacePrefix_] :=
    Module[{fieldTypeTraitDefinition = ""},
           fieldTypeTraitDefinition = CreateFieldTypeTraitDefinition[field, namespacePrefix];
           fieldTypeTraitDefinition <> "\n"
          ];

CreateFieldTraitsDefinitions[fields_, namespacePrefix_:""] :=
    StringJoin[CreateFieldTraitsDefinition[#, namespacePrefix]& /@ fields];

(* adjacencyMatrix must be undirected (i.e symmetric) *)
FeynmanDiagramsOfType[adjacencyMatrix_List,externalFields_List] :=
  Module[{externalVertices = externalFields[[All,1]],
          internalVertices,externalRules,
          internalFieldCouplings,
          unspecifiedEdgesLess,unspecifiedEdgesEqual,
          insertFieldRulesLess,insertFieldRulesGreater,insertFieldRulesEqual,
          fieldsToInsert,
          unresolvedFieldCouplings,resolvedFields,resolvedFieldCouplings,
          diagrams = {}},
   internalVertices = Complement[Table[k,{k,Length[adjacencyMatrix]}],externalVertices];
   externalRules = Flatten @ ({{_,#,_} :> SARAH`AntiField[# /. externalFields],
                               {#,_,_} :> SARAH`AntiField[# /. externalFields]} & /@ externalVertices);

   internalFieldCouplings = (Flatten[(Flatten @ Position[adjacencyMatrix[[#]],Except[0],{1},Heads -> False]
                                /. {i_Integer :> Table[{#,i,k},{k,adjacencyMatrix[[#,i]]}]}),1] &
                             /@ internalVertices) /. externalRules;

   unspecifiedEdgesLess = Cases[internalFieldCouplings,{i_,j_,_} /; i < j,{2}];
   unspecifiedEdgesEqual = Cases[internalFieldCouplings,{i_,i_,_},{2}];

   If[unspecifiedEdgesLess === {} && unspecifiedEdgesEqual === {},
      diagrams = Table[k, {k, 1, Length[adjacencyMatrix]}] /. externalFields /.
          {Thread[Rule[internalVertices, internalFieldCouplings]]};
      ,
      insertFieldRulesLess = MapIndexed[#1 -> SARAH`FieldToInsert[#2[[1]]] &,unspecifiedEdgesLess];
      insertFieldRulesGreater = (insertFieldRulesLess /. {Rule[{i_,j_,k_},field_] :> Rule[{j,i,k},SARAH`AntiField[field]]});
      insertFieldRulesEqual = MapIndexed[#1 -> {SARAH`FieldToInsert[#2[[1]]+Length[insertFieldRulesLess]],
                                               SARAH`AntiField[SARAH`FieldToInsert[#2[[1]]+Length[insertFieldRulesLess]]]} &,
                                         unspecifiedEdgesEqual];
      fieldsToInsert = Table[SARAH`FieldToInsert[k],
                {k,Length[insertFieldRulesLess] + Length[insertFieldRulesEqual]}];

      unresolvedFieldCouplings = internalFieldCouplings
        /. insertFieldRulesLess /. insertFieldRulesGreater /. insertFieldRulesEqual;
      resolvedFields = SARAH`InsFields[{C @@@ unresolvedFieldCouplings,
                                        fieldsToInsert}][[All,2]];
      If[resolvedFields =!= {},
         resolvedFieldCouplings = unresolvedFieldCouplings /.
           ((Rule @@@ Transpose[{fieldsToInsert,#}]) & /@ resolvedFields);

         diagrams = Table[k,{k,Length[adjacencyMatrix]}] /. externalFields /.
           ((Rule @@@ Transpose[{internalVertices,#}]) & /@ resolvedFieldCouplings);
        ];
     ];

   DeleteDuplicates[diagrams,
     AllTrue[Cases[Transpose[{#1,#2}],{{___},{___}}], (* External lines *)
             (Sort[#[[1]]] === Sort[#[[2]]]&)] &]
  ];

VerticesForDiagram[diagram_] := Select[diagram,Length[#] > 1 &]

CreateMassFunctions[fieldsNamespace_:""] :=
  Module[{massiveFields,
          ghostMappings = SelfEnergies`ReplaceGhosts[FlexibleSUSY`FSEigenstates]},
    massiveFields = TreeMasses`GetParticles[];

    StringJoin @ Riffle[
      Module[{fieldInfo = FieldInfo[#], numberOfIndices},
             numberOfIndices = Length @ fieldInfo[[5]];

             "template<> inline\n" <>
             "double " <> FlexibleSUSY`FSModelName <> "_context_base::mass_impl<" <>
               TreeMasses`CreateFieldClassName[#, prefixNamespace -> fieldsNamespace] <>
             ">(const std::array<int, " <> ToString @ numberOfIndices <>
             ">&" <> If[TreeMasses`GetDimension[#] === 1, "", " indices"] <> ") const\n" <>
             "{ return model.get_M" <> TreeMasses`CreateFieldClassName[# /. ghostMappings] <>
             If[TreeMasses`GetDimension[#] === 1, "()", "(indices[0])"] <> "; }"
            ] & /@ massiveFields, "\n\n"]
        ]

CreatePhysicalMassFunctions[fieldsNamespace_:""] :=
  Module[{massiveFields,
          ghostMappings = SelfEnergies`ReplaceGhosts[FlexibleSUSY`FSEigenstates]},
    massiveFields = TreeMasses`GetParticles[];

    StringJoin @ Riffle[
      Module[{fieldInfo = FieldInfo[#], numberOfIndices},
             numberOfIndices = Length @ fieldInfo[[5]];

             "template<> inline\n" <>
             "double " <> FlexibleSUSY`FSModelName <> "_context_base::physical_mass_impl<" <>
               TreeMasses`CreateFieldClassName[#, prefixNamespace -> fieldsNamespace] <>
             ">(const std::array<int, " <> ToString @ numberOfIndices <>
             ">&" <> If[TreeMasses`GetDimension[#] === 1, "", " indices"] <> ") const\n" <>
             "{ return model.get_physical().M" <> TreeMasses`CreateFieldClassName[# /. ghostMappings] <>
             If[TreeMasses`GetDimension[#] === 1, "", "[indices[0]]"] <> "; }"
            ] & /@ massiveFields, "\n\n"]
        ]

CreateUnitCharge[] :=
  Module[{electron,photon,vertex,vertexBody,
          numberOfElectronIndices,numberOfPhotonIndices},
         electron = AtomHead @ TreeMasses`GetSMElectronLepton[];
         photon = SARAH`Photon;
         vertex = {SARAH`bar[electron], electron, photon};
         vertexBody = Vertices`VertexFunctionBodyForFields[vertex];
         numberOfElectronIndices = NumberOfFieldIndices[electron];
         numberOfPhotonIndices = NumberOfFieldIndices[photon];

         "static FFSVertex unit_charge(const " <> FlexibleSUSY`FSModelName <> "_context_base& context)\n" <>
         "{\n" <>
         TextFormatting`IndentText["using vertex_type = FFSVertex;"] <> "\n\n" <>
         TextFormatting`IndentText @
           ("std::array<int, " <> ToString @ numberOfElectronIndices <> "> electron_indices = {" <>
              If[TreeMasses`GetDimension[electron] =!= 1,
                 " " <> ToString @ (FieldInfo[electron][[2]]-1) <> (* Electron has the lowest index *)
                 If[numberOfElectronIndices =!= 1,
                    StringJoin @ Table[", 0", {numberOfElectronIndices-1}],
                    ""] <> " ",
                 If[numberOfElectronIndices =!= 0,
                    StringJoin @ Riffle[Table[" 0", {numberOfElectronIndices}], ","] <> " ",
                    ""]
                ] <>
            "};\n") <>
         TextFormatting`IndentText @
           ("std::array<int, " <> ToString @ numberOfPhotonIndices <> "> photon_indices = {" <>
               If[TreeMasses`GetDimension[photon] =!= 1,
                 " " <> ToString @ (FieldInfo[photon][[2]]-1) <>
                 If[numberOfPhotonIndices =!= 1,
                    StringJoin @ Table[", 0", {numberOfPhotonIndices-1}],
                    ""] <> " ",
                 If[numberOfPhotonIndices =!= 0,
                    StringJoin @ Riffle[Table[" 0", {numberOfPhotonIndices}], ","] <> " ",
                    ""]
                ] <>
            "};\n") <>
         TextFormatting`IndentText @
           ("std::array<int, " <> ToString[
              Total[NumberOfFieldIndices /@ {photon,electron,electron}]] <>
            "> indices = " <>
              "concatenate(photon_indices, electron_indices, electron_indices);\n\n") <>
         TextFormatting`IndentText @ vertexBody <> "\n" <>
         "}"
  ]

CreateStrongCoupling[] :=
  Module[{downquark,gluon,vertex,
          vertexBody, numberOfdownquarkIndices,numberOfgluonIndices},
         downquark = AtomHead @ TreeMasses`GetSMDownQuark[];
         gluon = SARAH`Gluon;
         vertex = {gluon, downquark, SARAH`bar[downquark]};
         vertexBody = Vertices`VertexFunctionBodyForFields[vertex];
         numberOfdownquarkIndices = NumberOfFieldIndices[downquark];
         numberOfgluonIndices = NumberOfFieldIndices[gluon];

         "static FFSVertex strong_coupling(const " <> FlexibleSUSY`FSModelName <> "_context_base& context)\n" <>
         "{\n" <>
         TextFormatting`IndentText["using vertex_type = FFSVertex;"] <> "\n\n" <>
         TextFormatting`IndentText @
           ("std::array<int, " <> ToString @ numberOfdownquarkIndices <> "> downquark_indices = {" <>
              If[TreeMasses`GetDimension[downquark] =!= 1,
                 " " <> ToString @ (FieldInfo[downquark][[2]]-1) <> (* downquark has the lowest index *)
                 If[numberOfdownquarkIndices =!= 1,
                    StringJoin @ Table[", 0", {numberOfdownquarkIndices-1}],
                    ""] <> " ",
                 If[numberOfdownquarkIndices =!= 0,
                    StringJoin @ Riffle[Table[" 0", {numberOfdownquarkIndices}], ","] <> " ",
                    ""]
                ] <>
            "};\n") <>
         TextFormatting`IndentText @
           ("std::array<int, " <> ToString @ numberOfgluonIndices <> "> gluon_indices = {" <>
               If[TreeMasses`GetDimension[gluon] =!= 1,
                 " " <> ToString @ (FieldInfo[gluon][[2]]-1) <>
                 If[numberOfgluonIndices =!= 1,
                    StringJoin @ Table[", 0", {numberOfgluonIndices-1}],
                    ""] <> " ",
                 If[numberOfgluonIndices =!= 0,
                    StringJoin @ Riffle[Table[" 0", {numberOfgluonIndices}], ","] <> " ",
                    ""]
                ] <>
            "};\n") <>
         TextFormatting`IndentText @
           ("std::array<int, " <> ToString[
              Total[NumberOfFieldIndices /@ {gluon,downquark,downquark}]] <>
            "> indices = " <>
              "concatenate(gluon_indices, downquark_indices, downquark_indices);\n\n") <>

         TextFormatting`IndentText @ vertexBody <> "\n" <>
         "}"
  ]

NumberOfFieldIndices[field_] := Length @ FieldInfo[field][[5]]

IndexBoundsForField[field_] :=
  Module[{fieldInfo = FieldInfo[field]},
    If[NumberOfFieldIndices[field] === 0,
       Return[{{},{}}]];
    If[Length @ Cases[fieldInfo[[5]],{SARAH`generation,_}] === 0,
       Transpose[{1,#[[2]]} & /@ fieldInfo[[5]]],
       Transpose[Prepend[
         {1,#[[2]]} & /@ DeleteCases[fieldInfo[[5]],{SARAH`generation,_}],
         {fieldInfo[[2]],fieldInfo[[3]]}]]]]

End[];
EndPackage[];
