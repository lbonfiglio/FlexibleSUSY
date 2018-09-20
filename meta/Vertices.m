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

BeginPackage["Vertices`", {
    "SARAH`",
    "SelfEnergies`",
    "Parameters`",
    "TreeMasses`",
    "CConversion`",
    "LatticeUtils`"}]

FSVertexTypes = { SSSVertex, SSSSVertex, FFSVertex, FFVVertex, SSVVertex, SSVVVertex, SVVVertex, VVVVertex, VVVVVertex, GGSVertex, GGVVertex };
VertexTypes::usage="";
VertexTypeForFields::usage="Returns the vertex type for a vertex with a given list of fields.";

VertexRules::usage;
ToCpPattern::usage="ToCpPattern[cp] converts field indices inside cp to patterns, e.g. ToCpPattern[Cp[bar[UFd[{gO1}]], Sd[{gI1}], Glu[{1}]][PL]] === Cp[bar[UFd[{gO1_}]], Sd[{gI1_}], Glu[{1}]][PL].";
ToCp::usage="ToCp[cpPattern] converts field index patterns inside cpPattern to symbols, e.g. ToCp@Cp[bar[UFd[{gO1_}]], Sd[{gI1_}], Glu[{1}]][PL] === Cp[bar[UFd[{gO1}]], Sd[{gI1}], Glu[{1}]][PL].";
FieldIndexList::usage;
SortCp::usage="SortCp[cp] sorts fields in cp into SARAH internal order.";
SortCps::usage="SortCps[nPointFunctions] sorts all SARAH`Cp[] and SARAH`Cp[][] in nPointFunctions.";
EnforceCpColorStructures::usage;
EnforceCpColorStructures::cpext="Fixing positions of external field `1` within `2`.  This might happen with SARAH version 4.1.0 or earlier.  Please report to us if you see this message with a newer version of SARAH.";

GetLorentzStructure::usage;
GetParticleList::usage;
IsUnrotated::usage;
ToRotatedField::usage;
ReplaceUnrotatedFields::usage;
StripGroupStructure::usage="Removes group generators and Kronecker deltas.";
StripFieldIndices::usage;

CreateVertexData::usage="";
CreateVertices::usage="";
VertexFunctionBodyForFields::usage="";

Begin["`Private`"]

VertexTypes[] := FSVertexTypes;

IsSSSVertex[fields_List] :=
    Module[{},
           If[Count[fields, _?TreeMasses`IsFermion] != 0 ||
              Count[fields, _?TreeMasses`IsVector]  != 0 ||
              Count[fields, _?TreeMasses`IsGhost]   != 0,
              Return[False];
             ];
           Count[fields, _?TreeMasses`IsScalar] == 3
          ];

IsSSSSVertex[fields_List] :=
    Module[{},
           If[Count[fields, _?TreeMasses`IsFermion] != 0 ||
              Count[fields, _?TreeMasses`IsVector]  != 0 ||
              Count[fields, _?TreeMasses`IsGhost]   != 0,
              Return[False];
             ];
           Count[fields, _?TreeMasses`IsScalar] == 4
          ];

IsFFSVertex[fields_List] :=
    Module[{fermionCount},
           If[Count[fields, _?TreeMasses`IsGhost]  != 0 ||
              Count[fields, _?TreeMasses`IsVector] != 0,
              Return[False];
             ];
           fermionCount = Count[fields, _?TreeMasses`IsFermion];
           If[fermionCount != 2,
              Return[False];
             ];
           Count[fields, _?TreeMasses`IsScalar] == 1
          ];

IsFFVVertex[fields_List] :=
    Module[{fermionCount},
           If[Count[fields, _?TreeMasses`IsGhost]  != 0 ||
              Count[fields, _?TreeMasses`IsScalar] != 0,
              Return[False];
             ];
           fermionCount = Count[fields, _?TreeMasses`IsFermion];
           If[fermionCount != 2,
              Return[False];
             ];
           Count[fields, _?TreeMasses`IsVector] == 1
          ];

IsSSVVertex[fields_List] :=
    Module[{scalarCount, vectorCount},
           If[Count[fields, _?TreeMasses`IsFermion] != 0 || Count[fields, _?TreeMasses`IsGhost] != 0,
              Return[False];
             ];
           scalarCount = Count[fields, _?TreeMasses`IsScalar];
           vectorCount = Count[fields, _?TreeMasses`IsVector];
           scalarCount == 2 && vectorCount == 1
          ];

IsSVVVertex[fields_List] :=
    Module[{scalarCount, vectorCount},
           If[Count[fields, _?TreeMasses`IsFermion] != 0 || Count[fields, _?TreeMasses`IsGhost] != 0,
              Return[False];
             ];
           scalarCount = Count[fields, _?TreeMasses`IsScalar];
           vectorCount = Count[fields, _?TreeMasses`IsVector];
           vectorCount == 2 && scalarCount == 1
          ];

IsSSVVVertex[fields_List] :=
    Module[{scalarCount, vectorCount},
           If[Count[fields, _?TreeMasses`IsFermion] != 0 || Count[fields, _?TreeMasses`IsGhost] != 0,
              Return[False];
             ];
           scalarCount = Count[fields, _?TreeMasses`IsScalar];
           vectorCount = Count[fields, _?TreeMasses`IsVector];
           vectorCount == 2 && scalarCount == 2
          ];

IsVVVVertex[fields_List] :=
    Module[{},
           If[Count[fields, _?TreeMasses`IsFermion] != 0 ||
              Count[fields, _?TreeMasses`IsScalar]  != 0 ||
              Count[fields, _?TreeMasses`IsGhost]   != 0,
              Return[False];
             ];
           Count[fields, _?TreeMasses`IsVector] == 3
          ];

IsVVVVVertex[fields_List] :=
    Module[{},
           If[Count[fields, _?TreeMasses`IsFermion] != 0 ||
              Count[fields, _?TreeMasses`IsScalar]  != 0 ||
              Count[fields, _?TreeMasses`IsGhost]   != 0,
              Return[False];
             ];
           Count[fields, _?TreeMasses`IsVector] == 4
          ];

IsGGSVertex[fields_List] :=
    Module[{ghostCount},
           If[Count[fields, _?TreeMasses`IsFermion] != 0 ||
              Count[fields, _?TreeMasses`IsVector]  != 0,
              Return[False];
             ];
           ghostCount = Count[fields, _?TreeMasses`IsGhost];
           If[ghostCount != 2,
              Return[False];
             ];
           Count[fields, _?TreeMasses`IsScalar] == 1
          ];

IsGGVVertex[fields_List] :=
    Module[{ghostCount},
           If[Count[fields, _?TreeMasses`IsFermion] != 0 ||
              Count[fields, _?TreeMasses`IsScalar]  != 0,
              Return[False];
             ];
           ghostCount = Count[fields, _?TreeMasses`IsGhost];
           If[ghostCount != 2,
              Return[False];
             ];
           Count[fields, _?TreeMasses`IsVector] == 1
          ];

VertexTypeForFields[fields_List /; IsSSSVertex[fields]] := SSSVertex;
VertexTypeForFields[fields_List /; IsSSSSVertex[fields]] := SSSSVertex;
VertexTypeForFields[fields_List /; IsFFSVertex[fields]] := FFSVertex;
VertexTypeForFields[fields_List /; IsFFVVertex[fields]] := FFVVertex;
VertexTypeForFields[fields_List /; IsSSVVertex[fields]] := SSVVertex;
VertexTypeForFields[fields_List /; IsSSVVVertex[fields]] := SSVVVertex;
VertexTypeForFields[fields_List /; IsSVVVertex[fields]] := SVVVertex;
VertexTypeForFields[fields_List /; IsVVVVertex[fields]] := VVVVertex;
VertexTypeForFields[fields_List /; IsVVVVVertex[fields]] := VVVVVertex;
VertexTypeForFields[fields_List /; IsGGSVertex[fields]] := GGSVertex;
VertexTypeForFields[fields_List /; IsGGVVertex[fields]] := GGVVertex;
VertexTypeForFields[fields_List] := "UnknownVertexType" <> ToString[fields];

(* There is a sign ambiguity when SARAH`Vertex[] factors an SSV-type
   vertex into a coefficient and a Lorentz part even though the
   product of the two is always the same.  SARAH seems to have an
   internal convention such that the Lorentz part always takes the
   form

      Mom[S1[{gt1,___}], _] - Mom[S2[{gt2,___}], _]

   if it is the result of SARAH`Vertex[{S1, S2, V}] where S1 and S2
   are scalars without (generation) indices.  This alone does not yet
   determine how the scalars in Cp[] are mapped to S1 and S2 above:
   which of conj[Sd] and Sd becomes S1 in Cp[VG, conj[Sd[{gI1}]],
   Sd[{gI2}]] for instance.  This seems to be fixed by another
   internal rule: all SARAH/SPheno subroutines to calculate couplings
   are generated with the fields sorted by SortCoup[] in
   SARAH/Package/deriveModel.m.  The same ordering is performed by
   Vertices`SortCp[] and Vertices`SortCps[] in FlexibleSUSY, which
   normalize e.g. the above coupling to Cp[Sd[{gI2}], conj[Sd[{gI1}]],
   VG].  VertexRules[] then passes the arguments of Cp[] to
   SARAH`Vertex[] with the indices gI2 and gI1 omitted, so that the
   mapping is determined to be {Sd -> S1, conj[Sd] -> S2}. *)

VertexRules[nPointFunctions_, massMatrices_] := Block[{
        UnitaryMatrixQ,
        nCpPatterns,
        cpPatterns = DeleteRedundantCpPatterns[
            ToCpPattern /@ RenumberCpIndices /@ Cases[
                nPointFunctions, _SARAH`Cp|_SARAH`Cp[_], {0, Infinity}]]
    },
    (UnitaryMatrixQ[#] = True)& /@
        Flatten@DeleteCases[massMatrices[[All,3]], Null];
    UnitaryMatrixQ[_] := False;
    nCpPatterns = Length[cpPatterns];
    MapIndexed[
        DoneLn[#1 -> VertexExp[#1, nPointFunctions, massMatrices],
               "[",First[#2],"/",nCpPatterns,"] calculating ", #1, "... "]&,
        cpPatterns]
];

SortCps[nPointFunctions_List] := Module[{
        exprs = nPointFunctions[[All,2]]
    },
    Fold[
        Module[{sortedCp = SortCp[#2]},
            If[sortedCp =!= #2, #1 /. #2 -> sortedCp, #1]] &,
        nPointFunctions,
        Union @ Select[Cases[exprs, _SARAH`Cp|_SARAH`Cp[_], Infinity],
                       UnresolvedColorFactorFreeQ[#, exprs] &]]
];

SortCp[SARAH`Cp[fields__]] :=
    SARAH`Cp @@ SortFieldsInCp @ StripExtraFieldIndices[{fields}];

SortCp[SARAH`Cp[fields__][lor_]] := SortCp[SARAH`Cp[fields]][lor];

(* see OrderVVVV[] in SARAH/Package/SPheno/SPhenoFunc.m *)
SortCp[cp : SARAH`Cp[vectors__][lor_Integer]] /; CpType[cp] === VVVV :=
Module[{
        vs = StripExtraFieldIndices[{vectors}],
        svs, lors,
        sortedVectors,
        ssvs, sortedLors,
        map
    },
    sortedVectors = SortFieldsInCp[vs];
    svs  = StripFieldIndices[vs];
    ssvs = StripFieldIndices[sortedVectors];
    lors = {
        SARAH`g[ svs[[1]],  svs[[2]]] SARAH`g[ svs[[3]],  svs[[4]]],
        SARAH`g[ svs[[1]],  svs[[3]]] SARAH`g[ svs[[2]],  svs[[4]]],
        SARAH`g[ svs[[1]],  svs[[4]]] SARAH`g[ svs[[2]],  svs[[3]]]
    };
    sortedLors = {
        SARAH`g[ssvs[[1]], ssvs[[2]]] SARAH`g[ssvs[[3]], ssvs[[4]]],
        SARAH`g[ssvs[[1]], ssvs[[3]]] SARAH`g[ssvs[[2]], ssvs[[4]]],
        SARAH`g[ssvs[[1]], ssvs[[4]]] SARAH`g[ssvs[[2]], ssvs[[3]]]
    };
    map = Ordering[sortedLors, 3, OrderedQ[First@Position[lors, #]& /@ {##}]&];
    (SARAH`Cp @@ sortedVectors)[map[[lor]]]
];

(* see WriteFermionProp[] in SARAH/Package/SPheno/SPhenoLoopMasses *)
SortCp[cp : SARAH`Cp[fields__][lor:PL|PR]] /; CpType[cp] === FFV := Module[{
        fs = StripExtraFieldIndices[{fields}],
        sorted,
        fermions, sortedFermions
    },
    sorted = SortFieldsInCp[fs];
    fermions       = Select[fs     , GetFieldType@ToRotatedField[#] === F &];
    sortedFermions = Select[sorted , GetFieldType@ToRotatedField[#] === F &];
    If[First[fermions] === First[sortedFermions],
          (SARAH`Cp @@ sorted)[lor],
        - (SARAH`Cp @@ sorted)[First @ Complement[{PL, PR}, {lor}]]]
];

SortFieldsInCp[fields_List] :=
    SortBy[fields, (GetTypeSort[#][#]& @ ToRotatedField[#]) &];

(* Same as SARAH`getTypeSort but works when
   SARAH`CurrentStates =!= FSEigenstates *)
GetTypeSort[Susyno`LieGroups`conj[x_]] :=
    Switch[GetFieldType[x],
        S, Szc,
        V, Vc ,
        A, Ab ];

GetTypeSort[SARAH`bar[x_]] :=
    Switch[GetFieldType[x],
        F, Fb,
        G, Gb];

GetTypeSort[x_ /; SARAH`bar[x] === x] /;
    GetFieldType[x] === F :=
           Fm;

GetTypeSort[x_] :=
    Switch[GetFieldType[x],
        F, Fn,
        S, Sn,
        V, Vn,
        G, Gn,
        A, An];

GetFieldType[x_] := SARAH`getType[x, False, FlexibleSUSY`FSEigenstates];

EnforceCpColorStructures[nPointFunctions_List] :=
    EnforceCpColorStructures /@ nPointFunctions;

EnforceCpColorStructures[
    SelfEnergies`FSSelfEnergy[p : f_Symbol[__] | f_Symbol, expr__]] :=
    SelfEnergies`FSSelfEnergy@@Prepend[EnforceCpColorStructures[f, {expr}], p];

EnforceCpColorStructures[nPointFunction_] := nPointFunction;

EnforceCpColorStructures[externalField_, exprs_List] := Fold[
    EnforceCpColorStructure[externalField, #1, #2] &,
    exprs,
    Select[Cases[exprs, _SARAH`Cp|_SARAH`Cp[_], Infinity],
           !UnresolvedColorFactorFreeQ[#, exprs] &]
];

EnforceCpColorStructure[extField_, exprs_List, cp:SARAH`Cp[fields__]] :=
    CpColorStructureFunction[extField, {fields}, cp, SARAH`Cp@@# &] @
    exprs;

EnforceCpColorStructure[extField_, exprs_List, cp:SARAH`Cp[fields__][lor_]] :=
    CpColorStructureFunction[extField, {fields}, cp, (SARAH`Cp@@#)[lor] &] @
    exprs;

CpColorStructureFunction[extField_, fields_List, cp_, wrap_] := Module[{
        reordered = PullExternalFieldsToLeft[extField, fields]
    },
    If[fields === reordered,
       Identity,
       Message[EnforceCpColorStructures::cpext, extField, cp];
       # /. cp -> wrap[reordered] &]
];

(* look for external fields preferably in unrotated basis *)

PullExternalFieldsToLeft[f_, lst:{a_?IsUnrotated, b_?IsUnrotated, ___}] /;
    FieldHead@ToRotatedField[a] === FieldHead@ToRotatedField[b] === f &&
    StripFieldIndices[a] === Susyno`LieGroups`conj@StripFieldIndices[b] :=
        lst;

PullExternalFieldsToLeft[
    f_,
    {x___, a_?IsUnrotated, y___, b_Susyno`LieGroups`conj?IsUnrotated, z___} |
    {x___, b_Susyno`LieGroups`conj?IsUnrotated, y___, a_?IsUnrotated, z___}] /;
    FieldHead@ToRotatedField[a] === FieldHead@ToRotatedField[b] === f &&
    StripFieldIndices[a] === Susyno`LieGroups`conj@StripFieldIndices[b] :=
        {a, b, x, y, z}

(* in case f is colored real scalar *)
PullExternalFieldsToLeft[
    f_, {x___, a_?IsUnrotated, y___, b_?IsUnrotated, z___}] /;
    FieldHead@ToRotatedField[a] === FieldHead@ToRotatedField[b] === f &&
    StripFieldIndices[a] === StripFieldIndices[b] :=
        {a, b, x, y, z}

(* fall back on external fields in rotated basis *)

PullExternalFieldsToLeft[f_, lst:{a_, b_, ___}] /;
    FieldHead[a] === FieldHead[b] === f &&
    StripFieldIndices[a] === Susyno`LieGroups`conj@StripFieldIndices[b] :=
        lst;

PullExternalFieldsToLeft[
    f_, {x___, a_, y___, b_Susyno`LieGroups`conj, z___} |
        {x___, b_Susyno`LieGroups`conj, y___, a_, z___}] /;
    FieldHead[a] === FieldHead[b] === f &&
    StripFieldIndices[a] === Susyno`LieGroups`conj@StripFieldIndices[b] :=
        {a, b, x, y, z};

(* in case f is colored real scalar *)
PullExternalFieldsToLeft[f_, {x___, a_, y___, b_, z___}] /;
    FieldHead[a] === FieldHead[b] === f &&
    StripFieldIndices[a] === StripFieldIndices[b] :=
        {a, b, x, y, z};

PullExternalFieldsToLeft[f_, lst_] := (
    Print["Vertices`Private`PullExternalFieldsToLeft[",
          f, ", ", lst, "] failed."];
    Abort[]
);

StripExtraFieldIndices[fields_List] := StripExtraFieldIndices /@ fields;

StripExtraFieldIndices[field_] /; !FreeQ[field, _[{}]] :=
    StripFieldIndices[field];

StripExtraFieldIndices[field_] /; !FreeQ[field, _[_?VectorQ]] &&
    SARAH`getIndizes @ FieldHead[field] === {} := StripFieldIndices[field];

StripExtraFieldIndices[field_] /; !FreeQ[field, _[{_Integer, ___}?VectorQ]] &&
    !MatchQ[SARAH`getIndizes @ FieldHead[field], {SARAH`generation, ___}] :=
    StripExtraFieldIndices[
        field /.
            head_Symbol[{_Integer, indices___}?VectorQ] :> head[{indices}]];

StripExtraFieldIndices[field_] := field;

DeleteRedundantCpPatterns[cpPatterns_] :=
    First @ Sort[#, MatchQ[#2, #1]&]& /@
    Gather[cpPatterns, MatchQ[#1, #2] || MatchQ[#2, #1]&];

VertexExp[cpPattern_, nPointFunctions_, massMatrices_] := Module[{
        cp = ToCp[cpPattern],
        rotatedCp, fieldsInRotatedCp,
        sarahVertex,
        fields, vertices,
        lorentzTag, lorentz, vertex,
        strippedIndices,
        contraction,
        factor
    },
    rotatedCp = ReplaceUnrotatedFields[cp];
    fieldsInRotatedCp = GetParticleList[rotatedCp];
    sarahVertex = SARAHVertex[fieldsInRotatedCp];
    fields = First[sarahVertex];
    vertices = Rest[sarahVertex];
    lorentzTag = GetLorentzStructure[rotatedCp];
    {vertex, lorentz} = FindVertexWithLorentzStructure[vertices, lorentzTag];
    strippedIndices = Complement[Flatten[FieldIndexList /@ fields],
                                 Flatten[FieldIndexList /@ fieldsInRotatedCp]];
    vertex = StripGroupStructure[
        ResolveColorFactor[
            vertex, fields, cpPattern, nPointFunctions[[All,2]]],
        strippedIndices];
    contraction = Block[{
            SARAH`sum
            (* corrupts a polynomial (monomial + monomial + ...) summand *)
        },
        ExpandSarahSum @ SimplifyContraction @
        InTermsOfRotatedVertex[
            vertex, lorentz,
            GetParticleList[cp], massMatrices]];
    (* see SPhenoCouplingList[] in SARAH/Package/SPheno/SPhenoCoupling.m
       for the following sign factor *)
    factor = If[GetFieldType /@ fieldsInRotatedCp === {S,S,V}, -1, 1];
    -I factor TreeMasses`ReplaceDependencies[contraction] /.
        Parameters`ApplyGUTNormalization[]
];

SARAHVertex[fieldsInRotatedCp_List] := Module[{
        sarahVertex = SARAH`Vertex @ StripFieldIndices[fieldsInRotatedCp],
        fields,
        restoreIndicesRules
    },
    Assert[MatchQ[sarahVertex, {_, __}]];
    fields = First[sarahVertex];
    Assert[StripFieldIndices[fields] === StripFieldIndices[fieldsInRotatedCp]];
    restoreIndicesRules = Flatten[
        Thread[Take[Last[#], Length @ First[#]] -> First[#]]& /@
        Transpose[{FieldIndexList /@ fieldsInRotatedCp,
                   FieldIndexList /@ fields}]];
    sarahVertex /. restoreIndicesRules
];

StripGroupStructure[expr_, indices_List] := Module[{
        indexPattern = Alternatives@@indices
    },
    expr /. {
        PatternWithDistinctIndices[SARAH`Delta, 2, indexPattern] -> 1,
        SARAH`epsTensor[__] -> 1,
        SARAH`Lam[__] -> 2,
        SARAH`Sig[__] -> 2,
        SARAH`fSU2[__] -> 1,
        SARAH`fSU3[__] -> 1,
        SARAH`CG[SARAH`SU[3], {{1,0},{0,1},{1,0},{0,1}}][a_,b_,c_,d_] -> SARAH`Delta[a,b] SARAH`Delta[c,d],
        SARAH`CG[SARAH`SU[3], {{0,1},{1,0},{0,1},{1,0}}][a_,b_,c_,d_] -> SARAH`Delta[a,b] SARAH`Delta[c,d],
        SARAH`Generator[SARAH`SU[3],___][___] -> 2
    }
];

PatternWithDistinctIndices[head_, nIndices_, indexPattern_] :=
ReleaseHold[Block[{
        patternNames = Table[Unique[], {nIndices}],
        UnsameQ
    },
    With[{condition = UnsameQ @@ patternNames},
         head @@ (Pattern[#, indexPattern]& /@ patternNames) /;
         Hold[condition]]]];

FindVertexWithLorentzStructure[vertices_List, lorentz_Integer] :=
    vertices[[lorentz]];

FindVertexWithLorentzStructure[vertices_List, lorentz_] :=
    SingleCase[vertices, {_, str_ /; !FreeQ[str, lorentz]}];

simplifyContractionDispatch = Dispatch[{
    HoldPattern[
        SARAH`sum[i_, 1, d_, x_ z_?UnitaryMatrixQ[i_, j_]]
        /; HasContractionQ[x, Susyno`LieGroups`conj[z[i, _]]] &&
        SARAH`getDim[z] === d] :>
    (x /. Susyno`LieGroups`conj[z[i, k_]] :> SARAH`Delta[j, k]),

    HoldPattern[
        SARAH`sum[i_,1,d_, x_ Susyno`LieGroups`conj[z_?UnitaryMatrixQ[i_,j_]]]
        /; HasContractionQ[x, z[i, _]] &&
        SARAH`getDim[z] === d] :>
    (x /. z[i, k_] :> SARAH`Delta[j, k]),

    HoldPattern[
        SARAH`sum[i_, 1, d_, x_ z_?UnitaryMatrixQ[j_, i_]]
        /; HasContractionQ[x, Susyno`LieGroups`conj[z[_, i]]] &&
        SARAH`getDim[z] === d] :>
    (x /. Susyno`LieGroups`conj[z[k_, i]] :> SARAH`Delta[j, k]),

    HoldPattern[
        SARAH`sum[i_,1,d_, x_ Susyno`LieGroups`conj[z_?UnitaryMatrixQ[j_,i_]]]
        /; HasContractionQ[x, z[_, i]] &&
        SARAH`getDim[z] === d] :>
    (x /. z[k_, i] :> SARAH`Delta[j, k]),

    s:HoldPattern[
    (SARAH`sum[i_, 1, d_, x_ z_?UnitaryMatrixQ[i_, j_]]
     /; HasMixedContractionQ[x, Susyno`LieGroups`conj[z[i, _]]]) |
    (SARAH`sum[i_, 1, d_, x_ Susyno`LieGroups`conj[z_?UnitaryMatrixQ[i_, j_]]]
     /; HasMixedContractionQ[x, z[i, _]]) |
    (SARAH`sum[i_, 1, d_, x_ z_?UnitaryMatrixQ[j_, i_]]
     /; HasMixedContractionQ[x, Susyno`LieGroups`conj[z[_, i]]]) |
    (SARAH`sum[i_, 1, d_, x_ Susyno`LieGroups`conj[z_?UnitaryMatrixQ[j_, i_]]]
     /; HasMixedContractionQ[x, z[_, i]]) /;
    SARAH`getDim[z] === d] :> ExpandSarahSum[s]
}];

SimplifyContraction[expr_] := expr //. simplifyContractionDispatch;

HasContractionQ[x_, form_] :=
    !MatchQ[FindContraction[x, form], Unindexed|Indexed|Mixed];

HasMixedContractionQ[x_, form_] := FindContraction[x, form] === Mixed;

FindContraction[x_, form:_?UnitaryMatrixQ[___, i_Symbol,___]] :=
    FindContraction[x, form, i];

FindContraction[x_, form:HoldPattern@Susyno`LieGroups`conj
                [_?UnitaryMatrixQ[___, i_Symbol, ___]]] :=
    FindContraction[x, form, i];

FindContraction[x_Times, form_, index_] := Module[{
        structure = FindContraction[#, form, index]& /@ List@@x,
        match
    },
    Which[MatchQ[structure, {Unindexed..}], Unindexed,
          (match =
           Replace[structure,
                   {{Unindexed... , p:Except[Unindexed], Unindexed...} -> p,
                    _ -> False}]) =!= False, match,
          MatchQ[structure, {___, Mixed, ___}], Mixed,
          MatchQ[structure, {___, Except[Unindexed], ___}], Indexed,
          True, Print["Vertices`Private`FindContraction[",
                      x, ", ", form, ", ", index, "] failed."]; Abort[]]
];

FindContraction[x_Plus, form_, index_] := Module[{
        structure = FindContraction[#, form, index]& /@ List@@x,
        match
    },
    Which[MatchQ[structure, {Unindexed..}], Unindexed,
          MatchQ[structure, {(f:Except[Unindexed])..}], First[structure],
          MatchQ[structure, {___, Except[Unindexed], ___}], Mixed,
          True, Print["Vertices`Private`FindContraction[",
                      x, ", ", form, ", ", index, "] failed."]; Abort[]]
];

FindContraction[HoldPattern@SARAH`sum[_, _, _, x_], form_, index_] :=
    FindContraction[x, form, index];

FindContraction[x_, form_, index_] /; FreeQ[x, index] := Unindexed;

FindContraction[x_, form_, index_] /; MatchQ[x, form] := x;

FindContraction[x_, form_, index_] := Indexed;

ExpandSarahSum[expr_] := expr //.
    SARAH`sum[a_, b_, c_, x_] /; Head[Expand[x]] === Plus :>
    (SARAH`sum[a, b, c, #]& /@ Expand[x]);

InTermsOfRotatedVertex[vertex_, lorentz_, uFields_List, massMatrices_] :=
Block[{
        SARAH`bar
    },
    Fold[If[IsUnrotated[#2],
            RewriteUnrotatedField[
                #1, lorentz, #2,
                SingleCase[massMatrices,
                           _[_,FieldHead@ToRotatedField[#2],z_] :> z]],
            #1]&,
         vertex, RestoreBarOnMajorana[uFields, lorentz]]
];

RestoreBarOnMajorana[uFields_List, lorentz_] :=
    RestoreBarOnMajoranaUntil[2, uFields];

RestoreBarOnMajorana[uFields_List, LorentzProduct[gamma[_], PL|PR]] :=
    RestoreBarOnMajoranaUntil[1, uFields];

RestoreBarOnMajoranaUntil[lastPos_Integer, uFields_List] :=
MapIndexed[If[First[#2] <= lastPos && MajoranaQ@FieldHead@ToRotatedField[#1],
              SARAH`bar[#1], #1]&,
           uFields];

RewriteUnrotatedField[
    expr_, _,
    uField:_Symbol[{__}], z_Symbol] :=
    ContractMixingMatrix[expr, uField, z, Identity];

RewriteUnrotatedField[
    expr_, _,
    uField:HoldPattern@Susyno`LieGroups`conj[_Symbol[{__}]], z_Symbol] :=
    ContractMixingMatrix[expr, uField, z, Susyno`LieGroups`conj];

RewriteUnrotatedField[
    expr_, LorentzProduct[gamma[_], PL],
    uField:_Symbol[{__}], {u_, v_}] :=
    ContractMixingMatrix[expr, uField, u, Identity];

RewriteUnrotatedField[
    expr_, LorentzProduct[gamma[_], PL],
    uField:HoldPattern@SARAH`bar[_Symbol[{__}]], {u_, v_}] :=
    ContractMixingMatrix[expr, uField, u, Susyno`LieGroups`conj];

RewriteUnrotatedField[
    expr_, LorentzProduct[gamma[_], PR],
    uField:_Symbol[{__}], {u_, v_}] :=
    ContractMixingMatrix[expr, uField, v, Susyno`LieGroups`conj];

RewriteUnrotatedField[
    expr_, LorentzProduct[gamma[_], PR],
    uField:HoldPattern@SARAH`bar[_Symbol[{__}]], {u_, v_}] :=
    ContractMixingMatrix[expr, uField, v, Identity];

RewriteUnrotatedField[
    expr_, PL,
    uField:_Symbol[{__}], {u_, v_}] :=
    ContractMixingMatrix[expr, uField, u, Identity];

RewriteUnrotatedField[
    expr_, PL,
    uField:HoldPattern@SARAH`bar[_Symbol[{__}]], {u_, v_}] :=
    ContractMixingMatrix[expr, uField, v, Identity];

RewriteUnrotatedField[
    expr_, PR,
    uField:_Symbol[{__}], {u_, v_}] :=
    ContractMixingMatrix[expr, uField, v, Susyno`LieGroups`conj];

RewriteUnrotatedField[
    expr_, PR,
    uField:HoldPattern@SARAH`bar[_Symbol[{__}]], {u_, v_}] :=
    ContractMixingMatrix[expr, uField, u, Susyno`LieGroups`conj];

ContractMixingMatrix[expr_, uField_, z_, op_] := Module[{
        uIndices = FieldIndexList[uField],
        field = ToRotatedField[uField],
        uIndex, i = Unique["gl"]
    },
    (* CHECK: does the first index always denote the flavour? *)
    uIndex = First[uIndices];
    SARAH`sum[i, SARAH`getGenStart[field], SARAH`getGen[field],
              (expr /. uIndex -> i) op @ z[i, uIndex]]
];

FieldHead[Susyno`LieGroups`conj[field_]] := FieldHead[field];

FieldHead[SARAH`bar[field_]] := FieldHead[field];

FieldHead[field_Symbol[{__}]] := field;

FieldHead[field_Symbol] := field;

ToCpPattern[cp : _SARAH`Cp|_SARAH`Cp[_]] := cp /.
    ((# -> (# /. Thread[(# -> If[Head[#] === Symbol, Pattern[#, _], #])& /@
                        FieldIndexList[#]]))& /@
     GetParticleList[cp]);

ToCp[cpPattern : _SARAH`Cp|_SARAH`Cp[_]] := cpPattern /. p_Pattern :> First[p];

CpType[cp : _SARAH`Cp|_SARAH`Cp[_]] := RotatedCpType @
    ReplaceUnrotatedFields[cp];

RotatedCpType[SARAH`Cp[fields__] | SARAH`Cp[fields__][_]] :=
    SARAH`VType @@ GetFieldType /@ {fields};

RenumberCpIndices[SARAH`Cp[fields__]] :=
    SARAH`Cp @@ RenumberFieldIndices[{fields}]

RenumberCpIndices[SARAH`Cp[fields__][l_]] :=
    (SARAH`Cp @@ RenumberFieldIndices[{fields}])[l]

RenumberFieldIndices[fields_List] := Block[{
        UsedSarahIndexQ
    },
    UsedSarahIndexQ[_] := False;
    RenumberFieldIndices /@ fields
];

RenumberFieldIndices[field_] :=
    field /. ((# -> RenumberSarahIndex[#])& /@ FieldIndexList[field]);

RenumberSarahIndex[index_Symbol?UsedSarahIndexQ] :=
    RenumberSarahIndex @ Symbol @
        StringReplace[ToString[index],
                      RegularExpression["[[:digit:]]+$"] :>
                      ToString[ToExpression["$0"]+1]];

RenumberSarahIndex[index_Symbol] := (UsedSarahIndexQ[index] = True; index);

RenumberSarahIndex[index_] := index;

ResolveColorFactor[vertex_, fields_, cpPattern_, exprs_] /;
    UnresolvedColorFactorFreeQ[cpPattern, exprs] := vertex;

(* see SARAH`sumOverNonAbelianIndizes[] *)
ResolveColorFactor[vertex_, fields_, cpPattern_, exprs_] := Module[{
        loopArgs,
        internalColorIndices = InternalColorIndices[fields],
        externalColorIndices = ExternalColorIndices[fields],
        sumIdx, sumRange
    },
    Assert[Length[externalColorIndices] === 2];
    (* Q: does one need to sum also over external color indices?
       A: it is a way to strip the color structure of this class
       of vertices *)
    Length[internalColorIndices] === 2 || (
        Print["Vertices`Private`ResolveColorFactor[",
              vertex, ", ", fields, ", ", cpPattern, ", ", exprs, "] failed."];
        Abort[]);
    sumRange = ColorIndexRange[First @ internalColorIndices, fields];
    loopArgs = Prepend[{#, 1}& /@ externalColorIndices, {sumIdx, sumRange}];
    Sum @@
    Prepend[loopArgs,
            vertex /. Alternatives@@internalColorIndices -> sumIdx] //.
      SARAH`sum[a_,b_,c_,x_] /; !FreeQ[x, (SARAH`fSU3|SARAH`Lam)[___,a,___]] :>
      Sum[x, {a, b, c}]
];

InternalColorIndices[fields_List] :=
    Union@Cases[Drop[FieldIndexList /@ fields, 2],
                _?SarahColorIndexQ, Infinity];

ExternalColorIndices[fields_List] :=
    Union@Cases[Take[FieldIndexList /@ fields, 2],
                _?SarahColorIndexQ, Infinity];

FieldIndexList[field_] := Flatten@Cases[field, _?VectorQ, {0, Infinity}];

StripFieldIndices[field_] := field /. head_[_?VectorQ] :> head;

StripLorentzIndices[p_Symbol] := p;
StripLorentzIndices[SARAH`bar[p_]] := SARAH`bar[StripLorentzIndices[p]];
StripLorentzIndices[Susyno`LieGroups`conj[p_]] := Susyno`LieGroups`conj[StripLorentzIndices[p]];
StripLorentzIndices[p_] :=
    Module[{remainingIndices},
           remainingIndices = Select[p[[1]], (!Parameters`IsLorentzIndex[#] &)];
           If[Length[remainingIndices] === 0,
              Head[p],
              Head[p][remainingIndices]
             ]
          ];

ColorIndexRange[colorIndex_, fields_] :=
    SingleCase[
        SARAH`getIndizesWI[First @ Select[fields, !FreeQ[#, colorIndex] &]],
        {color, n_} :> n];

UnresolvedColorFactorFreeQ[cpPattern_, exprs_] := Module[{
        fstPos = First@Position[exprs, cpPattern],
        cpInstance,
        exprInstance
    },
    cpInstance = Extract[exprs, fstPos];
    exprInstance = Extract[exprs, Take[fstPos, 1]];
    FreeQ[Coefficient[exprInstance //. SARAH`sum[__, ex_] :> ex, cpInstance],
          C]
];

(* CHECK: are the following right semantics of SARAH indices? *)
SarahExternalGenerationIndexQ[index_Symbol] :=
    StringMatchQ[ToString[index], RegularExpression["gO[[:digit:]]+"]];

SarahInternalGenerationIndexQ[index_Symbol] :=
    StringMatchQ[ToString[index], RegularExpression["gI[[:digit:]]+"]];

SarahColorIndexQ[index_Symbol] :=
    StringMatchQ[ToString[index], RegularExpression["ct[[:digit:]]+"]];

GetLorentzStructure[SARAH`Cp[__]] := 1;

GetLorentzStructure[SARAH`Cp[__][a_]] := a;

GetParticleList[SARAH`Cp[a__]] := {a};

GetParticleList[SARAH`Cp[a__][_]] := {a};

IsUnrotated[SARAH`bar[field_]] := IsUnrotated[field];

IsUnrotated[Susyno`LieGroups`conj[field_]] := IsUnrotated[field];

IsUnrotated[field_[__]] := IsUnrotated[field];

IsUnrotated[field_Symbol] := StringTake[ToString[field],1] === "U";

ToRotatedField[field_Symbol] :=
    Symbol[StringReplace[ToString[field], StartOfString ~~ "U" ~~ rest_ :> rest]];

ToRotatedField[SARAH`bar[field_]] := SARAH`bar[ToRotatedField[field]];

ToRotatedField[Susyno`LieGroups`conj[field_]] := Susyno`LieGroups`conj[ToRotatedField[field]];

ToRotatedField[field_List] := ToRotatedField /@ field;

ToRotatedField[field_[indices__]] := ToRotatedField[field][indices];

ReplaceUnrotatedFields[SARAH`Cp[p__]] :=
    SARAH`Cp @@ ToRotatedField[{p}];

ReplaceUnrotatedFields[SARAH`Cp[p__][lorentz_]] :=
    ReplaceUnrotatedFields[SARAH`Cp[p]][lorentz];

LoadVerticesIfNecessary[] :=
    If[Head[SARAH`VertexList3] === Symbol || Length[SARAH`VertexList3] === 0,
       SA`CurrentStates = FlexibleSUSY`FSEigenstates;
       SARAH`InitVertexCalculation[FlexibleSUSY`FSEigenstates, False];
       SARAH`partDefinition = ParticleDefinitions[FlexibleSUSY`FSEigenstates];
       SARAH`Particles[SARAH`Current] = SARAH`Particles[FlexibleSUSY`FSEigenstates];
       SARAH`ReadVertexList[FlexibleSUSY`FSEigenstates, False, False, True];
       SARAH`MakeCouplingLists;
      ];

CreateVertexData[fields_List, fieldsNamespace_:""] :=
    Module[{dataClassName,indexBounds,parsedVertex},
           dataClassName = "VertexData<" <> StringJoin[Riffle[
               TreeMasses`CreateFieldClassName[#, prefixNamespace -> fieldsNamespace] & /@ fields,
               ", "]] <> ">";

           "template<> struct " <> dataClassName <> "\n" <> "{\n" <>
           TextFormatting`IndentText[
               "using vertex_type = " <> SymbolName[VertexTypeForFields[fields]] <>
               ";"] <> "\n" <> "};"
          ];

(* Returns the necessary c++ code corresponding to the vertices that need to be calculated. *)
CreateVertices[vertices_List, fieldsNamespace_:""] :=
    StringJoin @\[NonBreakingSpace]Riffle[CreateVertex[#, fieldsNamespace] & /@ DeleteDuplicates[vertices], "\n\n"]

(* Creates the actual c++ code for a vertex with given fields.
 You should never need to change this code! *)
CreateVertex[fields_List, fieldsNamespace_:""] :=
    Module[{parsedVertex, functionClassName},
           LoadVerticesIfNecessary[];
           functionClassName = "Vertex<" <> StringJoin @ Riffle[
           TreeMasses`CreateFieldClassName[#, prefixNamespace -> fieldsNamespace] & /@ fields, ", "] <> ">";

           "template<> template <class EvaluationContext> inline\n" <>
           functionClassName <> "::vertex_type\n" <>
           functionClassName <> "::evaluate(const indices_type& indices, const EvaluationContext& context)\n" <>
           "{\n" <>
           TextFormatting`IndentText @ VertexFunctionBodyForFields[fields] <> "\n" <>
           "}"
          ];

VertexFunctionBodyForFields[fields_List] :=
    Switch[Length[fields],
           3, VertexFunctionBodyForFieldsImpl[fields, SARAH`VertexList3],
           4, VertexFunctionBodyForFieldsImpl[fields, SARAH`VertexList4],
           _, "non-(3,4)-point vertex"
          ];

GetIndexedFieldsForVertex[fields_, vertex_] :=
    Module[{sortedIndexedFields, sortedFields,
            fieldsOrdering, sortedFieldsOrdering,
            inverseFOrdering, fOrderingWRTSortedF},
           sortedIndexedFields = vertex[[1]];
           sortedFields = StripFieldIndices[sortedIndexedFields];

           (* Mathematica 7 does not know about permutations... :'-( *)
           fieldsOrdering = Ordering[fields];
           sortedFieldsOrdering = Ordering[sortedFields];

           inverseFOrdering = Ordering[fieldsOrdering];
           fOrderingWRTSortedF = sortedFieldsOrdering[[inverseFOrdering]];

           sortedIndexedFields[[fOrderingWRTSortedF]]
          ];

CreateZeroVertex[fields_?IsSSSVertex] := "return vertex_type(0);";
CreateZeroVertex[fields_?IsFFSVertex] := "return vertex_type(0, 0);";
CreateZeroVertex[fields_?IsSSVVertex] :=
    "return vertex_type(0, " <>
    StringJoin[Riffle[
        ToString /@ Flatten[Position[fields,
                                     field_ /; TreeMasses`IsScalar[field] || TreeMasses`IsGhost[field],
                                     {1}, Heads -> False] - 1], ", "]] <> ");";
CreateZeroVertex[fields_?IsSVVVertex] := "return vertex_type(0);";

(* Creates local declarations of field indices, whose values are taken
   from the elements of `arrayName'.
 *)
DeclareIndices[indexedFields_List, arrayName_String] :=
    Module[{p, total = 0, fieldIndexList, decl = ""},
           DeclareIndex[idx_, num_Integer, an_String] := (
               "const int " <> CConversion`ToValidCSymbolString[idx] <>
               " = " <> an <> "[" <> ToString[num] <> "];\n");
           For[p = 1, p <= Length[indexedFields], p++,
               fieldIndexList = FieldIndexList[indexedFields[[p]]];
               decl = decl <> StringJoin[DeclareIndex[#, total++, arrayName]& /@ fieldIndexList];
              ];
           Assert[total == Total[Length[FieldIndexList[#]]& /@ indexedFields]];
           decl
          ];

(* @todo implement or remove this! *)
SarahToFSVertexConventions[sortedFields_, expr_] := expr;

CreateScalarVertexFunctionBody[fields_, vertex_, stripGroupStructureRules_] :=
    Module[{sortedIndexedFields, sortedFields, indexedFields,
            vertexRules, expr, resultType},
           sortedIndexedFields = vertex[[1]];
           sortedFields = StripFieldIndices[sortedIndexedFields];
           indexedFields = GetIndexedFieldsForVertex[fields, vertex];

           vertexRules = {(SARAH`Cp @@ sortedIndexedFields) ->
                          FindVertexWithLorentzStructure[Rest[vertex], 1][[1]]};

           expr = CanonicalizeCoupling[SARAH`Cp @@ fields,
                  sortedFields, sortedIndexedFields] /. vertexRules /. stripGroupStructureRules;

           expr = SarahToFSVertexConventions[sortedFields, expr];
           expr = TreeMasses`ReplaceDependenciesReverse[expr];

           resultType = CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]];

           DeclareIndices[StripLorentzIndices /@ indexedFields, "indices"] <>
           Parameters`CreateLocalConstRefs[expr] <> "\n" <>
           "const " <> resultType <> " result = " <>
           Parameters`ExpressionToString[expr] <> ";\n\n" <>
           "return vertex_type(result);"
          ];

CreateChiralVertexFunctionBody[fields_, vertex_, stripGroupStructureRules_] :=
    Module[{sortedIndexedFields, sortedFields, indexedFields,
            vertexRules, exprL, exprR, resultTypeL, resultTypeR},
           sortedIndexedFields = vertex[[1]];
           sortedFields = StripFieldIndices[sortedIndexedFields];
           indexedFields = GetIndexedFieldsForVertex[fields, vertex];

           vertexRules = {
                          (SARAH`Cp @@ sortedIndexedFields)[SARAH`PL] ->
                          FindVertexWithLorentzStructure[Rest[vertex], SARAH`PL][[1]],
                          (SARAH`Cp @@ sortedIndexedFields)[SARAH`PR] ->
                          FindVertexWithLorentzStructure[Rest[vertex], SARAH`PR][[1]]
                         };

           exprL = CanonicalizeCoupling[(SARAH`Cp @@ fields)[SARAH`PL],
                                        sortedFields, sortedIndexedFields] /. vertexRules /. stripGroupStructureRules;
           exprR = CanonicalizeCoupling[(SARAH`Cp @@ fields)[SARAH`PR],
                                        sortedFields, sortedIndexedFields] /. vertexRules /. stripGroupStructureRules;

           exprL = SarahToFSVertexConventions[sortedFields, exprL];
           exprR = SarahToFSVertexConventions[sortedFields, exprR];
           exprL = TreeMasses`ReplaceDependenciesReverse[exprL];
           exprR = TreeMasses`ReplaceDependenciesReverse[exprR];

           resultTypeL = CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]];
           resultTypeR = CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]];

           DeclareIndices[StripLorentzIndices /@ indexedFields, "indices"] <>
           Parameters`CreateLocalConstRefs[{exprL, exprR}] <> "\n" <>
           "const " <> resultTypeL <> " left = " <>
           Parameters`ExpressionToString[exprL] <> ";\n\n" <>
           "const " <> resultTypeR <> " right = " <>
           Parameters`ExpressionToString[exprR] <> ";\n\n" <>
           "return vertex_type(left, right);"
          ];

CreateSSVVertexFunctionBody[fields_, vertex_, stripGroupStructureRules_] :=
    Module[{sortedIndexedFields, sortedFields, indexedFields,
            incomingScalar, outgoingScalar, vertexRules,
            expr, resultType},
           sortedIndexedFields = vertex[[1]];
           sortedFields = StripFieldIndices[sortedIndexedFields];
           indexedFields = GetIndexedFieldsForVertex[fields, vertex];

           {incomingScalar, outgoingScalar} = Replace[vertex[[2,2]],
                                                      SARAH`Mom[is_,_] - SARAH`Mom[os_,_] :> {is, os}];
           vertexRules = {(SARAH`Cp @@ sortedIndexedFields) -> vertex[[2,1]]};

           expr = CanonicalizeCoupling[SARAH`Cp @@ fields,
                                       sortedFields, sortedIndexedFields] /. vertexRules /. stripGroupStructureRules;

           expr = SarahToFSVertexConventions[sortedFields, expr];
           expr = TreeMasses`ReplaceDependenciesReverse[expr];

           resultType = CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]];

           "int minuend_index = " <>
           ToString[Position[indexedFields, incomingScalar][[1,1]] - 1] <> ";\n" <>
           "int subtrahend_index = " <>
           ToString[Position[indexedFields, outgoingScalar][[1,1]] - 1] <> ";\n\n" <>
           DeclareIndices[StripLorentzIndices /@ indexedFields, "indices"] <>
           Parameters`CreateLocalConstRefs[expr] <> "\n" <>
           "const " <> resultType <> " result = " <>
           Parameters`ExpressionToString[expr] <> ";\n\n" <>
           "return vertex_type(result, minuend_index, subtrahend_index);"
          ];

CreateSSVVVertexFunctionBody[fields_, vertex_, stripGroupStructureRules_] := "";

CreateSVVVertexFunctionBody[fields_, vertex_, stripGroupStructureRules_] :=
    Module[{sortedIndexedFields, sortedFields, indexedFields,
            vertexRules, expr, resultType},
           sortedIndexedFields = vertex[[1]];
           sortedFields = StripFieldIndices[sortedIndexedFields];
           indexedFields = GetIndexedFieldsForVertex[fields, vertex];

           vertexRules = {(SARAH`Cp @@ sortedIndexedFields) -> vertex[[2,1]]};

           expr = CanonicalizeCoupling[SARAH`Cp @@ fields,
                                       sortedFields, sortedIndexedFields] /. vertexRules /. stripGroupStructureRules;

           expr = SarahToFSVertexConventions[sortedFields, expr];
           expr = TreeMasses`ReplaceDependenciesReverse[expr];

           resultType = CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]];

           DeclareIndices[StripLorentzIndices /@ indexedFields, "indices"] <>
           Parameters`CreateLocalConstRefs[expr] <> "\n" <>
           "const " <> resultType <> " result = " <>
           Parameters`ExpressionToString[expr] <> ";\n\n" <>
           "return vertex_type(result);"
          ];

CreateVVVVertexFunctionBody[fields_, vertex_, stripGroupStructureRules_] := "";

CreateVVVVVertexFunctionBody[fields_, vertex_, stripGroupStructureRules_] :=
    Module[{sortedIndexedFields, sortedFields, indexedFields,
            vertexRules, expr1234, expr1324, expr1423,
            resultType1234, resultType1324, resultType1423},
           sortedIndexedFields = vertex[[1]];
           sortedFields = StripFieldIndices[sortedIndexedFields];
           indexedFields = GetIndexedFieldsForVertex[fields, vertex];

           vertexRules = {
                          (SARAH`Cp @@ sortedIndexedFields)[SARAH`g[SARAH`lt1, SARAH`lt2] SARAH`g[SARAH`lt3, SARAH`lt4]]
                              -> FindVertexWithLorentzStructure[Rest[vertex], SARAH`g[SARAH`lt1, SARAH`lt2] SARAH`g[SARAH`lt3, SARAH`lt4]][[1]],
                          (SARAH`Cp @@ sortedIndexedFields)[SARAH`g[SARAH`lt1, SARAH`lt3] SARAH`g[SARAH`lt2, SARAH`lt4]]
                              -> FindVertexWithLorentzStructure[Rest[vertex], SARAH`g[SARAH`lt1, SARAH`lt3] SARAH`g[SARAH`lt2, SARAH`lt4]][[1]],
                          (SARAH`Cp @@ sortedIndexedFields)[SARAH`g[SARAH`lt1, SARAH`lt4] SARAH`g[SARAH`lt2, SARAH`lt3]]
                              -> FindVertexWithLorentzStructure[Rest[vertex], SARAH`g[SARAH`lt1, SARAH`lt4] SARAH`g[SARAH`lt2, SARAH`lt3]][[1]]
                         };

           expr1234 = CanonicalizeCoupling[(SARAH`Cp @@ fields)[SARAH`g[SARAH`lt1, SARAH`lt2] SARAH`g[SARAH`lt3, SARAH`lt4]],
                                           sortedFields, sortedIndexedFields] /. vertexRules /. stripGroupStructureRules;
           expr1324 = CanonicalizeCoupling[(SARAH`Cp @@ fields)[SARAH`g[SARAH`lt1, SARAH`lt3] SARAH`g[SARAH`lt2, SARAH`lt4]],
                                           sortedFields, sortedIndexedFields] /. vertexRules /. stripGroupStructureRules;
           expr1423 = CanonicalizeCoupling[(SARAH`Cp @@ fields)[SARAH`g[SARAH`lt1, SARAH`lt4] SARAH`g[SARAH`lt2, SARAH`lt3]],
                                           sortedFields, sortedIndexedFields] /. vertexRules /. stripGroupStructureRules;

           expr1234 = SarahToFSVertexConventions[sortedFields, expr1234];
           expr1324 = SarahToFSVertexConventions[sortedFields, expr1324];
           expr1423 = SarahToFSVertexConventions[sortedFields, expr1423];
           expr1234 = TreeMasses`ReplaceDependenciesReverse[expr1234];
           expr1324 = TreeMasses`ReplaceDependenciesReverse[expr1324];
           expr1423 = TreeMasses`ReplaceDependenciesReverse[expr1423];

           resultType1234 = CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]];
           resultType1324 = CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]];
           resultType1423 = CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]];

           "const int i = " <>
           ToString[Position[indexedFields, sortedIndexedFields[[1]]][[1,1]] - 1] <> ";\n" <>
           "const int j = " <>
           ToString[Position[indexedFields, sortedIndexedFields[[2]]][[1,1]] - 1] <> ";\n" <>
           "const int k = " <>
           ToString[Position[indexedFields, sortedIndexedFields[[3]]][[1,1]] - 1] <> ";\n" <>
           "const int l = " <>
           ToString[Position[indexedFields, sortedIndexedFields[[4]]][[1,1]] - 1] <> ";\n" <>
           DeclareIndices[StripLorentzIndices /@ indexedFields, "indices"] <>
           Parameters`CreateLocalConstRefs[{expr1234, expr1324, expr1423}] <> "\n" <>
           "const " <> resultType1234 <> " c1 = " <>
           Parameters`ExpressionToString[expr1234] <> ";\n\n" <>
           "const " <> resultType1324 <> " c2 = " <>
           Parameters`ExpressionToString[expr1324] <> ";\n\n" <>
           "const " <> resultType1423 <> " c3 = " <>
           Parameters`ExpressionToString[expr1423] <> ";\n\n" <>
           "return vertex_type(i, j, k, l, c1, c2, c3);"
          ];

CreateGGVVertexFunctionBody[fields_, vertex_, stripGroupStructureRules_] :=
    Module[{sortedIndexedFields, sortedFields, indexedFields,
            vertexRules, expr, resultType},
           sortedIndexedFields = vertex[[1]];
           sortedFields = StripFieldIndices[sortedIndexedFields];
           indexedFields = GetIndexedFieldsForVertex[fields, vertex];

           vertexRules = {(SARAH`Cp @@ sortedIndexedFields) ->
                          FindVertexWithLorentzStructure[Rest[vertex], SARAH`Mom[_, _]][[1]]};

           expr = CanonicalizeCoupling[SARAH`Cp @@ fields,
                                       sortedFields, sortedIndexedFields] /. vertexRules /. stripGroupStructureRules;

           expr = SarahToFSVertexConventions[sortedFields, expr];
           expr = TreeMasses`ReplaceDependenciesReverse[expr];

           resultType = CConversion`CreateCType[CConversion`ScalarType[CConversion`complexScalarCType]];

           DeclareIndices[StripLorentzIndices /@ indexedFields, "indices"] <>
           Parameters`CreateLocalConstRefs[expr] <> "\n" <>
           "const " <> resultType <> " result = " <>
           Parameters`ExpressionToString[expr] <> ";\n\n" <>
           "return vertex_type(result);"
          ];

VertexFunctionBodyForFieldsImpl[fields_List, vertexList_List] :=
    Module[{sortedFields,
            vertex, vertexType = VertexTypeForFields[fields],
            stripGroupStructure = {SARAH`Lam[__] -> 2, SARAH`fSU3[__] -> 1}},
           sortedFields = SortFieldsInCp[fields];

           vertex = Select[vertexList, StripFieldIndices[#[[1]]] === sortedFields &, 1];

           If[vertex === {},
              Return[CreateZeroVertex[fields]];
             ];

           vertex = vertex[[1]];

           Switch[vertexType,
                  SSSVertex, CreateScalarVertexFunctionBody[fields, vertex, stripGroupStructure],
                  SSSSVertex, CreateScalarVertexFunctionBody[fields, vertex, stripGroupStructure],
                  FFSVertex, CreateChiralVertexFunctionBody[fields, vertex, stripGroupStructure],
                  FFVVertex, CreateChiralVertexFunctionBody[fields, vertex, stripGroupStructure],
                  SSVVertex, CreateSSVVertexFunctionBody[fields, vertex, stripGroupStructure],
                  SSVVVertex, CreateSSVVVertexFunctionBody[fields, vertex, stripGroupStructure],
                  SVVVertex, CreateSVVVertexFunctionBody[fields, vertex, stripGroupStructure],
                  VVVVertex, CreateVVVVertexFunctionBody[fields, vertex, stripGroupStructure],
                  VVVVVertex, CreateVVVVVertexFunctionBody[fields, vertex, stripGroupStructure],
                  GGSVertex, CreateScalarVertexFunctionBody[fields, vertex, stripGroupStructure],
                  GGVVertex, CreateGGVVertexFunctionBody[fields, vertex, stripGroupStructure]
                 ]
          ];

(* Get a mathematical expression of the requested vertex
   in terms of its canonically ordered vertex. *)
CanonicalizeCoupling[coupling_, sortedFields_List, sortedIndexedFields_List] :=
    Module[{expr},
           (* Note: Cannot use indexed fields, because their internal
              SARAH ordering is different... *)
           expr = Vertices`SortCp[coupling];

           (* Index the fields *)
           expr = expr /. {SARAH`Cp @@ sortedFields -> SARAH`Cp @@ sortedIndexedFields}
          ];

End[] (* `Private` *)

EndPackage[]
