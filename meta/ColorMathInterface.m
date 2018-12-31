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

(* This is a FlexibleSUSY interface to ColorMath package by Malin Sj̈odahl
   [http://inspirehep.net/record/1201957] *)

BeginPackage["ColorMathInterface`",
   {"SARAH`", "TreeMasses`", (*IsColorIndex is there*)"Parameters`", "ColorMath`"}
];

FSCalcColorFactor::usage = "";
SARAHToColorMathSymbols::usage = "";
GetFieldColorIndex::usage = "";
ColorN::usage = "Evaluate numerically the analytic form of the color factor.";
ConnectColorLines::usage = "";
StripSU3Generators::usage = "";
SortColorDeltas::usage = "Sort ";
TakeOnlyColor::usage = "";

Begin["`Private`"];

(*  given a field will return it's indices,
    e.g. Fd[{a,b}] or bar[Fd[{a,b}] or conj[Sd[{a,b}]] will return {a,b} *)
GetFieldIndices[field_] :=
    field /. SARAH`bar | Susyno`LieGroups`conj -> Identity /. _[x_List] :> x;

(*  given a field with indices will discard indices,
    e.g. bar[Fd[{a,b}] -> bar[Fd] *)
DropFieldIndices[field_] :=
    field /. f_[_List] -> f;


GetFieldColorIndex[field_/;TreeMasses`ColorChargedQ[field]]:=
  Module[{res},
    res = GetFieldIndices[field];
    res = Select[res, IsColorIndex];
    Assert[Length[res] === 1];
    res[[1]]
  ];

(* Calculate overal color factor from a list of color connected
   SARAH'Vertex 'objects' *)
FSCalcColorFactor[vertex_List] :=
   Module[{return},
      return =
          (SARAHToColorMathSymbols /@ vertex) // DropColorles;
      If[return === {}, Return[1]];
      return =
         TakeOnlyColor @ return;
      return = Times @@ return;
      (* CSimplify[1] doesn't evaluate *)
      If[return === 1,
         1,
         CSimplify[return]
      ]
   ];

ColorStructureFreeQ[el_] :=
   FreeQ[el,
      Subscript[Superscript[Superscript[ColorMath`CMt, List[__]], _], _] |
         Superscript[ColorMath`CMf, List[__]] |
         Superscript[ColorMath`CMd, List[__]] |
         Superscript[ColorMath`CMo, List[__]] |
         Subscript[Superscript[ColorMath`CM\[Delta], _], _] |
         Superscript[ColorMath`CM\[CapitalDelta], List[_, _]]
   ];

DropColorles::notes = "Drop colorles vertices from the list of Vertex objets  ";
DropColorles[vertices_List] :=  
   DeleteCases[vertices,
      el_ /; ColorStructureFreeQ[el]
   ];


TakeOnlyColor[vvvv__] :=
    Module[{result},
      (* the generic structure of the Vertex "object" is 
         {{ParticleList},{{Coefficient 1, Lorentz 1},{Coefficient 2, Lorentz 2},...} *)
      (* drop ParticleList *)
      result = Drop[#, 1]& /@ vvvv;
      result = (Transpose @ Drop[Transpose[#], -1])& /@ result;
      result = result //.
         __?ColorStructureFreeQ el_ /; !ColorStructureFreeQ[el] :> el //.
          (el_  /; !ColorStructureFreeQ[el]) __?ColorStructureFreeQ  :> el;
      result = DeleteCases[#, {0}]& /@ result;
      (* There should be one, overal color factor for an entire vertex.
         For example, both left and right handed parts should have it the same. *)
      Assert[CountDistinct[#] === 1]& /@ result;
      result = DeleteDuplicates /@ result;
      result = Flatten[result, 2];
      Print["Final ", result];
      result = result /. Subscript[Superscript[Superscript[ColorMath`CMt, List[a__]], b_], c_] :> 2 * Subscript[Superscript[Superscript[ColorMath`CMt, List[a]], b], c];
      Print["Final ", result];
      result
    ];

AntiFieldQ[field_] :=
    If[Head[field] === bar || Head[field] === conj, True, False];

(* Convert color symbols in a single SARAH`Vertex object to ColorMAth convention *)
SARAHToColorMathSymbols[vertex_List] :=
    Module[{result},

   result = vertex //.
      SARAH`Lam[colIdx1_, colIdx2_, colIdx3_] :> 2 ColorMath`CMt[{colIdx1}, colIdx2, colIdx3] //.
      SARAH`fSU3[colSeq__] :> ColorMath`CMf[colSeq];

   (* if the result has Delta, we need to find out if it's adj. or fundamental *)
    result = result /. SARAH`Delta[c1_?IsColorIndex, c2_?IsColorIndex] :>
        If[getColorRep[Select[vertex[[1]], ! FreeQ[#, c1] &][[1]]] === T,

           ColorMath`CMdelta @@ If[!AntiFieldQ[Select[vertex[[1]], ! FreeQ[#, c1] &][[1]]], {c2, c1}, {c1,c2}]
        ];
   result = result /. SARAH`Delta[c1_?IsColorIndex, c2_?IsColorIndex] :>
       If[getColorRep[Select[vertex[[1]], ! FreeQ[#, c1] &][[1]]] === O,
          ColorMath`CMDelta[c2, c1]
       ];

       result
   ];

(* for SU(3) *)
ColorN[expr_] :=
   expr /. ColorMath`Nc -> 3 /. ColorMath`TR -> 1/2;

(* connect color indices of field1 and field2 *)
ConnectColorLines[field1_, field2_] :=
   Module[{r1 = getColorRep[field1], r2 = getColorRep[field2]},
      Assert[r1 === r2];
      Switch[r1,
         S, 1,
         T, ColorMath`CMdelta @@ (GetFieldColorIndex /@ {field1, field2}),
         O, ColorMath`CMDelta @@ (GetFieldColorIndex /@ {field1, field2}),
         (* for now we only deal with color scalars, triplets and octets *)
         _, Abort[]
      ]
   ];


(* coefficient of the generator *)
StripSU3Generators[inP_, outP_, spec_, c_] :=
   Module[{temp, temp2, temp3},
      (* SSS *)
      If[getColorRep[inP] === S && getColorRep[outP] === S && getColorRep[spec] === S,
         Return[c]
      ];

      (* TTS *)
      If[getColorRep[inP] === T && getColorRep[outP] === T && getColorRep[spec] === S,
         temp =  Coefficient[c, ColorMath`CMdelta @@ (GetFieldColorIndex /@ {outP, inP})];
         If[temp === 0, Abort[]];
         Return[temp];
      ];

      (* TTO *)
      If[getColorRep[inP] === T && getColorRep[outP] === T && getColorRep[spec] === O,
         temp = 
            Coefficient[c, ColorMath`CMt[{GetFieldColorIndex[spec]}, GetFieldColorIndex[outP], GetFieldColorIndex[inP]]];
         If[temp === 0, Abort[]];
         Return[temp];
      ];

      (* OOO *)
      (* closed color triplet line gives
         o^abc = Tr(Ta Tb Tc) = 1/2 TR (dabc + I fabc) *)
      If[getColorRep[inP] === O && getColorRep[outP] === O && getColorRep[spec] === O,
         (* expect f^{out in V} *)
         temp =
            Coefficient[
               Expand[c] /.  a_ Superscript[ColorMath`CMo, {c1_,c2_,c3_}] +
                  b_ Superscript[ColorMath`CMo, {c1_,c3_,c2_}] /; a===-b :>
a TR Superscript[ColorMath`CMf, {c1,c2,c3}]
            (*/. Superscript[ColorMath`CMo, l_List /; Length[l]===3] :> 1/2 ColorMath`TR (
                  Superscript[ColorMath`CMd, l] + I Superscript[ColorMath`CMf, l])*),
                  Superscript[ColorMath`CMf, GetFieldColorIndex /@ {outP, inP, spec}]
            ];
         If[temp === 0,

            temp = Coefficient[
               c,
               Superscript[ColorMath`CMo, GetFieldColorIndex /@ {outP, spec, inP}]
            ];

         ];
         If[temp === 0,

            temp = Coefficient[
               c,
               Superscript[ColorMath`CMo, GetFieldColorIndex /@ {outP, inP, spec}]
            ];

         ];
            If[temp === 0, Abort[]];
            Return[temp];
      ];

      Print["Unhandled color combination of external particles in form factor"];
      Abort[];
   ];

FlipDeltaIdxRule[field1_, field2_] :=
    Module[{colIdx1 = GetFieldColorIndex[field1], colIdx2 = GetFieldColorIndex[field2],
            colRep1 = getColorRep[field1], colRep2 = getColorRep[field2], repl},

       Assert[colRep1 === colRep2];
       repl = Switch[colRep1,
          T, ColorMath`CMdelta[colIdx1, colIdx2] -> ColorMath`CMdelta[colIdx2, colIdx1],
          O, ColorMath`CMDelta[colIdx1, colIdx2] -> ColorMath`CMDelta[colIdx2, colIdx1],
          _, Abort[]
       ];
       repl
    ];

SortColorDeltas[inP_, outP_, V_, Fin_, Fout_, SIn_, Sout_] :=
    Module[{
       inPCCharge = getColorRep[inP],
       outPCCharge = getColorRep[outP],
       VCCharge = getColorRep[V],
       FCCharge = getColorRep[Fout],
       SCCharge = getColorRep[Sout],
       rule
    },
       rule = Switch[{inPCCharge, VCCharge, SCCharge, FCCharge},
          {T, S, T, S}, {
             FlipDeltaIdxRule[inP, SIn],
             FlipDeltaIdxRule[Sout, outP],
             FlipDeltaIdxRule[SIn, Sout]
          },
          {T, S, S, T}, {
             FlipDeltaIdxRule[inP, Fin],
             FlipDeltaIdxRule[Fout, outP]
          },
          (* closed color loop *)
          {S, S, T, T}, {
             FlipDeltaIdxRule[SIn, Fin],
             FlipDeltaIdxRule[Fout, Sout],
             FlipDeltaIdxRule[Sout, SIn],
             FlipDeltaIdxRule[Fin, Fout]
          },
          {T, O, S, T}, {
             FlipDeltaIdxRule[inP, Fin],
             FlipDeltaIdxRule[Fout, outP]
          },
          {T, O, T, S},  {
             FlipDeltaIdxRule[inP, SIn],
             FlipDeltaIdxRule[Sout, outP],
             (* this is fishy *)
             ColorMath`CMt[{GetFieldColorIndex[V]}, Sequence @@ (GetFieldColorIndex /@ {SIn, Sout})] ->
                 ColorMath`CMt[{GetFieldColorIndex[V]}, Sequence @@ (GetFieldColorIndex /@ {Sout, SIn})]
          },
          {T, O, T, O}, {
             ColorMath`CMt[{GetFieldColorIndex[V]}, Sequence @@ (GetFieldColorIndex /@ {SIn, Fin})] ->
                 ColorMath`CMt[{GetFieldColorIndex[V]}, Sequence @@ (GetFieldColorIndex /@ {Fin, SIn})],
             ColorMath`CMt[{GetFieldColorIndex[V]}, Sequence @@ (GetFieldColorIndex /@ {Fout, Sout})] ->
                 ColorMath`CMt[{GetFieldColorIndex[V]}, Sequence @@ (GetFieldColorIndex /@ {Sout, Fout})],
             ColorMath`CMt[{GetFieldColorIndex[V]}, Sequence @@ (GetFieldColorIndex /@ {SIn, Sout})] ->
                 ColorMath`CMt[{GetFieldColorIndex[V]}, Sequence @@ (GetFieldColorIndex /@ {Sout, SIn})],
             ColorMath`CMf[GetFieldColorIndex /@ {Fin, Fout, V}] -> ColorMath`CMf[GetFieldColorIndex /@ {Fout, Fin, V}]
          },
         {O, O, T, T}, {
             FlipDeltaIdxRule[Sout, SIn],
             FlipDeltaIdxRule[Fin, Fout],
             ColorMath`CMt[{GetFieldColorIndex[inP]}, Sequence @@ (GetFieldColorIndex /@ {SIn, Fin})] ->
                 ColorMath`CMt[{GetFieldColorIndex[inP]}, Sequence @@ (GetFieldColorIndex /@ {Fin, SIn})],
             ColorMath`CMt[{GetFieldColorIndex[outP]}, Sequence @@ (GetFieldColorIndex /@ {Fout, Sout})] ->
                 ColorMath`CMt[{GetFieldColorIndex[outP]}, Sequence @@ (GetFieldColorIndex /@ {Sout, Fout})],
             ColorMath`CMt[{GetFieldColorIndex[V]}, Sequence @@ (GetFieldColorIndex /@ {Sout, SIn})] ->
                ColorMath`CMt[{GetFieldColorIndex[V]}, Sequence @@ (GetFieldColorIndex /@ {SIn, Sout})],
             ColorMath`CMt[{GetFieldColorIndex[V]}, Sequence @@ (GetFieldColorIndex /@ {Fin, Fout})] ->
                ColorMath`CMt[{GetFieldColorIndex[V]}, Sequence @@ (GetFieldColorIndex /@ {Fout, Fin})]
         },
          _, {}
       ];
       rule
    ];

(*NumberOfExternalParticles[l_List] :=*)
   (*Count[l, Except[_List], 1];*)

(*NumberOfVertices[l_List] :=*)
   (*Length[l] - NumberOfExternalParticles[l];*)
(*
RegenerateIndices[l_List, graph_]:=
    Module[{keys, extFields, particlesInVertices,
            vertices, vertex, vertex1, vertex2,
            extField, inField,
            fieldsInVertex, fieldsInVertex1, fieldsInVertex2},

        extFields = TakeWhile[l, (Head[#]=!=List)&];
        keys = GenerateUniqueColorAssociationsForExternalParticles[l];
        particlesInVertices = Drop[l, Length@extFields];
        (* change to vertices = SARAH`Vertex /@ particlesInVertices; *)
        vertices = SARAH`Vertex[#]& /@ particlesInVertices;

        (* STEP1: set colors of external particles in vertices to values in keys *)
        (* loop over external particles *)
        For[extIdx=1, extIdx <= Length[extFields], extIdx++,
            extField = extFields[[extIdx]];
            (* skip if un-colored *)
            If[!TreeMasses`ColorChargedQ[extField], Continue[]];
            (* loop over vertices *)
            For[vertIdx=1, vertIdx<=Length[vertices], vertIdx++,
                vertex = vertices[[vertIdx]];
                fieldsInVertex = vertex[[1]];
                (* check if enternal field is connected to the vertex *)
                If[graph[[extIdx, vertIdx+Length[extFields]]] === 0, Continue[]];
                (* loop over particles in the vertex *)
                For[vertFieldIdx=1, vertFieldIdx<=Length[fieldsInVertex], vertFieldIdx++,
                    inField = fieldsInVertex[[vertFieldIdx]];
                    If[!TreeMasses`ColorChargedQ[inField], Continue[]];
                    If[AntiField[extField] =!= (inField /. f_[_List] -> f), Continue[]];
                    If[MemberQ[Values[keys], GetFieldColorIndex[inField]], Continue[]];
                    vertices = MapAt[
                        (# //. GetFieldColorIndex[inField] :> keys[extIdx])&,
                        vertices,
                        vertIdx
                    ];
                ];
            ];
        ];

        symbol = {};
        (* loop over vertices pairs *)
        vertex1 = vertices[[1]];
        fieldsInVertex1 = vertex1[[1]];
        field1ColorIndexOld = GetFieldColorIndex[fieldsInVertex1[[2]]];
        field1ColorIndexNew = Unique["c"];
        vertices = MapAt[
        (# //. field1ColorIndexOld -> field1ColorIndexNew)&,
        vertices, 1
        ];
        vertices = MapAt[
          (# //. GetFieldColorIndex[vertices[[2,1,2]]] -> field1ColorIndexNew)&,
          vertices, 2
        ];
        field1ColorIndexOld = GetFieldColorIndex[fieldsInVertex1[[3]]];
        field1ColorIndexNew = Unique["c"];
        vertices = MapAt[
          (# //. field1ColorIndexOld -> field1ColorIndexNew)&,
          vertices, 1
        ];
        vertices = MapAt[
          (# //. GetFieldColorIndex[vertices[[3,1,2]]] -> field1ColorIndexNew)&,
          vertices, 3
        ];
        vertices = MapAt[
          (# //. GetFieldColorIndex[vertices[[2,1,3]]] -> GetFieldColorIndex[vertices[[3,1,3]]])&,
          vertices, 3
        ];

        (*For[vertIdx1 = 1, vertIdx1<=Length[vertices], vertIdx1++,*)
            (*vertex1 = vertices[[vertIdx1]];*)
            (*fieldsInVertex1 = vertex1[[1]];*)
            (** loop over fields in vertex1 *)
            (*For[v1i=1, v1i<=Length[fieldsInVertex1], v1i++,*)
                (*field1 = fieldsInVertex1[[v1i]];*)
                (*If[!TreeMasses`ColorChargedQ[field1], Continue[]];*)
                (*field1ColorIndexOld = GetFieldColorIndex[field1];*)
                (** cycle if external particle *)
                (*If[MemberQ[Values[keys], field1ColorIndexOld], Continue[]];*)
                (*field1ColorIndexNew = Unique["c"];*)
                (*vertices = MapAt[*)
                    (*(# //. field1ColorIndexOld -> field1ColorIndexNew)&,*)
                    (*vertices, vertIdx1*)
                (*];*)
                (*AppendTo[symbol, field1ColorIndexNew];*)
                (*For[vertIdx2=vertIdx1+1, vertIdx2<=Length[vertices], vertIdx2++,*)
                    (** cycle if two vertices are not connected *)
                    (*If[graph[[vertIdx1+Length[extFields], vertIdx2+Length[extFields]]] === 0, Continue[]];*)
                    (*vertex2 = vertices[[vertIdx2]];*)
                    (*fieldsInVertex2 = vertex2[[1]];*)
                    (*For[v2i=1, v2i<=Length[fieldsInVertex2], v2i++,*)
                        (*field2 = fieldsInVertex2[[v2i]];*)
                        (*If[!TreeMasses`ColorChargedQ[field2], Continue[]];*)
                        (*field2ColorIndex = GetFieldColorIndex[field2];*)
                        (*If[MemberQ[Values[keys], field2ColorIndex], Continue[]];*)
                        (*If[MemberQ[symbol, field2ColorIndex], Continue[]];*)
                        (** we want to make a propagator < field bar[field]> *)
                        (*If[(field1 /. f_[_List] -> f) =!= (AntiField[field2] /. f_[_List] -> f), Continue[]];*)
                        (*If[vertIdx1 === 4 && v1i === 2 && vertIdx2 =!= 5, Continue[]];*)
                        (*vertices = MapAt[*)
                            (*(# //. field2ColorIndex -> field1ColorIndexNew)&,*)
                            (*vertices, vertIdx2*)
                        (*];*)
                    (*];*)
                (*];*)
            (*];*)
        (*];*)
        vertices
   ];
*)
(* input
    {Fe,bar[Fe],VP,{bar[Fe],Ah,Fe},{Fe,Ah,bar[Fe]},{VP,bar[Fe],Fe}}
    *)
(*
GenerateUniqueColorAssociationsForExternalParticles::notes=
  "Generates unique color indices for external particles"    
GenerateUniqueColorAssociationsForExternalParticles[vvvv_List]:=
  Module[{inOutParticles,inOutColoredParticles,a},
    inOutParticles=TakeWhile[vvvv,(Head[#]=!=List)&];
    (* generate a unique color index for every external particle *)
    a = Association[{}];
    inOutParticlesWithColorIndices = 
    MapIndexed[
      If[TreeMasses`ColorChargedQ[#1],AssociateTo[a,#2[[1]]->Unique["c"]] ]&,
inOutParticles
];
a
    ]
*)

End[];

EndPackage[];
