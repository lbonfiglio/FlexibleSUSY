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

AppendTo[$Path, FileNameJoin[{Directory[], "meta"}]];

Needs["SARAH`"];
Needs["FlexibleSUSY`", "FlexibleSUSY.m"];
Needs["TestSuite`", "TestSuite.m"];
Needs["TreeMasses`", "TreeMasses.m"];

SARAH`SARAH[OutputDirectory] = FileNameJoin[{Directory[], "Output"}];
SARAH`SARAH[InputDirectories] = {
    FileNameJoin[{Directory[], "sarah"}],
    ToFileName[{$sarahDir, "Models"}]
};

Start["MRSSM"];

Print["testing IsMassless[] ..."];

TestEquality[TreeMasses`IsMassless[gG], True];
TestEquality[TreeMasses`IsMassless[VG], True];
TestEquality[TreeMasses`IsMassless[gP], True];
TestEquality[TreeMasses`IsMassless[VP], True];
TestEquality[TreeMasses`IsMassless[gZ], False];
TestEquality[TreeMasses`IsMassless[VZ], False];
TestEquality[TreeMasses`IsMassless[gWp], False];
TestEquality[TreeMasses`IsMassless[VWp], False];

Print["testing getters for particle collection..."];

TestEquality[TreeMasses`GetSusyParticles[], {Glu, SRdp, SRum, sigmaO, phiO, Sd, Sv, Su, Se, hh,
   Ah, Rh, SARAH`Hpm, Chi, Cha1, Cha2}];
TestEquality[TreeMasses`GetColoredParticles[], {VG, gG, Glu, sigmaO, phiO, Sd, Su, Fd, Fu}];
TestEquality[TreeMasses`GetVectorBosons[], {VG, VP, VZ, VWm}];

Print["testing getters for specific particles..."];

TestEquality[TreeMasses`GetPhoton[], VP];
TestEquality[TreeMasses`GetGluon[], VG];
TestEquality[TreeMasses`GetZBoson[], VZ];
TestEquality[TreeMasses`GetWBoson[], VWm];
TestEquality[TreeMasses`GetHiggsBoson[], hh];
TestEquality[TreeMasses`GetChargedHiggsBoson[], Hpm];
TestEquality[TreeMasses`GetPseudoscalarHiggsBoson[], Ah];


TestEquality[Select[GetParticles[], IsScalar], {SRdp, SRum, sigmaO, phiO, Sd, Sv, Su, Se, hh, Ah, Rh, Hpm}];

PrintTestSummary[];
