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

ParticleDefinitions[GaugeES] = {
      {SdL,  {  Description -> "Left Down-Squarks"}},
      {SdR,  { Description -> "Right Down-Squarks"}},
      {SuL,  { Description -> "Left Up-Squarks"}}, 
      {SuR,  { Description -> "Right Up-Squarks" }},                
      {SeL,  { Description -> "Left Selectron"}}, 
      {SeR,  { Description -> "Right Selectron"}},                         
      {SvL,  { Description -> "Left Sneutrino"}},                         
      {SHd0, { Description -> "Neutral Down-Higgs"}},                         
      {SHdm, { Description -> "Charged Down-Higgs"}}, 
      {SHu0, { Description -> "Neutral Up-Higgs"}},                         
      {SHup, { Description -> "Charged Up-Higgs"}},
      {SsR,  { Description -> "Singlet"}},            
      {SDxL,        {Description -> "Left SExotics",
                     FeynArtsNr -> 666,
	             LaTeX -> "\\tilde{Dx}_L",
                     OutputName -> "SDxL"}},
      {SDxbarR,     {Description -> "Right SExotics",
	             FeynArtsNr -> 667,
	             LaTeX -> "\\tilde{Dx}_R",
                     OutputName -> "SDxbR"}},
      {SH1I0,       {Description -> "Neutral Inert-Down-Higgs",
  		     FeynArtsNr -> 101,
  		     LaTeX -> "h^{0Inert}_{1}",
                     OutputName -> "SH1I0" }},
      {SH1Im,  {     Description -> "Charged Inert-Down-Higgs",
                     FeynArtsNr -> 102,
                     LaTeX -> "h^{-Inert}_{1}",
                     OutputName -> "SH1Im"}},
      {SH2I0,  {     Description -> "Neutral Inert-Up-Higgs",
  		     FeynArtsNr -> 103,
  		     LaTeX -> "h^{0Inert}_{2}",
                     OutputName -> "SH2I0"}},
      {SH2Ip,  { Description -> "Charged Inert-Up-Higgs",
  		     FeynArtsNr -> 104,
  		     LaTeX -> "h^{+Inert}_{2}",
                     OutputName -> "SH2Ip"}},
      {SsIR,   { Description -> "Inert-Singlet",
  		     FeynArtsNr -> 105,
  		     LaTeX -> "s^{Inert}",
                     OutputName -> "SSI"}},
       {SHpd0,   { Description -> "Neutral Prime-Higgs",
      		     FeynArtsNr -> 900,
      		     LaTeX -> "H^{'0}",
                     OutputName -> "SHpd0"}},
      {SHpdm,   {    Description -> "Negative charged Prime-Higgs",
      		     FeynArtsNr -> 901,
      		     LaTeX -> "H^{'-}",
                     OutputName -> "SHpdm"}},
      {SHpup,{ Description -> "Positive charged BarPrime-Higgs",
      		     FeynArtsNr -> 902,
      		     LaTeX -> "\\bar{H}^{'+}",
                     OutputName -> "SHpup"}},
      {SHpu0,{ Description -> "Neutral BarPrime-Higgs",
      		     FeynArtsNr -> 903,
      		     LaTeX -> "\\bar{H}^{'0}",
                     OutputName -> "SHpu0"}},    
      {VB,   { Description -> "B-Boson"}},             
      {VG,   { Description -> "Gluon"}},          
      {VWB,  { Description -> "W-Bosons"}},          
      {gB,   { Description -> "B-Boson Ghost"}},   
      {gG,   { Description -> "Gluon Ghost" }},          
      {gWB,  { Description -> "W-Boson Ghost"}},    
      {Glu,  { Description -> "Gluino"}},
      {Wino, { Description -> "Wino"}},
      {Bino, { Description -> "Bino"}},                          
      {H01,        { Description -> "Neutral down Higgsinos",
		     FeynArtsNr -> 122,
		     LaTeX -> "\\tilde{H}^0_1",
		     OutputName -> "H01"}},
      {HC1,         { Description -> "Charged down Higgsinos",
		      FeynArtsNr -> 123,
		      LaTeX -> "\\tilde{H}^-_1",
		      OutputName -> "HC1"}},
      {H02,         { Description -> "Neutral up Higgsinos",
		      FeynArtsNr -> 124,
		      LaTeX -> "\\tilde{H}^0_2",
		      OutputName -> "H02"}},
      {HC2,         { Description -> "Charged up Higgsinos",
		      FeynArtsNr -> 125,
		      LaTeX -> "\\tilde{H}^+_2",
		      OutputName -> "HC2"}},
      (* This does not work in the E6SSM because these can only be combined 
       after the U(1)N symmetry has been broken *) 
      (* {H0,   { Description -> "Neutral Higgsinos"}},
   {HC,   { Description -> "Charged Higgsinos"}},*)           
      (*  {FS,   { Description -> "Singlino" }}, *)
      {H0I1,    { Description -> "Neutral Inert-down-Higgsinos",
                     Width -> 0,
                     FeynArtsNr -> 721,
                     LaTeX -> "\\tilde{h}^{0,Inert}_1",
		    OutputName -> "HNI1"}},
      {HCI1,    { Description -> "Charged Inert-down-Higgsinos",
                      Width -> 0,
                     FeynArtsNr -> 722,
                     LaTeX -> "\\tilde{h}^{-,Inert}_1",
		    OutputName -> "HCI1"}},
      {H0I2,    { Description -> "Neutral  Inert-up-Higgsinos",
                     Width -> 0,
                     FeynArtsNr -> 723,
                     LaTeX -> "\\tilde{h}^{0,Inert}_2",
		    OutputName -> "HNI2"}},
      {HCI2,    { Description -> "Charged Inert-up-Higgsinos",
                     Width -> 0,
                     FeynArtsNr -> 724,
                     LaTeX -> "\\tilde{h}^{+}_2",
		    OutputName -> "HCI2"}}, 
      {Hp0,    { Description -> "Neutral Prime-Higgsinos",
                     Width -> 0,
                     FeynArtsNr -> 950,
                     LaTeX -> "\\tilde{h}^{'0}",
		    OutputName -> "HpN"}},
      {HpC,    { Description -> "Charged Prime-Higgsinos",
                 Width -> 0,
                 FeynArtsNr -> 951,
                 LaTeX -> "\\tilde{h}^{'-}",
                 OutputName -> "HpC"}},
      {Hp01,    { Description -> "Neutral Prime-Higgsinos",
                     Width -> 0,
                     FeynArtsNr -> 950,
                     LaTeX -> "\\tilde{h}^{'0}_1",
		    OutputName -> "HNP1"}},
      {HpC1,    { Description -> "Charged Prime-Higgsinos",
                     Width -> 0,
                     FeynArtsNr -> 951,
                     LaTeX -> "\\tilde{h}^{'-}_1",
		    OutputName -> "HCP1"}},
      {Hp02,    { Description -> "Neutral  Prime-Bar-Higgsinos",
                     Width -> 0,
                     FeynArtsNr -> 952,
                     LaTeX -> "\\tilde{\\bar{h}}^{'0}_2",
		    OutputName -> "HNP2"}},
      {HpC2,    { Description -> "Charged Prime-Bar-Higgsinos",
                     Width -> 0,
                     FeynArtsNr -> 953,
                     LaTeX -> "\\tilde{\\bar{h}}^{'+}_2",
		    OutputName -> "HCP2"}},
      (*  {FDX,        { Description -> "Dirac Exotics",
        	     LaTeX -> "x_1",
		     FeynArtsNr -> 660,
                         OutputName -> "FDx"}}, *)
      {FDx1,        { Description -> "Dirac Left Exotics",
        	     LaTeX -> "Dx_1",
		     FeynArtsNr -> 660,
		     OutputName -> "Dx1"}},
      {FDx2,        { Description -> "Dirac Right Exotics",
        	     LaTeX -> "Dx_2",
		     FeynArtsNr -> 661,
		     OutputName -> "Dx2"}},
      {FS1,        { Description -> "Dirac Left Singlino",
		     FeynArtsNr -> 806,
  		     LaTeX -> "\\tilde{s}_L",
		     OutputName -> "FSL"}},
      {FS2,        { Description -> "Dirac Right Singlino",
		     FeynArtsNr -> 807,
  		     LaTeX -> "\\tilde{s}_R",
		     OutputName -> "FSR"}},
      {FSI1,    { Description -> "Inert-Singlino-1",
                      FeynArtsNr -> 905,
                      LaTeX -> "\\tilde{S}^{Inert}_1",
                      OutputName -> "FSI1"}},
      {FSI2,    { Description -> "Inert-Singlino-2",
                      FeynArtsNr -> 905,
                      LaTeX -> "\\tilde{S}^{Inert}_2",
                      OutputName -> "FSI2"}},
      {Fd1,  { Description -> "Dirac Left Down-Quark"}},
      {Fd2,  { Description -> "Dirac Right Down-Quark"}},
      {Fu1,  { Description -> "Dirac Left Up-Quark"}},
      {Fu2,  { Description -> "Dirac Right Up-Quark"}},
      {Fe1,  { Description -> "Dirac Left Electron"}},
      {Fe2,  { Description -> "Dirac Right Electron"}},
      {Fv,   { Description -> "Neutrinos" }},
      {VBp,   { LaTeX -> "Bp",
               OutputName -> "VBp",
               FeynArtsNr -> 10 }},
      {gBp,   { LaTeX -> "g_Bp",
               OutputName -> "gBp",
               FeynArtsNr -> 10 }},
      {FBp,  { LaTeX -> "\\tilde{Bp}_1",
               FeynArtsNr -> 10,
               OutputName -> "FBp"}}
      };
      
      
  ParticleDefinitions[TEMP] = {
    {VZ,   { Description -> "Z-Boson" }}, 
    {gZ,   { Description -> "Z-Boson Ghost" }}
    };
      
      
  ParticleDefinitions[EWSB] = {

     {Sd ,  { Description -> "Down-Squarks"}},
     {Su ,  { Description -> "Up-Squarks"}},   
     {Se ,  { Description -> "Sleptons"}}, 
     {Sv ,  { Description -> "Sneutrinos"}}, 
     {SDX,   { Description -> "SExotics",
              PDG -> {1000051,2000051,1000052,2000052,1000053,2000053},
              Mass -> Automatic,
              ElectricCharge -> -1/3,
              FeynArtsNr -> 666,
              LaTeX -> "\\tilde{x}",
              OutputName -> "SX"}},
     {ChiI, { Description -> "Inert Neutralinos",
              PDG -> {1000081,1000082,1000083,1000084},
              Mass -> Automatic,
              FeynArtsNr -> 13,
              ElectricCharge -> 0,
              LaTeX -> "\\tilde{\\chi}^{0,Inert}",
              OutputName -> "NI"}},
     {ChaI, { Description -> "Inert Charginos",
              PDG -> {1000085,1000086},
              Mass -> Automatic,
              ElectricCharge -> -1,
              FeynArtsNr -> 14,
              LaTeX -> {"\\tilde{\\chi}^{-,Inert}",
                        "\\tilde{\\chi}^{+,Inert}"},
              OutputName -> "AI"}},
     {FSI,  { Description -> "Inert-Singlino",
              PDG -> {1000089,1000090},
              FeynArtsNr -> 105,
              ElectricCharge -> 0,
              LaTeX -> "\\tilde{S}^{Inert}",
              OutputName -> "FSI"}},
     {hh   ,  {  Description -> "Higgs", 
                 PDG -> {25, 35,45},
                 PDG.IX ->{101000001,101000002,101000003} }}, 
     {Ah   ,  {    Description -> "Pseudo-Scalar Higgs",
                 PDG -> {0, 0, 36},
                 PDG.IX ->{0,0, 102000001} }},                

  (* {hh ,  { Description -> "Higgs", 
              PDG -> {25, 35,45},
              PDG.IX ->{101000001,101000002,101000003} }}, 
     {Ah ,  { Description -> "Pseudo-Scalar Higgs", 
              LaTeX -> "A",
              PDG -> {0, 0, 36},
              PDG.IX ->{0,0, 102000001},
              FeynArtsNr -> 2 }}, *)                      
      {Hpm,  { Description -> "Charged Higgs"}},   
      {SHI0, { Description -> "Neutral Inert-Higgs",
               PDG -> {82,86,84,88},
               FeynArtsNr -> 7,
               ElectricCharge -> 0,
               LaTeX -> "h^{0,Inert}",
               OutputName -> "HNI"}},
      {SHIp, { Description -> "Charged Inert-Higgs",
               PDG -> {81,85,83,87},
               FeynArtsNr -> 8,
               ElectricCharge -> -1 ,
               LaTeX -> {"h^{-,Inert}","h^{+,Inert}"},
               OutputName -> "HCI"}}, 
     {SSI0,   { Description -> "Inert-Singlet",
                    PDG -> {89,90},
                    FeynArtsNr -> 105,
                    ElectricCharge -> 0,
                    LaTeX -> "s^{Inert}",
                    OutputName -> "SSI"}},              
     (* {SsIR,  { Description -> "Inert-Singlet",
               PDG -> {89,90},
               FeynArtsNr -> 105,
               ElectricCharge -> 0,
               LaTeX -> "s^{Inert}",
               OutputName -> "SSI"}},*)
      {SHp0,   { Description -> "Neutral Prime-Higgs",
                     PDG -> {92,94},
      		     FeynArtsNr -> 900,
		     ElectricCharge -> 0,
      		     LaTeX -> "H^{'0}",
		     OutputName -> "HPN"}},
      {SHpp,   { Description -> "Charged Prime-Higgs",
                     PDG -> {91,93},
		     ElectricCharge -> -1,
      		     FeynArtsNr -> 901,
      		     LaTeX -> {"H^{'-}","H^{'+}"},
		     OutputName -> "HPA"}},
      {VP,   { Description -> "Photon"}},
      {VZp,    { Description -> "Z'-Boson"}},  
      {gZp,    { Description -> "Z'-Ghost" }},  
      {VZ,   { Description -> "Z-Boson" }}, 
      {gZ,   { Description -> "Z-Boson Ghost" }},  
      {VG,   { Description -> "Gluon" }},          
      {VWm,  { Description -> "W-Boson" }},         
      {gP,   { Description -> "Photon Ghost"}}, 
      {gWm,  { Description -> "Negative W-Boson Ghost"}}, 
      {gWmC, { Description -> "Positive W-Boson Ghost" }}, 
      {gG,   { Description -> "Gluon Ghost" }},                 
      {Fd,   { Description -> "Down-Quarks"}},   
      {Fu,   { Description -> "Up-Quarks"}},   
      {Fe,   { Description -> "Leptons" }},
      {Fv,   { Description -> "Neutrinos" }}, 
      {FDX,  { Description->"Exotics",
    		   PDG -> {51,52,53},
    		   Width -> 0,
    		   Mass -> Automatic,
		   ElectricCharge -> -1/3,
    		   FeynArtsNr -> 666,
    		   LaTeX -> "x",
    		   OutputName -> "FDX"}},
      {Glu,  { Description -> "Gluino" }},
      {Chi,  { Description -> "Neutralinos",
               PDG -> {1000022,1000023,1000025,1000035,1000045,1000055},
               PDG.IX ->{211000001,211000002,211000003,211000004,211000005,211000006},
               FeynArtsNr -> 11 }},
     {Cha,  { Description -> "Charginos",
              FeynArtsNr -> 12}},
      {ChiP,   { Description -> "Prime Neutralinos",
                     PDG -> {1000092,1000094},
		     ElectricCharge -> 0,
      		     FeynArtsNr -> 900,
      		     LaTeX -> "\\tilde{\\chi}^{'0}",
		     OutputName -> "NP"}},
      {ChaP,   { Description -> "Prime Chargino",
                     PDG -> {1000091},
		     ElectricCharge -> -1,
    		     (*Mass -> Automatic,*)
      		     FeynArtsNr -> 901,
      		     LaTeX -> {"\\tilde{\\chi}^{'-}",
      			       "\\tilde{\\chi}^{'+}"},
		     OutputName -> "AP"}}                            
        };    
        
        
        
 WeylFermionAndIndermediate = {
       {FHd0,  { Description -> "Neutral Down-Higgsino"}},      
       {FHu0,  { Description -> "Neutral Up-Higgsino" }}, 
       {FHdm,  { Description -> "Charged Down-Higgsino"}}, 
       {FHup,  { Description -> "Charged Up-Higgsino"}},
       {FHpd0, { Description -> "Neutral Higgsino-Prime",
                     FeynArtsNr -> 1053,
                     LaTeX -> "\\tilde{h'}^{0}"}},
       {FHpu0, { Description -> "Neutral Higgsino-Bar-Prime",
                 FeynArtsNr -> 1063,
                 LaTeX -> "\\tilde{\\bar{h'}}^{0}"}},
       {FHpdm, { Description -> "Charged Higgsino-Prime",
                 FeynArtsNr -> 1054,
                 LaTeX -> "\\tilde{h'}^{0}"}},
       {FHpup, { Description -> "Charged Higgsino-Bar-Prime",
                 FeynArtsNr -> 1064,
                 LaTeX -> "\\tilde{\\bar{h'}^}{+}"}},
       {FH1I0, { Description -> "Neutral Inert-down-Higgsino",
                     FeynArtsNr -> 821,
                     LaTeX -> "\\tilde{h}^{0,Inert}_{d}"}},
       {FH2I0, { Description -> "Neutral Inert-up-Higgsino",
                     FeynArtsNr -> 823,
                     LaTeX -> "\\tilde{h}^{0,Inert}_{u}"}},
       {FH1Im, { Description -> "Charged Inert-down-Higgsino",
                     FeynArtsNr -> 822,
                     LaTeX -> "\\tilde{h}^{-,Inert}_{d}"}},
       {FH2Ip, { Description -> "Charged Inert-up-Higgsino",
                     FeynArtsNr -> 824,
                     LaTeX -> "\\tilde{h}^{+,Inert}_{u}"}},
       {FsIR,    { Description -> "Weyl-Inert-Singlino",
                       PDG -> {1000089,1000090},
                       FeynArtsNr -> 205,
                       LaTeX -> "\\tilde{S}^{Inert}",
                       OutputName -> "FS0I"}},
       {L0,    { Description -> "Neutralino Weyl-Spinor"}},
       {Lm,    { Description -> "Negative Chargino Weyl-Spinor"}},
       {Lp,    { Description -> "Positive Chargino Weyl-Spinor"}}, 
       {L0I,   { Description -> "Neutralino Inert-Weyl-Spinor"}},
       {LmI,    { Description -> "Negative Chargino Inert-Weyl-Spinor"}},
       {LpI,    { Description -> "Positive Chargino Inert-Weyl-Spinor"}},
       {L0p,    { Description -> "Neutralino Prime-Weyl-Spinor"}},
       {LS0I,    { Description -> "Neutralino inert-singlino-Weyl-Spinor"}},
       
       {fG,    { Description ->"Gluino Weyl-Spinor"}},
       {fWB,   { Description ->"Wino Weyl-Spinor"}},
       {fW0,   { Description ->"Neutral Wino" }},
       {fWm,   { Description ->"Negative Wino"}},                 
       {fWp,   { Description ->"Positive Wino"}},                 
       {fB,    { Description ->"Bino Weyl-Spinor"}},    
       {phid,  { Description -> "Scalar Down"}},                                                                       
       {phiu,  { Description -> "Scalar Up"}}, 
       {sigmad,   { Description -> "Pseudo Scalar Down"}},                                                                                        
       {sigmau,   { Description -> "Pseudo Scalar Up" }},
       {SHd,  { Description -> "Down-Higgs"}},
       {SHu,  { Description -> "Up-Higgs"}},
       {Sl,   { Description -> "Left Slepton" }},
       {Sq,   { Description -> "Left Squark" }},
       {FeL,  { Description -> "Left Electron" }},
       {FeR,  { Description -> "Right Electron" }} ,
       {FdL,  { Description -> "Left Down-Quark" }},
       {FdR,  { Description -> "Right Down-Quark" }},              
       {FuL,  { Description -> "Left Up-Quark" }},
       {FuR,  { Description -> "Right Up-Quark" }},
       {FEL,  { Description -> "Rotated Left Electron" }},
       {FER,  { Description -> "Rotated Right Electron" }} ,
       {FDL,  { Description -> "Rotated Left Up-Quark" }},
       {FDR,  { Description -> "Rotated Right Up-Quark" }},              
       {FUL,  { Description -> "Rotated Left Down-Quark"}},
       {FUR,  { Description -> "Rotated Right Down-Quark" }},                
       {FHd,  { Description -> "Down-Higgsino" }},
       {FHu,  { Description -> "Up-Higgsino" }},
       {Fl,   { Description -> "Left Leptons"}},
       {Fq,   { Description -> "Left Quarks"}}, 
       {FvL,  { Description -> "Left Neutrino"}}, 
       {FsR,   { Description -> "Weyl Spinor of Singlino"}}, 
       
       {fBp,   { LaTeX -> "\\lambda_Bp" }},
       {FDXL,   { Description -> "Rotated Left Exotics",
                 LaTeX -> "X_L"}},
       {FDXR,   { Description->"Rotated Right Exotics",
                 LaTeX -> "X_L"}},         
 {sigmaS,      {Description -> "Scalar Singlet" }}  ,
                 
 {phiS,      { Description -> "Pseudo Scalar Singlet"}},


       {e,    { Description -> "Right Electron Superfield" }},
       {d,    { Description -> "Right Down-Quark Superfield" }},                 
       {q,    { Description -> "Left Quark Superfield" }},                 
       {u,    { Description -> "Right Up-Quark Superfield" }},                 
       {l,    { Description -> "left Lepton Superfield" }},  
       {Hd,   { Description -> "Down-Higgs Superfield"}},                                         
       {Hu,   { Description -> "Up-Higgs Superfield" }},                 
       {s,    { Description -> "Singlet Superfield" }}, 
       {G,    { Description -> "Gluon Superfield" }},                 
       {B,    { Description -> "B Superfield" }},                 
       {WB,   { Description -> "W Superfield" }},
       {Bp,    { Description -> "U(1)' gauge Superfield", 
                 LaTeX -> "\\hat{Bp}"}},
      {H1I,    { Description -> "Inert-Down-Higgs Superfield",
		     LaTeX -> "\\hat{H}^{Inert}_{1}"}},
      {H1,    { Description -> "Inert-Up-Higgs Superfield",
		     LaTeX -> "\\hat{H}^{Inert}_{2}"}},
      {sI,     { Description -> "Inert-Singlet Superfield",
		     LaTeX -> "\\hat{S}^{Inert}"}},    
       {MassZp,    { LaTeX -> "M_{Zp}",
                    Form -> Scalar }},
      {SH1I,   { Description-> "Scalar Inert Down Higgs",
		     LaTeX -> "H^{Inert}_1"}},
      {SH2I,   { Description-> "Scalar Inert Up Higgs",
		     LaTeX -> "H^{Inert}_2"}},
      {SHp,   { Description-> "Scalar Prime Higgs",
		     LaTeX -> "H'"}},
      {SHpbar,   { Description-> "Scalar Bar Prime Higgs",
		     LaTeX -> "\\bar{H'}"}},
      {Dx,          { Description -> "Left Exotics Superfield",
		     LaTeX -> "\\hat{Dx}"}},
      {Dxbar,       { Description -> "Right Exotics Superfield",
		     LaTeX -> "\\hat{\\bar{Dx}}"}},
        {Hp,          { Description -> "HPrime Superfield",
                     LaTeX -> "\\hat{H'}"}},
      {Hpbar,       { Description -> "HPrimeBar Superfield",
                     LaTeX -> "\\hat{\\bar{H'}}"}}
};       


