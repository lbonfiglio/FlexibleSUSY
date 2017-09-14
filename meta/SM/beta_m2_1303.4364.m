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

(* Beta function for mass parameter*) 
   
bms = h*((-9*al1)/20 - (9*al2)/4 + atau + 12*lam + ab*NR + at*NR) + 
     h^2*((-9*atau^2)/4 - 120*lam^2 - (9*ab^2*NR)/4 - (7*ab*at*NR)/2 - 
       (9*at^2*NR)/4 + lam*((72*al1)/5 + 72*al2 - 24*atau - 24*ab*NR - 
         24*at*NR) + (3*al1*((15*al2)/16 + (25*atau)/8 + (25*ab*NR)/72 + 
          (85*at*NR)/72))/5 + al2*((15*atau)/8 + (15*ab*NR)/8 + 
         (15*at*NR)/8) + al3*(5*ab*cR*NR + 5*at*cR*NR) + 
       (9*al1^2*(157/96 + (5*NG)/8 + (55*NG*NR)/216))/25 + 
       al2^2*(-385/32 + (5*NG)/8 + (5*NG*NR)/8)) + 
     h^3*(8208*lam^3 + 24*ab^2*atau*NR + 24*ab*atau^2*NR + 
       lam^2*(198*atau + 198*ab*NR + 198*at*NR + al2*(-252 - 432*z3) + 
         (3*al1*(-84 - 144*z3))/5) + atau^3*(-233/16 + 15*z3) + 
       ab^3*((-617*NR)/16 + 24*NR^2 + 15*NR*z3) + 
       at^3*((-617*NR)/16 + 24*NR^2 + 15*NR*z3) + 
       (27*al1^3*(839/108 + (1145*NG)/144 + (35*NG^2)/48 + 
          (7195*NG*NR)/3888 + (385*NG^2*NR)/648 + (4235*NG^2*NR^2)/34992 - 
          (51*z3)/16 - (27*NG*z3)/4 - (137*NG*NR*z3)/108))/125 + 
       al2^3*(-39415/576 + (2867*NG)/288 + (35*NG^2)/144 + (2867*NG*NR)/288 + 
         (35*NG^2*NR)/72 + (35*NG^2*NR^2)/144 + (711*z3)/16 + (45*NG*z3)/4 + 
         (45*NG*NR*z3)/4) + at^2*(24*atau*NR + ab*((29*NR)/16 + (55*NR^2)/2 + 
           12*NR*z3)) + at*((7*ab*atau*NR)/2 + 24*atau^2*NR + 
         ab^2*((29*NR)/16 + (55*NR^2)/2 + 12*NR*z3)) + 
       al3*(ab^2*((447*cR*NR)/8 - 90*cR*NR*z3) + 
         at^2*((447*cR*NR)/8 - 90*cR*NR*z3) + ab*at*((41*cR*NR)/4 - 
           12*cR*NR*z3)) + al3^2*(ab*((77*cA*cR*NR)/2 - (119*cR^2*NR)/4 - 
           16*cR*NG*NR*TF - 18*cA*cR*NR*z3 + 36*cR^2*NR*z3) + 
         at*((77*cA*cR*NR)/2 - (119*cR^2*NR)/4 - 16*cR*NG*NR*TF - 
           18*cA*cR*NR*z3 + 36*cR^2*NR*z3)) + 
       al2^2*(atau*(-255/128 - (21*NG)/32 - (21*NG*NR)/32 - (81*z3)/4) + 
         ab*((-255*NR)/128 - (21*NG*NR)/32 - (21*NG*NR^2)/32 - 
           (81*NR*z3)/4) + at*((-255*NR)/128 - (21*NG*NR)/32 - 
           (21*NG*NR^2)/32 - (81*NR*z3)/4) + al3*((135*cR*NG*NR)/16 - 
           9*cR*NG*NR*z3)) + 
       (9*al1^2*(atau*(-3607/128 - (117*NG)/32 - (143*NG*NR)/96 - 
            (15*z3)/4) + at*((-123103*NR)/10368 - (127*NG*NR)/96 - 
            (1397*NG*NR^2)/2592 - (149*NR*z3)/36) + 
          ab*((-79207*NR)/10368 - (31*NG*NR)/96 - (341*NG*NR^2)/2592 - 
            (35*NR*z3)/36) + al2*(1053/32 + (63*NG)/16 + (47*NG*NR)/48 - 
            (207*z3)/16 - (9*NG*z3)/4 - (NG*NR*z3)/4) + 
          al3*((55*cR*NG*NR)/16 - (11*cR*NG*NR*z3)/3)))/25 + 
       (3*al1*((-3*ab*atau*NR)/2 + atau^2*(291/16 - 36*z3) + 
          at^2*((-323*NR)/48 - (3*NR^2)/4 + 4*NR*z3) + 
          ab^2*((-959*NR)/48 - (3*NR^2)/4 + 24*NR*z3) + 
          al2^2*(2691/64 + (87*NG)/32 + (27*NG*NR)/32 - (405*z3)/16 - 
            (9*NG*z3)/4 - (NG*NR*z3)/4) + at*((-3*atau*NR)/2 + 
            ab*((-605*NR)/72 - (3*NR^2)/2 - (NR*z3)/3)) + 
          al2*((-865*ab*NR)/192 + atau*(-2331/64 + 45*z3) + 
            at*((-3277*NR)/192 + (39*NR*z3)/2)) + 
          al3*(ab*((-991*cR*NR)/144 + 5*cR*NR*z3) + 
            at*((-2419*cR*NR)/144 + 17*cR*NR*z3))))/5 + 
       al2*((-9*ab*atau*NR)/2 + atau^2*(-987/16 + 54*z3) + 
         ab^2*((-951*NR)/16 - (9*NR^2)/4 + 54*NR*z3) + 
         at^2*((-951*NR)/16 - (9*NR^2)/4 + 54*NR*z3) + 
         at*((-9*atau*NR)/2 + ab*((27*NR)/8 - (9*NR^2)/2 - 9*NR*z3)) + 
         al3*(ab*((-489*cR*NR)/16 + 27*cR*NR*z3) + 
           at*((-489*cR*NR)/16 + 27*cR*NR*z3))) + 
       lam*(-72*ab*atau*NR + al2^2*(11511/16 - (153*NG)/4 - (153*NG*NR)/4 - 
           324*z3) + (9*al1^2*(-1077/16 - (153*NG)/4 - (187*NG*NR)/12 - 
            36*z3))/25 + atau^2*(261/2 + 144*z3) + 
         ab^2*((333*NR)/2 - 36*NR^2 + 144*NR*z3) + 
         at^2*((333*NR)/2 - 36*NR^2 + 144*NR*z3) + 
         al2*(atau*(189/4 - 216*z3) + ab*((189*NR)/4 - 216*NR*z3) + 
           at*((189*NR)/4 - 216*NR*z3)) + at*(-72*atau*NR + 
           ab*(111*NR - 72*NR^2 - 144*NR*z3)) + 
         (3*al1*(al2*(-1701/8 + 72*z3) + atau*(-549/4 + 72*z3) + 
            ab*((131*NR)/4 - 88*NR*z3) + at*((-73*NR)/4 - 40*NR*z3)))/5 + 
         al3*(ab*(-306*cR*NR + 288*cR*NR*z3) + at*(-306*cR*NR + 
             288*cR*NR*z3))));
