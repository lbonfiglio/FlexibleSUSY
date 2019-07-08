#!/usr/bin/env python
# base slha file
FS_slha_setup = """
Block MODSEL                 # Select model
    6   0                    # flavour violation
Block FlexibleSUSY
    0   1.000000000e-05      # precision goal
    1   0                    # max. iterations (0 = automatic)
    2   0                    # algorithm (0 = two_scale, 1 = lattice)
    3   1                    # calculate SM pole masses
    4   1                    # pole mass loop order
    5   2                    # EWSB loop order
    6   3                    # beta-functions loop order
    7   2                    # threshold corrections loop order
    8   1                    # Higgs 2-loop corrections O(alpha_t alpha_s)
    9   1                    # Higgs 2-loop corrections O(alpha_b alpha_s)
   10   1                    # Higgs 2-loop corrections O((alpha_t + alpha_b)^2)
   11   1                    # Higgs 2-loop corrections O(alpha_tau^2)
   12   0                    # force output
   13   1                    # Top quark 2-loop corrections QCD
   14   0                    # Higgs logarithmic resummation
   15   1.000000000e-11      # beta-function zero threshold
"""

SPheno_slha_setup = """
Block MODSEL      #
1 0               #  1/0: High/low scale input
2 1              # Boundary Condition
6 0               # Generation Mixing
12  {Ms: .16E}     #
Block SPhenoInput   # SPheno specific input
  1 -1              # error level
  2  0              # SPA conventions
  7  0              # Skip 2-loop Higgs corrections
  8  3              # Method used for two-loop calculation
  9  1              # Gaugeless limit used at two-loop
 10  0              # safe-mode used at two-loop
 11 0               # calculate branching ratios
 13 1               # include 3-Body decays
 12 1.000E-04       # write only branching ratios larger than this value
 15 1.000E-30       # write only decay if width larger than this value
 31 -1              # fixed GUT scale (-1: dynamical GUT scale)
 32 0               # Strict unification
 34 1.000E-04       # Precision of mass calculation
 35 80              # Maximal number of iterations
 36 40              # Minimal number of iterations before discarding point
 37 1               # Set Yukawa scheme
 38 2               # 1- or 2-Loop RGEs
 50 1               # Majorana phases: use only positive masses
 51 0               # Write Output in CKM basis
 52 1               # Write spectrum in case of tachyonic states
 55 1               # Calculate one loop masses
 57 0               # Calculate low energy constraints
 65 1               # Solution tadpole equation
 75 1               # Write WHIZARD files
 76 1               # Write HiggsBounds file
 86 0.              # Maximal width to be counted as invisible in Higgs decays; -1: only LSP
510 0.              # Write tree level values for tadpole solutions
515 0               # Write parameter values at GUT scale
520 1.              # Write effective Higgs couplings (HiggsBounds blocks)
525 0.              # Write loop contributions to diphoton decay of Higgs
530 1.              # Write Blocks for Vevacious
"""

slha_sminputs="""
BLOCK SMINPUTS               # Standard Model inputs
    1   1.279440000e+02      # alpha^(-1) SM MSbar(MZ)
    2   1.166380000e-05      # G_Fermi
    3   1.184000000e-01      # alpha_s(MZ) SM MSbar
    4   9.118760000e+01      # MZ(pole)
    5   4.180000000e+00      # mb(mb) SM MSbar
    6   1.733400000e+02      # mtop(pole)
    7   1.777000000e+00      # mtau(pole)
    8   0.000000000e+00      # mnu3(pole)
   11   5.109989020e-04      # melectron(pole)
   12   0.000000000e+00      # mnu1(pole)
   13   1.056583570e-01      # mmuon(pole)
   14   0.000000000e+00      # mnu2(pole)
   21   4.750000000e-03      # md(2 GeV) MS-bar
   22   2.400000000e-03      # mu(2 GeV) MS-bar
   23   1.040000000e-01      # ms(2 GeV) MS-bar
   24   1.270000000e+00      # mc(mc) MS-bar"""
FS_model_inputs = """
Block EXTPAR
    61  {lam: .16E}    # lambda (MSUSY)
    62  {kap: .16E}    # kappa (MSUSY)
Block Ms
    {Ms: .16E}             # SUSY scale
Block TanBeta
    {tbMSUSY: .16E}        # tan(Beta) at the SUSY scale
Block Xtt
    {Xt: .16E}             # Xt
"""
spheno_model_inputs = """
BLOCK EXTPAR	 	 #
0	{Ms: .16E}	 # MSUSY
1	{Ms: .16E}	 # M1
2	{Ms: .16E}	 # M2
3	{Ms: .16E}	 # M3
25	{tb: .16E}	 # tb(MZ) - check this because this breaks SLHA
31	{Ms: .16E}   	 #
32	{Ms: .16E}	 #
33	{Ms: .16E}	 #
34	{Ms: .16E}	 #
35	{Ms: .16E}	 #
36	{Ms: .16E}	 #
41	{Ms: .16E}	 #
42	{Ms: .16E}	 #
43	{Ms: .16E}	 #
44	{Ms: .16E}	 #
45	{Ms: .16E}	 #
46	{Ms: .16E}	 #
47	{Ms: .16E}	 #
48	{Ms: .16E}	 #
49	{Ms: .16E}	 #
61	{lam: .16E}	 #
62	{kap: .16E}	 #
63	{Al: .16E}	 #
64	{Ak: .16E}	 #
65	{Ms: .16E}	 #
1011    {Tu11In: .16E},  #
1012    {Td11In: .16E},  #
1013    {Te11In: .16E},  #
1021    {Tu22In: .16E},  #
1022    {Td22In: .16E},  #
1023    {Te22In: .16E},  #
1031    {Tu33In: .16E},  #
1032    {Td33In: .16E},  #
1033    {Te33In: .16E}	 #
"""

nmssmtools_setup = """
BLOCK MODSEL
	3	1		# NMSSM particle content
	1	0		# IMOD (0=general NMSSM, 1=mSUGRA, 2=GMSB)
	10	0		# ISCAN (0=no scan, 1=grid scan, 2=random scan, 3=MCMC)
	9	0		# |OMGFLAG|=0: no (default), =1: relic density only,
#				  =2: dir. det. rate, =3: indir. det. rate, =4: both,
#				  OMGFLAG>0: 0.107<OMG<0.131, <0: OMG<0.131
	8       2               # Precision for Higgs masses (default 0: as before,
#                                 1: full 1 loop + full 2 loop from top/bot Yukawas
#				  2: as 1 + pole masses - 1&2 by courtesy of P. Slavich)
	13      0               # 1: Sparticle decays via NMSDECAY (default 0)
	11      0               # Constraints on (g-2)_muon (1=yes, 0=no, default=1)
	14      0               # 0: H-> VV,VV* (default); 1: H->VV,VV*,V*V*
	15	0		# Precision for micromegas (defalt=0):
#				  +0/1: fast computation on/off
#				  +0/2: Beps=1d-3, 1d-6
#				  +0/4: virtual Ws off/on
#
#
"""
nmssmtools_model_inputs = """
BLOCK MINPAR
	3	{tb: .16E}	# TANB at MZ
#
BLOCK EXTPAR
	0	{Ms: .16E} 	# MSUSY (If =/= SQRT(2*MQ1+MU1+MD1)/2)
        1	{Ms: .16E}	# M1 (If =/= M2/2)
	2	{Ms: .16E}	# M2
        3	{Ms: .16E}	# M3
	11	{At: .16E}	# AU3
	12	{Ab: .16E}	# AD3
	13	{Atau: .16E}	# AE3
        16	{Amu: .16E}	# AE2
        31	{Ms: .16E}   	#
        32	{Ms: .16E}	 #
	33	{Ms: .16E}	 #
        34	{Ms: .16E}	 #
	35	{Ms: .16E}	 #
        36	{Ms: .16E}	 #
	41	{Ms: .16E}	 #
	42	{Ms: .16E}	 #
	43	{Ms: .16E}	 #
        44	{Ms: .16E}	 #
        45	{Ms: .16E}	 #
        46	{Ms: .16E}	 #
	47	{Ms: .16E}	 #
        48	{Ms: .16E}	 #
	49	{Ms: .16E}	 #
        61	{lam: .16E}	 #
	62	{kap: .16E}	 #
	63	{Al: .16E}	 #
	64	{Ak: .16E}	 #
    	65	{Ms: .16E}	 # """
softsusy_setup="""
Block MODSEL                 # Select model
    6    0                   # flavour violation
    1    0                   # mSUGRA
    3    1                   # NMSSM
Block SOFTSUSY               # SOFTSUSY specific inputs
    1   1.000000000e-04      # tolerance
    2   2.000000000e+00      # up-quark mixing (=1) or down (=2)
    5   1.000000000E+00      # 2-loop running
    3   1.000000000E+00      # printout
    7   2
   15   1.000000000E+00      # NMSSMTools compatible output (default: 0)
   16   0.000000000E+00      # Call micrOmegas (default: 0 = no,
                             # 1 = relic density only,
                             # 2 = direct detection + relic density,
                             # 3 = indirect detection + relic density,
                             # 4 = all)
   17   1.000000000E+00      # Sparticle decays via NMSDECAY (default: 0)
   18   1.000000000E+00      # use soft Higgs masses as EWSB output
"""
nmssmcalc_model_inputs = """
BLOCK MINPAR
	3	{tbMSUSY: .16E}	# TANB at MSUSY
#
BLOCK EXTPAR
0	{Ms: .16E} 	# MSUSY (If =/= SQRT(2*MQ1+MU1+MD1)/2)
1	{Ms: .16E}	# M1 (If =/= M2/2)
2	{Ms: .16E}	# M2
3	{Ms: .16E}	# M3
11	{At: .16E}	# AU3
12	{Ab: .16E}	# AD3
13	{Atau: .16E}	# AE3
16	{Amu: .16E}	# AE2
25      {tbMSUSY: .16E}
31	{Ms: .16E}   	#
32	{Ms: .16E}	 #
33	{Ms: .16E}	 #
34	{Ms: .16E}	 #
35	{Ms: .16E}	 #
36	{Ms: .16E}	 #
41	{Ms: .16E}	 #
42	{Ms: .16E}	 #
43	{Ms: .16E}	 #
44	{Ms: .16E}	 #
45	{Ms: .16E}	 #
46	{Ms: .16E}	 #
47	{Ms: .16E}	 #
48	{Ms: .16E}	 #
49	{Ms: .16E}	 #
61	{lam: .16E}	 #
62	{kap: .16E}	 #
63	{Al: .16E}	 #
64	{Ak: .16E}	 #
65	{Ms: .16E}	 #"""

nmssmcalc_setup = """Block MODSEL
  3  1   # NMSSM
  5  0   # 0: CP-conserving; 2: general CP-violation
  6  2   # loop level 1: one 2:two
  7  1   # 1: DRbar scheme for top/stop-sector; 2: OS scheme for top/stop-sector
  8  9   #  """

nmssmcalcOS_setup ="""Block MODSEL
  3  1   # NMSSM
  5  0   # 0: CP-conserving; 2: general CP-violation
  6  2   # loop level 1: one 2:two
  7  1   # 1: DRbar scheme for top/stop-sector; 2: OS scheme for top/stop-sector
  8  9   #"""

def fill_block_3b3(line, myfile, M):
#    OK for trilinears in SOFTSUSY not safe in genera\l
    while(line.split()[0] != "3"  or line.split()[1] != "3"  ) :
        line = next(myfile)
        k = int(line.split()[0])
        j = int(line.split()[1])
        M[k-1][j-1] = line.split()[2]

def fill_light_higgs(slhaout,mh1):
    for line in slhaout:
        if 'Block MASS' in line:
            #skip spheno's comment
            next(slhaout)
            nextl = next(slhaout)
            while(nextl.split()[0] != "25"):
                nextl = next(slhaout)
            mh1 = nextl.split()[1]
            break


fs_slha=FS_slha_setup+slha_sminputs+FS_model_inputs
sp_slha=SPheno_slha_setup+slha_sminputs+spheno_model_inputs
nt_slha=nmssmtools_setup+slha_sminputs+nmssmtools_model_inputs
ss_slha=softsusy_setup+slha_sminputs+nmssmtools_model_inputs
nc_slha=nmssmcalc_setup+slha_sminputs+nmssmcalc_model_inputs
ncOS_slha=nmssmcalc_setup+slha_sminputs+nmssmcalc_model_inputs

# parameters we fix
lambdaMS = 0.4
kappaMS = 0.6
tanbMS = 5.0
Xtop = 0.0

# parameters for simple MSUSY scan
startMS = 100
endMS = 30000
numMS = 1000
endXt = 8000.0
startXt = -8000.0

specgen_call = 'models/NMSSMtower/run_NMSSMtower.x'
input_option = ' --slha-input-file='
input_file = 'LesHouches.in.NMSSMtower'
#output_option = ' >& '

for i in range(0,numMS+1):
  #  MSUSY =  (endMS - startMS)  /  (numMS + 0.0) * i + startMS
    MSUSY = 2000.0
    Xtop =  startXt + (endXt - startXt)  /  (numMS + 0.0) * i

    # print "*** FOR POINT ", i , " ****"
    output_file = 'scan_outputs/NMSSMtowerNEWTEST/SLHA.BMMs' + 'Xt'+str(Xtop) +'.out'
#    cmd_line_call = specgen_call + input_option + input_file + output_option + output_file

    cmd_line_call = specgen_call + input_option + input_file

    #    print cmd_line_call
    f = open(input_file, 'w')
    f.write(fs_slha.format( Ms = MSUSY, tbMSUSY = tanbMS, lam = lambdaMS,
                            kap=kappaMS, Xt=Xtop ))
    f.close()
    import shlex, subprocess
    args = shlex.split(cmd_line_call)
    fout = open(output_file, "w")
    proc = subprocess.Popen(args, stdin = subprocess.PIPE, stdout = fout)
    stdout, stderr = proc.communicate('dir c:\\')

    slhaout = open(output_file)

    # No ERROR unless we find problem flag
    badword = "Problems"
    fs_error = 0
    TYu = [[0 for x in range(3)] for x in range(3)]
    TYe = [[0 for x in range(3)] for x in range(3)]
    TYd = [[0 for x in range(3)] for x in range(3)]
    Yu = [[0 for x in range(3)] for x in range(3)]
    Ye = [[0 for x in range(3)] for x in range(3)]
    Yd = [[0 for x in range(3)] for x in range(3)]
    Au = [[0 for x in range(3)] for x in range(3)]
    Ae = [[0 for x in range(3)] for x in range(3)]
    Ad = [[0 for x in range(3)] for x in range(3)]
    Alam = 0.0
    Akap = 0.0
    for line in slhaout:
        if badword in line :
            ourline = 'Xt'+str(Xtop) + line
            fs_error = 1
            continue
        if 'Block TanBetaAtMZ' in line:
            tbline = next(slhaout)
            if tbline.split()[0] == "1" :
                tbMZ = tbline.split()[1]
#                print "now tbMZ =", tbMZ
        if 'Block Tu' in line:
            fill_block_3b3(line, slhaout, TYu)
        if 'Block Te' in line:
            fill_block_3b3(line, slhaout, TYe)
        if 'Block Td' in line:
            fill_block_3b3(line, slhaout, TYd)
        if 'Block Yu' in line:
            fill_block_3b3(line, slhaout, Yu)
        if 'Block Ye' in line:
            fill_block_3b3(line, slhaout, Ye)
        if 'Block Yd' in line:
            fill_block_3b3(line, slhaout, Yd)
        if 'Block Au' in line:
            fill_block_3b3(line, slhaout, Au)
        if 'Block Ae' in line:
            fill_block_3b3(line, slhaout, Ae)
        if 'Block Ad' in line:
            fill_block_3b3(line, slhaout, Ad)
        if 'Block Alambda' in line:
            Alamline = next(slhaout)
            Alam = Alamline.split()[0]
        if 'Block Akappa' in line:
            Akapline = next(slhaout)
            Akap = Akapline.split()[0]
        if 'Block MASS' in line:
            nl = next(slhaout)
            while(nl.split()[0] != "25"):
                nl = next(slhaout)
            fs_mh1 = nl.split()[1]

    if fs_error == 0 :
        # print "trilinears", TYu[2][2], TYd[2][2], TYe[2][2]
        # print "Yukawas", Yu[2][2], Yd[2][2], Ye[2][2]
        # print "trilinears/yukawas ", float(TYu[2][2]) / float(Yu[2][2]), float(TYd[2][2]) / float(Yd[2][2]), float(TYe[2][2]) / float(Ye[2][2])
        # print "soft A terms       ", Au[2][2], Ad[2][2], Ae[2][2]
        # print "Alam, Akap = ", Alam, Akap
        spheno_exe = '/home/pathron/software/SPheno-3.3.6/bin/SPhenoNMSSM'
        spheno_input_file = 'LesHouches.in.NMSSM'

        SPin = open(spheno_input_file, 'w')

        # print "MSUSY =", MSUSY
        # print "tbMZ = ", tbMZ
        # print "lambdaMS = ", lambdaMS
        # print "kappaMS = ", kappaMS
        # print "At =", Au[2][2]
        # print "Ab =", Ad[2][2]
        # print "Atau =", Ae[2][2]
        SPin.write(sp_slha.format( Ms = MSUSY, tb = float(tbMZ), lam = lambdaMS,
                                   kap = kappaMS, At = float(Au[2][2]),
                                   Ab = float(Ad[2][2]),
                                   Atau = float(Ae[2][2]),
                                   Al = float(Alam), Ak=float(Akap),
                                   Tu11In = float(TYu[0][0]),
                                   Td11In = float(TYd[0][0]),
                                   Te11In = float(TYe[0][0]),
                                   Tu22In = float(TYu[1][1]),
                                   Td22In = float(TYd[1][1]),
                                   Te22In = float(TYe[1][1]),
                                   Tu33In = float(TYu[2][2]),
                                   Td33In = float(TYd[2][2]),
                                   Te33In = float(TYe[2][2])  ))

        SPin.close()
        import shutil
        spheno_input_save = "SPheno_NMSSM_scan/LesHouches.in.SLHABM"+'Xt'+str(Xtop)
        shutil.copyfile(spheno_input_file, spheno_input_save)

        proc = subprocess.Popen(spheno_exe, stdin = subprocess.PIPE, stdout = subprocess.PIPE)
        stdout, stderr = proc.communicate('dir c:\\')

        filename = "SPheno_NMSSM_scan/SLHA.BM"+ 'Xt'+str(Xtop) +'.out'
        import os
        os.rename("SPheno.spc.NMSSM", filename)

        SP_slhaout = open(filename)
        for line in SP_slhaout:
            if 'Block MASS' in line:
                #skip spheno's comment
                next(SP_slhaout)
                nextl = next(SP_slhaout)
                while(nextl.split()[0] != "25"):
                    nextl = next(SP_slhaout)
                sp_mh1 = nextl.split()[1]
                break

        nmssmtools_exe = '/home/pathron/software/NMSSMTools_4.8.2/run'
        nmssmtools_input_file = '../tom_steudtner/nmssmtools_scan/SLHA.'+'Xt'+str(Xtop)+'.inp.dat'
        nmssmtools_output_file = '../tom_steudtner/nmssmtools_scan/SLHA.'+'Xt'+str(Xtop)+'.spectr.dat'
        NTin = open(nmssmtools_input_file, 'w')

        # print "MSUSY =", MSUSY
        # print "tbMZ = ", tbMZ
        # print "lambdaMS = ", lambdaMS
        # print "kappaMS = ", kappaMS
        # print "At =", Au[2][2]
        # print "Ab =", Ad[2][2]
        # print "Atau =", Ae[2][2]
        NTin.write(nt_slha.format( Ms = MSUSY, tb = float(tbMZ), lam = lambdaMS,
                                   kap = kappaMS, At = float(Au[2][2]),
                                   Ab = float(Ad[2][2]),
                                   Atau = float(Ae[2][2]),
                                   Amu = float(Ae[1][1]),
                                   Al = float(Alam), Ak=float(Akap)))

        NTin.close()

        nt_cmd_line = nmssmtools_exe + ' ' + nmssmtools_input_file
        nt_args = shlex.split(nt_cmd_line)
        os.chdir('/home/pathron/software/NMSSMTools_4.8.2/')

        proc = subprocess.Popen(nt_args, stdin = subprocess.PIPE,
                                stdout = subprocess.PIPE)
        stdout, stderr = proc.communicate('dir c:\\')

        os.chdir('/home/pathron/software/tom_steudtner')
        NT_slhaout = open(nmssmtools_output_file)

        for line in NT_slhaout :
            if 'BLOCK MASS' in line:
                #skip spheno's comment
                next(NT_slhaout)
                nextl = next(NT_slhaout)
                while(nextl.split()[0] != "25"):
                    nextl = next(NT_slhaout)
                nt_mh1 = nextl.split()[1]
                break


        ss_call = '/home/pathron/software/softsusy-3.6.2/softpoint.x'
        ss_input_file ='softsusy_scan/NMSSM/SLHA.BMMs' + 'Xt'+str(Xtop) +'.in'
        ss_output_file ='softsusy_scan/NMSSM/SLHA.BMMs' + 'Xt'+str(Xtop) +'.out'
        input_options = ' leshouches'
        ss_cmd_line = ss_call + input_options

        SSin = open(ss_input_file, 'w')
        SSin.write(ss_slha.format( Ms = MSUSY, tb = float(tbMZ), lam = lambdaMS,
                                   kap = kappaMS, At = float(Au[2][2]),
                                   Ab = float(Ad[2][2]),
                                   Atau = float(Ae[2][2]),
                                   Amu = float(Ae[1][1]),
                                   Al = float(Alam), Ak=float(Akap)))
        SSin.close()
        ss_args = shlex.split(ss_cmd_line)
        ss_out = open(ss_output_file, "w")
        ss_in = open(ss_input_file)
        proc = subprocess.Popen(ss_args, stdin = ss_in, stdout = ss_out)
        stdout, stderr = proc.communicate('dir c:\\')

        SS_slhaout = open(ss_output_file)

        for line in SS_slhaout :
            if 'Block MASS' in line:
                #skip spheno's comment
                next(SS_slhaout)
                nextl = next(SS_slhaout)
                while(nextl.split()[0] != "25"):
                    nextl = next(SS_slhaout)
                ss_mh1 = nextl.split()[1]
                break

        print Xtop, "  fs_mh1= ", fs_mh1, "  sp_mh1= ", sp_mh1, "nt_mh1= ", nt_mh1, "ss_mh1=", ss_mh1
