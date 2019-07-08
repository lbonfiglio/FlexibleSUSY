#!/usr/bin/env python
#import CNMSSM input parameters from previous scan_outputs

specgen_call = 'models/CNSSMEFTHiggs/run_CNMSSMEFTHiggs.x' #calling the spectrum generator
input_option = ' --slha-input-file='
input_file = 'models/CNMSSMEFTHiggs/Misc/LesHouches.in.CNMSSMEFTHiggs'


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
