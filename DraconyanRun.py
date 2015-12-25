### Ulrike Hager, Apr 2010 ###

import Draconyan

drac = Draconyan.Draconyan()
stepCheck = drac.read_input_file()

drac.determine_nuTarget()

if stepCheck == 0:
    dataCheck = drac.read_rundata()
    if dataCheck == 1: exit(1)

    if len(drac.runs)==0:
        print "no runs found"
        exit()
    
    sbCheck = drac.read_sbdata()
    if sbCheck == 1: exit(1)

    drac.calc_Ecm()
    drac.convert_fc4()
    drac.calc_R_run()
    drac.calc_R1_run()
    drac.calc_IoverN_run()

    drac.write_runfile_step1()
        
elif stepCheck==1:
    if (len(drac.runs)) == 0:
        print "no runs found"
        exit(1)
    
elif stepCheck == 2:
    if (len(drac.sets)) == 0:
        print "no sets found"
        exit(1)
    drac.calc_beam_in_target()
    recoilCheck = 0
    
if (stepCheck==0 or stepCheck==1):
    drac.initialise_sets()
    drac.calc_set_averages()
    drac.sort_set_runNumbers()
    #drac.calc_R_set()
    drac.calc_elastics()
    drac.calc_beam_in_target()
    drac.calc_target_densities()
    drac.calc_stopping_powers()
    geantCheck = drac.read_geantdata()
    if geantCheck == 1:
        print "Error reading GEANT data, using default values: Transm = 0.97+/-0.03, BGOeff = 0.6+/-0.1"
    recoilCheck = drac.read_recoildata()

drac.sort_sets_by_Ecm()
drac.write_setfile_step2()

drac.calc_transmissions()

drac.plot_individual_Rs()
drac.plot_individual_R1s()

if recoilCheck == 1:
    print "no recoil data"
    drac.write_file_step3()
    exit()

drac.calc_yields()
drac.calc_sigmas()
drac.calc_wgs()
drac.calc_Sfactors()

drac.write_file_step3()

#drac.plot_individual_IoverN()
#drac.plot_yields()
#drac.plot_wg_coinc()
