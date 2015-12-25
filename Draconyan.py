### Ulrike Hager, Apr 2010 ###

from math import sqrt
from RunSet import Set, Run
import GetMass
try:
    import matplotlib.pylab as plt
except ImportError:
    has_matplotlib = False
else: has_matplotlib = True


class Draconyan(object):
    """Main analysis class with functions to read and write files and manage runs and sets.

    General setup from DraconyanIn.dat. Most variables have default values.
    String representation of projectile and target has to be in name format given in nubtab03.asc
    """
    def __init__(self):
        self.step = 0
        """Where to start.

        If step = 0, the program will start from the beginning
        if step = 1, the program will load a file (file name \'runFile\' required in DraconyanIn.dat) to read in the runs and continue the calculation from there. This file has to have the format of DraconyanOut_runs.dat
        if step = 2, the program will load a file (file name \'setFile\' required in DraconyanIn.dat) to read in the sets and continue the calculation from there. This file has to have the format of DraconyanOut_sets.dat"""
        self.qProj = 0
        """Charge state of incoming beam."""
        self.stoppingFactor = 0
        """Factor to convert stopping power from keV/(ug*cm^2) to eV/(10^15/cm^2)."""
        self.mcpEff = [0.79,0.03]
        """MCP efficiency and uncertainty"""
        self.mcpTrans = [0.95**2,0.05]
        """MCP transmission and uncertainty.

        The nominal transmission of the mesh is 95%, squared for 2 MCPs."""
        self.deltaE = 0.7
        """Maximum difference in incoming beam energy [keV/u] for two runs to be assigned to the same set."""
        self.deltaP = 1
        """Maximum difference in pressure for two runs to be assigned to the same set."""
        self.dataFileName = "runs.dat"
        """File containing run data
    
        Format:
        run    Ein     Eout    pressure      FC4before [eA]    FC4after [eA]    target  +/-  separator    +/-    CSF     +/-"""
        self.sbFileName = "Sb_lt.dat"
        """File containing SB and lifetime data.

        Format:
        Run    SB0   +/-   SB1   +/-   lt_HI  +/-   lt_BGO   +/-    SB0start/s   +/-     SB0end/s    +/-"""
        self.geantFileName = "geant.dat"
        """File containing GEANT transmission and BGO efficiency.

        Format:
        Ein[keV/u]   pTar[T]    transm  +/-    BGOeff  +/-"""
        self.recoilFileName = "recoils.dat"
        """File containing recoil data

        Format:
        Ein[keV/u]   pTar[T]    singles    coincidences"""
        self.setFileName = ""
        """File containing set data.

        Needed if step = 2, in which case the other input files are no longer needed.
        Format:
        Runs   Ein   Eout   Ecm  +/-  pressure  TargetDensity   +/-  stopping  +/-  tarTrans  +/-   sepTrans   +/-   HI-lt   +/-   CSF   +/-   MCPeff   +/-   geantTransm   +/-   geantBGOeff   +/-   SB0tot   +/-   R   +/-   beam   +/-   ( R1   +/-     FC4/SB0   +/-     FC4/SB1   +/- ( singles    coincidences  ))"""
        self.runFileName = ""
        """File containing pre-calculated run data.

        Needed if step = 1, the recoils file and the geant file are still needed.
        Format:
        Run   Ein   Eout   Ecm   +/-   pressure   FC4before   +/-   FC4after   +/-    tarTrans   +/-   sepTrans   +/-   HI_lt   +/-   CSF   +/-   SB0tot   +/-   SB0start   +/-   SB0end   +/-   Rstart   +/-   Rend   +/-   FC4/SB0start   +/-   FC4/SB0end   +/-   SB1start   +/-   SB1end   +/-   R1start   +/-   R1end   +/-   FC4/SB1start   +/-   FC4/SB1end   +/-"""
        self.targetLength = [12.3,0.5]
        """Length of DRAGON target in cm with uncertainty."""
        self.sets = []
        """List of sets of type Set"""
        self.runs = []
        """List of runs of type Run"""
        self.upperLimitSwitch = 1
        """Decides whether or not to use the Feldman & Cousins limits for small numbers

        \'1\' uses Feldman & Cousins for sets with up to 15 events, \'0\' uses normal error propagation for all sets."""

    def calc_beam_in_target(self):
        """Calls calc_beam_in_target for all sets."""
        for set in self.sets:
            set.calc_beam_in_target()

    def calc_beam_in_target_sb1(self):
        """Calls calc_beam_in_target for all sets."""
        for set in self.sets:
            set.calc_beam_in_target_sb1()

    def calc_Ecm(self):
        """Calls \'calc_Ecm\' for all runs and sets."""
        for run in self.runs:
            run.calc_Ecm(self.projectile['massU'], self.target['massU'])
        for set in self.sets:
            set.calc_Ecm(self.projectile['massU'], self.target['massU'])

    def calc_elastics(self):
        """Calculates elastics for SB0 and SB1 both Rutherford corrected (R, R1) and not corrected (IoverN, IoverN1).

        To calculate only the normalisation factor R, use calc_R_set()"""
        for set in self.sets:
            set.calc_elastics()
    
    def convert_fc4(self):
        """Calls convert_fc4 for all runs."""
        for run in self.runs:
            run.convert_fc4(self.qProj)

    def calc_IoverN_run(self):
        """Calls calc_IoverN and IoverN1 for all runs."""
        for run in self.runs:
            run.calc_IoverN()
            run.calc_IoverN1()

    def calc_R_run(self):
        """Calls calc_R for all runs."""
        for run in self.runs:
            run.calc_R()

    def calc_R_set(self):
        """Calls calc_R for all sets.

        For normalisation, this is sufficient. To analyse the elastics also using SB1, call calc_elastics()."""
        for set in self.sets:
            set.calc_R()

    def calc_R1_run(self):
        """Calls calc_R1 for all runs."""
        for run in self.runs:
            run.calc_R1()

    def calc_R1_set(self):
        """Calls calc_R1 for all sets.

        Normalisation factor R1 based on SB1 data. Not actually used for normalisation, use R (calc_R_set()) based on SB0 for that.
        To calculate both R and R1, and the ratios FC4/SBx, call calc_elastics()."""
        for set in self.sets:
            set.calc_R1()

    def calc_set_averages(self):
        """Calls calc_averages for all sets."""
        for set in self.sets:
            set.calc_averages()

    def sort_set_runNumbers(self):
        """Sorts run numbers in all sets."""
        for set in self.sets:
            set.runNumber.sort()

    def calc_Sfactors(self):
        """Calls calc_Sfactors for all sets.

        calc_sigmas has to be calles first."""
        reducedMass = (self.projectile['massU'] * self.target['massU']) / (self.projectile['massU'] + self.target['massU'])
        for set in self.sets:
            set.calc_Sfactor(self.projectile['Z'], self.target['Z'], reducedMass)
    
    def calc_sigmas(self):
        """Calls calc_sigma for all sets."""
        for set in self.sets:
            set.calc_sigma(self.targetLength)
            #print set.reacSigma

    def calc_stopping_powers(self):
        """Calls calc_stopping_power for all sets."""
        for set in self.sets:
            set.calc_stopping_power(self.projectile['A'], self.target['massKg'],self.targetLength)

    def calc_target_densities(self):
        """Calls calc_target_density for all sets."""
        for set in self.sets:
            set.calc_target_density(self.nuTarget)

    def calc_transmissions(self):
        """Calls calc_transmission for all sets."""
        for set in self.sets:
            set.calc_transmission()

    def calc_wgs(self):
        """Calls calc_wg for all sets."""
        for set in self.sets:
            check = set.calc_wg(self.projectile['A'], self.target['A'], self.stoppingFactor)
            if check == 1:
                print "Conversion factor for stopping power is 0, cannot calculate omega-gamma."
                return
 
    def calc_yields(self):
        """Calls calc_yield for all sets."""
        for set in self.sets:
            set.calc_yield()

    def determine_nuTarget(self):
        if self.target['name'] == '4He': self.nuTarget = 1
        elif self.target['name'] == '1H': self.nuTarget = 2
        elif self.target['name'] == '2H': self.nuTarget = 2
        else: self.nuTarget = 1

    def find_run_number(self,number):
        """Returns the run by runNumber, if not found, returns 0."""
        for run in self.runs:
            if run.runNumber == number: return run
        return 0

    def find_set(self,Ein,pTar,dE=0.5,dP=0.5):
        """Finds a set of runs based on energy and pressure.
        
        \'Ein\' is the incoming beam energy in keV/u, \'pTar\' the target pressure in torr
        \'dE\' and \'dP\' are the maximum deviations in energy and pressure allowed to be considered the same set."""
        for set in self.sets:
            if ( abs(set.eBeam[0]-Ein) < dE and abs(set.pTarget-pTar) < dP ): return set
        return 0

    def initialise_sets(self):
        """Sort runs into sets depending on incoming beam energy and presure.

        \'deltaE\' and \'deltaP\' are the maximum deviations in energy (keV/u) and pressure (torr) for a run to be considered as belonging to a set."""
        for run in self.runs:
            for set in self.sets:
                if ( abs(set.eBeam[0]-run.eBeam[0]) < self.deltaE and abs(set.pTarget-run.pTarget) < self.deltaP ):
                    set.runNumber.append(run.runNumber)
                    set.RList.append(run.R[0])
                    set.RList.append(run.R[1])
                    set.R1List.append(run.R1[0])
                    set.R1List.append(run.R1[1])
                    set.IoNList.append(run.IoverN[0])
                    set.IoNList.append(run.IoverN[1])
                    set.IoN1List.append(run.IoverN1[0])
                    set.IoN1List.append(run.IoverN1[1])
                    set.sb0Tot[0] += run.sb0Tot[0]
                    set.sb0Tot[1] = sqrt(set.sb0Tot[0])
                    set.sb1Tot[0] += run.sb1Tot[0]
                    set.sb1Tot[1] = sqrt(set.sb1Tot[0])
                    break
            else:
                set = Set().make_set_from_run(run)
                set.mcpEff = self.mcpEff
                set.mcpTrans = self.mcpTrans
                self.sets.append(set)

    def plot_individual_IoverN(self):
        """Plots FC4/SB0 for all runs.

        Requires matplotlib
        Can only run if running from the beginning, i.e. step = 0"""
        if has_matplotlib == False:
            print "Cannot use this function, matplotlib not available"
            return
        if len(self.runs) == 0:
            print "No run data available for plotting"
            return
        ## rEndList = []
        ## rStartList = []
        rList = []
        for run in self.runs:
            ## rStartList.append([run.runNumber-0.2,run.R[0][0],run.R[0][1]])
            ## rEndList.append([run.runNumber+0.2,run.R[1][0],run.R[1][1]])
            rList.append([run.runNumber-0.2,run.IoverN[0][0],run.IoverN[0][1]])
            rList.append([run.runNumber+0.2,run.IoverN[1][0],run.IoverN[1][1]])
        fig = plt.figure()
        ## eStart = fig.add_subplot(1,1,1)
        ## eStart.scatter([x for [x,y,z] in rStartList],[y for [x,y,z] in rStartList],s=20,alpha=0.7)
        ## eStart.errorbar([x for [x,y,z] in rStartList],[y for [x,y,z] in rStartList],[z for [x,y,z] in rStartList])
        ## eEnd = fig.add_subplot(1,1,1)
        ## eEnd.scatter([x for [x,y,z] in rEndList],[y for [x,y,z] in rEndList],s=20,marker=">",alpha=0.7,color='g')
        ## eEnd.errorbar([x for [x,y,z] in rEndList],[y for [x,y,z] in rEndList],[z for [x,y,z] in rEndList],ecolor='g')
        eComb = fig.add_subplot(1,1,1)
        eComb.scatter([x for [x,y,z] in rList],[y for [x,y,z] in rList],s=20,alpha=0.7)
        eComb.errorbar([x for [x,y,z] in rList],[y for [x,y,z] in rList],[z for [x,y,z] in rList])
        plt.xlabel('run number')
        plt.ylabel('FC4/SB0')
        plt.title('FC4/SB0 normalisation')
        
        plt.show()

    def plot_individual_IoverN1(self):
        """Plots FC4/SB0 for all runs.

        Requires matplotlib
        Can only run if running from the beginning, i.e. step = 0"""
        if has_matplotlib == False:
            print "Cannot use this function, matplotlib not available"
            return
        if len(self.runs) == 0:
            print "No run data available for plotting"
            return
        ## rEndList = []
        ## rStartList = []
        rList = []
        for run in self.runs:
            ## rStartList.append([run.runNumber-0.2,run.R[0][0],run.R[0][1]])
            ## rEndList.append([run.runNumber+0.2,run.R[1][0],run.R[1][1]])
            rList.append([run.runNumber-0.2,run.IoverN1[0][0],run.IoverN1[0][1]])
            rList.append([run.runNumber+0.2,run.IoverN1[1][0],run.IoverN1[1][1]])
        fig = plt.figure()
        ## eStart = fig.add_subplot(1,1,1)
        ## eStart.scatter([x for [x,y,z] in rStartList],[y for [x,y,z] in rStartList],s=20,alpha=0.7)
        ## eStart.errorbar([x for [x,y,z] in rStartList],[y for [x,y,z] in rStartList],[z for [x,y,z] in rStartList])
        ## eEnd = fig.add_subplot(1,1,1)
        ## eEnd.scatter([x for [x,y,z] in rEndList],[y for [x,y,z] in rEndList],s=20,marker=">",alpha=0.7,color='g')
        ## eEnd.errorbar([x for [x,y,z] in rEndList],[y for [x,y,z] in rEndList],[z for [x,y,z] in rEndList],ecolor='g')
        eComb = fig.add_subplot(1,1,1)
        eComb.scatter([x for [x,y,z] in rList],[y for [x,y,z] in rList],s=20,alpha=0.7)
        eComb.errorbar([x for [x,y,z] in rList],[y for [x,y,z] in rList],[z for [x,y,z] in rList])
        plt.xlabel('run number')
        plt.ylabel('FC4/SB0')
        plt.title('FC4/SB0 normalisation')
        
        plt.show()

    def plot_individual_Rs(self):
        """Plots normalisation factor R for all runs.

        Requires matplotlib
        Can only run if running from the beginning, i.e. step = 0"""
        if has_matplotlib == False:
            print "Cannot use this function, matplotlib not available"
            return
        if len(self.runs) == 0:
            print "No run data available for plotting"
            return
        ## rEndList = []
        ## rStartList = []
        rList = []
        for run in self.runs:
            ## rStartList.append([run.runNumber-0.2,run.R[0][0],run.R[0][1]])
            ## rEndList.append([run.runNumber+0.2,run.R[1][0],run.R[1][1]])
            rList.append([run.runNumber-0.2,run.R[0][0],run.R[0][1]])
            rList.append([run.runNumber+0.2,run.R[1][0],run.R[1][1]])
        fig = plt.figure()
        ## eStart = fig.add_subplot(1,1,1)
        ## eStart.scatter([x for [x,y,z] in rStartList],[y for [x,y,z] in rStartList],s=20,alpha=0.7)
        ## eStart.errorbar([x for [x,y,z] in rStartList],[y for [x,y,z] in rStartList],[z for [x,y,z] in rStartList])
        ## eEnd = fig.add_subplot(1,1,1)
        ## eEnd.scatter([x for [x,y,z] in rEndList],[y for [x,y,z] in rEndList],s=20,marker=">",alpha=0.7,color='g')
        ## eEnd.errorbar([x for [x,y,z] in rEndList],[y for [x,y,z] in rEndList],[z for [x,y,z] in rEndList],ecolor='g')
        eComb = fig.add_subplot(1,1,1)
        eComb.scatter([x for [x,y,z] in rList],[y for [x,y,z] in rList],s=20,alpha=0.7)
        eComb.errorbar([x for [x,y,z] in rList],[y for [x,y,z] in rList],[z for [x,y,z] in rList])
        plt.xlabel('run number')
        plt.ylabel('R')
        plt.title('Normalisation factor R')
        
        plt.show()

    def plot_individual_R1s(self):
        """Plots normalisation factor R for all runs.

        Requires matplotlib
        Can only run if running from the beginning, i.e. step = 0"""
        if has_matplotlib == False:
            print "Cannot use this function, matplotlib not available"
            return
        if len(self.runs) == 0:
            print "No run data available for plotting"
            return
        ## rEndList = []
        ## rStartList = []
        rList = []
        for run in self.runs:
            ## rStartList.append([run.runNumber-0.2,run.R[0][0],run.R[0][1]])
            ## rEndList.append([run.runNumber+0.2,run.R[1][0],run.R[1][1]])
            rList.append([run.runNumber-0.2,run.R1[0][0],run.R1[0][1]])
            rList.append([run.runNumber+0.2,run.R1[1][0],run.R1[1][1]])
        fig = plt.figure()
        ## eStart = fig.add_subplot(1,1,1)
        ## eStart.scatter([x for [x,y,z] in rStartList],[y for [x,y,z] in rStartList],s=20,alpha=0.7)
        ## eStart.errorbar([x for [x,y,z] in rStartList],[y for [x,y,z] in rStartList],[z for [x,y,z] in rStartList])
        ## eEnd = fig.add_subplot(1,1,1)
        ## eEnd.scatter([x for [x,y,z] in rEndList],[y for [x,y,z] in rEndList],s=20,marker=">",alpha=0.7,color='g')
        ## eEnd.errorbar([x for [x,y,z] in rEndList],[y for [x,y,z] in rEndList],[z for [x,y,z] in rEndList],ecolor='g')
        eComb = fig.add_subplot(1,1,1)
        eComb.scatter([x for [x,y,z] in rList],[y for [x,y,z] in rList],s=20,alpha=0.7)
        eComb.errorbar([x for [x,y,z] in rList],[y for [x,y,z] in rList],[z for [x,y,z] in rList])
        plt.xlabel('run number')
        plt.ylabel('R(SB1)')
        plt.title('Normalisation factor R using SB1')
        
        plt.show()


    def plot_yields(self): 
        """Plots yields for all sets

        Requires matplotlib."""
        if has_matplotlib == False:
            print "Cannot use plotting, matplotlib not available"
            return
        if len(self.sets) == 0:
            print "No set data available for plotting"
            return
        yieldL = []
        for set in self.sets: yieldL.append([set.Ecm[0],set.reacYield[0],set.reacYield[1]])
        #print yieldL
        fig = plt.figure() 
        yieldS = fig.add_subplot(1,1,1)
        yieldS.semilogy([E for [E,[Ys,dYs,dsysYs],[Yc,dYc,dsysYc]] in yieldL], [Ys for [E,[Ys,dYs,dsysYs],[Yc,dYc,dsysYc]] in yieldL],color='b')
        yieldS.errorbar([E for [E,[Ys,dYs,dsysYs],[Yc,dYc,dsysYc]] in yieldL], [Ys for [E,[Ys,dYs,dsysYs],[Yc,dYc,dsysYc]] in yieldL], [dYs for [E,[Ys,dYs,dsysYs],[Yc,dYc,dsysYc]] in yieldL],ecolor='b')
        yieldC = fig.add_subplot(1,1,1)
        yieldC.semilogy([E for [E,[Ys,dYs,dsysYs],[Yc,dYc,dsysYc]] in yieldL], [Yc for [E,[Ys,dYs,dsysYs],[Yc,dYc,dsysYc]] in yieldL],color='g')
        yieldC.errorbar([E for [E,[Ys,dYs,dsysYs],[Yc,dYc,dsysYc]] in yieldL], [Yc for [E,[Ys,dYs,dsysYs],[Yc,dYc,dsysYc]] in yieldL], [dYc for [E,[Ys,dYs,dsysYs],[Yc,dYc,dsysYc]] in yieldL],ecolor='g')
        plt.xlabel('E_cm [MeV]')
        plt.ylabel('yield')
        plt.title('Reaction yield')
        plt.show()

    def plot_wg_coinc(self): 
        """Plots wgs for all sets

        Requires matplotlib."""
        if has_matplotlib == False:
            print "Cannot use plotting, matplotlib not available"
            return
        if len(self.sets) == 0:
            print "No set data available for plotting"
            return
        wgL = []
        for set in self.sets: wgL.append([set.Ecm[0],set.reacWg[0],set.reacWg[1]])
        #print wgL
        fig = plt.figure() 
        wgC = fig.add_subplot(1,1,1)
        wgC.scatter([E for [E,[Ys,dYs,dsysYs],[Yc,dYc,dsysYc]] in wgL], [Yc for [E,[Ys,dYs,dsysYs],[Yc,dYc,dsysYc]] in wgL],color='g')
        #wgC.semilogy([E for [E,[Ys,dYs,dsysYs],[Yc,dYc,dsysYc]] in wgL], [Yc for [E,[Ys,dYs,dsysYs],[Yc,dYc,dsysYc]] in wgL],color='g')
        wgC.errorbar([E for [E,[Ys,dYs,dsysYs],[Yc,dYc,dsysYc]] in wgL], [Yc for [E,[Ys,dYs,dsysYs],[Yc,dYc,dsysYc]] in wgL], [dYc for [E,[Ys,dYs,dsysYs],[Yc,dYc,dsysYc]] in wgL],ecolor='g')
        plt.xlabel('E_cm [MeV]')
        plt.ylabel('wg')
        plt.title('Reaction wg')
        plt.show()

    def read_geantdata(self):
        """Read GEANT transmission and BGO efficiency data from file.

        Format:
        Ein[keV/u]   pTar[T]    transm  +/-    BGOeff  +/-
        """
        try:
            dataFile = open(self.geantFileName,'r')
        except IOError:
            print "could not open file", self.geantFileName 
            return 1
        else:
            for line in dataFile:
                if (line[0] == '#' or line == '\n'): continue
                try:
                    Ein = float(line.split()[0])
                    pTar = float(line.split()[1])
                    set = self.find_set(Ein,pTar,self.deltaE,self.deltaP)
                    if set != 0: set.get_geant_data(line)
                except IndexError:
                    print self.geantFileName , "has wrong format, required format:\nEin   pTar    transm  +/-    BGOeff  +/-"
                    dataFile.close()
                    return 1
            dataFile.close()
            return 0

    def read_input_file(self,fileName='DraconyanIn.dat'):
        try:
            dataFile = open(fileName,'r')
        except:
            print "Could not read input file", fileName
            exit(1)
        else:
            for line in dataFile:
                if (line[0] == '#' or line == '\n'): continue
                lList = line.split()
                if lList[0] == 'projectile': projectileName = lList[1]
                elif lList[0] == 'target': targetName = lList[1]
                elif lList[0] == 'projChargeState':  self.qProj = int(lList[1])
                elif lList[0] == 'stoppingFactor': self.stoppingFactor = float(lList[1])
                elif lList[0] == 'MCPeff':
                    self.mcpEff[0] = float(lList[1])
                    self.mcpEff[1] = float(lList[2])
                elif lList[0] == 'MCPtrans':
                    self.mcpTrans[0] = float(lList[1])
                    self.mcpTrans[1] = float(lList[2])
                elif lList[0] == 'deltaE': self.deltaE = float(lList[1])
                elif lList[0] == 'deltaP': self.deltaP = float(lList[1])
                elif lList[0] == 'dataFile': self.dataFileName = lList[1]
                elif lList[0] == 'geantFile': self.geantFileName = lList[1]
                elif lList[0] == 'recoilFile': self.recoilFileName = lList[1]
                elif lList[0] == 'SbFile': self.sbFileName = lList[1]
                elif lList[0] == 'setFile': self.setFileName = lList[1]
                elif lList[0] == 'runFile': self.runFileName = lList[1]
                elif lList[0] == 'step':  self.step = int(lList[1])
                elif lList[0] == 'targetLength':
                    self.targetLength[0] = float(lList[1])
                    self.targetLength[1] = float(lList[2])
                elif lList[0] == 'dEnergy': self.dEnergy = float(lList[1])
                elif lList[0] == 'dPressure': self.dPressure = float(lList[1])
                elif lList[0] == 'dFC4': self.dFc4 = float(lList[1])
                elif lList[0] == 'FeldmanCousinsLimits': self.upperLimitSwitch = int(lList[1])
            dataFile.close()
            self.projectile = GetMass.read_nubtab(projectileName)
            self.target = GetMass.read_nubtab(targetName)
            if (self.projectile == 0 or self.target == 0):
                print "Could not find projectile or target"
                exit(1)
            elif self.qProj == 0:
                print "Could not find charge state of projectile"
                exit(1)
            
        if self.step ==1:
            fileCheck = self.read_runs_from_file_step1()
            if fileCheck == 1:
                print "Could not run from step 1, running from the beginning"
                self.step = 0   
        elif self.step == 2:
            fileCheck = self.read_sets_from_file_step2()
            if fileCheck == 1:
                print "Could not run from step 2, running from the beginning"
                self.step = 0
        return self.step

    def read_recoildata(self):
        """Read recoil data from file.

        Format:
        Ein[keV/u]   pTar[T]    singles    coincidences"""
        try:
            dataFile = open(self.recoilFileName,'r')
        except IOError:
            print "could not open file", self.recoilFileName 
            return 1
        else:
            for line in dataFile:
                if (line[0] == '#' or line == '\n'): continue
                try:
                    Ein = float(line.split()[0])
                    pTar = float(line.split()[1])
                    set = self.find_set(Ein,pTar,self.deltaE, self.deltaP)
                    if set != 0: set.get_recoil_data(line, self.upperLimitSwitch)
                except IndexError:
                    print self.recoilFileName , "has wrong format, required format:\nEin   pTar  singles   coincidences"
                    dataFile.close()
                    return 1
            dataFile.close()
            return 0

    def read_rundata(self):
        """Reads run data from file.

        Format:
        run    Ein     Eout    pressure      FC4before [eA]    FC4after [eA]    target  +/-  separator    +/-    CSF     +/-
        and creates objects of type Run."""
        try:
            dataFile = open(self.dataFileName, 'r')
        except IOError:
            print "could not open file", self.dataFileName , ", required format\nrun    Ein     Eout    pressure      FC4before [eA]    FC4after [eA]    target  +/-  separator    +/-    CSF     +/-"
            return 1
        else:
            for line in dataFile:
                if (line[0] == '#' or line == '\n'): continue
                try:
                    self.runs.append(Run().make_run_from_input(line))
                except IndexError:
                    print self.dataFileName , "has wrong format, required format:\nrun    Ein     Eout    pressure      FC4before [eA]    FC4after [eA]    target  +/-  separator    +/-    CSF     +/-"
                    dataFile.close()
                    return 1
            dataFile.close()
            return 0
        
    def read_runs_from_file_step1(self):
        try:
            dataFile = open(self.runFileName, 'r')
        except IOError:
            print "could not open file", self.runFileName 
            return 1
        else:
            for line in dataFile:
                if (line[0] == '#' or line == '\n') : continue
                try:
                    self.runs.append(Run().make_run_from_file_step1(line))
                except IndexError:
                    print self.runFileName , "has wrong format, required format:\nRun   Ein   Eout   Ecm   +/-   pressure   FC4before   +/-   FC4after   +/-    tarTrans   +/-   sepTrans   +/-   HI_lt   +/-   CSF   +/-   SB0tot   +/-   SB0start   +/-   SB0end   +/-   Rstart   +/-   Rend   +/-   FC4/SB0start   +/-   FC4/SB0end   +/-   SB1start   +/-   SB1end   +/-   R1start   +/-   R1end   +/-   FC4/SB1start   +/-   FC4/SB1end   +/-"
                    dataFile.close()
                    return 1
            dataFile.close()
            return 0
            
    def read_sbdata(self):
        """Read SB and lt data from file.

        Format:
        Run    SB0   +/-   SB1   +/-   lt_HI  +/-   lt_BGO   +/-    SB0start/s   +/-     SB0end/s    +/-    SB1start/s   +/-     SB1end/s    +/- """
        try:
            dataFile = open(self.sbFileName, 'r')
        except IOError:
            print "could not open file", self.sbFileName , ", required format\nRun    SB0   +/-   SB1   +/-   lt_HI  +/-   lt_BGO   +/-    SB0start/s   +/-     SB0end/s    +/-    SB1start/s   +/-     SB1end/s    +/- "
            return 1
        else:
            for line in dataFile:
                if (line[0] == '#' or line == '\n'): continue
                try:
                    run = self.find_run_number(int(line.split()[0]))
                    if run != 0: run.get_sb_data(line)
                except IndexError:
                    print self.sbFileName , "has wrong format, required format:\nRun    SB0   +/-   SB1   +/-   lt_HI  +/-   lt_BGO   +/-    SB0start/s   +/-     SB0end/s    +/-    SB1start/s   +/-     SB1end/s    +/- "
                    dataFile.close()
                    return 1
            dataFile.close 
            return 0

    def read_sets_from_file_step2(self):
        try:
            dataFile = open(self.setFileName, 'r')
        except IOError:
            print "could not open file", self.setFileName 
            return 1
        else:
            for line in dataFile:
                if (line[0] == '#' or line == '\n') : continue
                try:
                    self.sets.append(Set().make_set_from_set_step2(line,self.upperLimitSwitch))
                except IndexError:
                    print self.setFileName , "has wrong format, required format:\nRuns   Ein   Eout   Ecm  +/-  pressure  TargetDensity   +/-  stopping  +/-  tarTrans  +/-   sepTrans   +/-   HI-lt   +/-   CSF   +/-   MCPeff   +/-   geantTransm   +/-   geantBGOeff   +/-   SB0tot   +/-   R   +/-   beam   +/-  ( R1   +/-     FC4/SB0   +/-     FC4/SB1   +/- ( singles   coincidences  ))"
                    dataFile.close()
                    return 1
            dataFile.close()
            return 0
            
    def return_run_list(self):
        runList = []
        for run in self.runs:
            runList.append(run.make_input_from_run())
        return runList

    def sort_sets_by_Ecm(self):
        """Sort sets by cm energy for nicer plots."""
        self.sets.sort(key=lambda x: x.Ecm[0])
        
    def write_file_step3(self):
        """Write file for sets with the results.

        Creates output file DraconyanOut_results.dat with yields and cross sections"""
        outFile = open("DraconyanOut_results.dat",'w')
        outFile.write(self.sets[0].summarise_results(1))
        for set in self.sets:
            outFile.write(set.summarise_results())
        outFile.close()

    def write_runfile_step1(self):
        """Write overview file for runs.

        Creates output files DraconianOut_runs.dat with the normalisation data"""
        outFile = open("DraconyanOut_runs.dat",'w')
        outFile.write(self.runs[0].get_normalisation(1))
        for run in self.runs:
            outFile.write(run.get_normalisation())
        outFile.close()

    def write_setfile_step2(self):
        """Write overview file for sets.

        Creates output files DraconianOut_sets.dat with the normalisation data.
        This format is used in the file needed when starting with step 1"""
        outFile = open("DraconyanOut_sets.dat",'w')
        outFile.write(self.sets[0].get_normalisation(1))
        for set in self.sets:
            outFile.write(set.get_normalisation())
        outFile.close()

