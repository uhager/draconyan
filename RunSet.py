### Ulrike Hager, Apr 2010 ###

from __future__ import division
from math import sqrt, isnan, exp

"""limits from Feldman&Cousins paper on low signals, use for sets with 15 or fewer events"""
lowStatLimits = ((0,1.29),(0.37,2.75),(0.74,4.25),(1.10,5.30),(2.34,6.78),(2.75,7.81),(3.82,9.28),(4.25,10.30),(5.30,11.32),(6.33,12.79),(6.78,13.81),(7.81,14.82),(8.83,16.29),(9.28,17.30),(10.30,18.32),(11.32,19.32))

class RunSettings:
    """Some values needed both for Set and Run, which are derived.

    Just to avoid having the code for the initialisation of those values twice."""
    def __init__(self):
        self.eBeam = [0,0]
        """Beam energy before and after the target in keV/u."""
        self.pTarget = 0
        """target pressure in torr"""
        self.Ecm = [0,0]
        """ Centre-of-mass energy at target centre in MeV."""
        self.tarTrans = [1,0]
        """Target transmission and uncertainty."""
        self.sepTrans = [1,0]
        """Separator transmission and uncertainty."""
        self.sb0Tot = [0,0]        
        self.sb1Tot = [0,0]        
        self.R = [0,0,0]
        """SB0 vs FC4 normalisation factor and statistic and systematic uncertainties."""  
        self.R1 = [0,0,0]
        """SB1 vs FC4 normalisation factor and uncertainty."""  
        self.dPressure = 2
        """Uncertainties of the target pressure in percent."""
        self.dEnergy = 1
        """Uncertainty of the beam energy in keV/u."""
        self.dFc4 = 3
        """FC4 reading uncertainty, in percent of the value

        Systematic uncertain"""
        self.hiLt = [1,0]
        """heavy ion life time and its uncertainty."""
        self.CSF = [1,0]
        """Charge state fraction of recoils and uncertainty."""
        self.mcpTrans = [0.95**2*0.98**3,0]
        self.mcpEff = [0.8,0.1]
        
    def calc_Ecm(self,massProj,massTarg):
        """Calculate the cm energy.

        \'massProj\' and \'massTarg\' are the projectile and target mass in amu, respectively."""
        self.Ecm[0] = ( (self.eBeam[0] + self.eBeam[1])/2.0 ) * massProj
        self.Ecm[0] = self.Ecm[0]/1000.0 * massTarg/(massTarg+massProj)
        self.Ecm[1] = sqrt(2 * (self.dEnergy**2))/2.0
        self.Ecm[1] = self.Ecm[1]/1000.0 * massTarg/(massTarg+massProj)


############## Set ######################
    
class Set(RunSettings):
    """Set of runs with same incoming beam energy and target pressure."""
    def __init__(self):
        RunSettings.__init__(self)
        self.runNumber = []
        self.runString = ''
        """A more readable list of the runs in the set"""
        self.pTargetList = []
        self.hiLtList = []
        self.geantTrans = [0.97,0.03]
        self.geantBGO = [0.6,0.1]
        self.recoilEvents = [0,0]
        """Observed single and coincidence events."""
        self.reacYield = [[0,0,0],[0,0,0]]
        """Yield for singles and coincidences with uncertainties.
        
        Statistical and systematic uncertainty separately"""
        self.reacSFactor =  [[0,0,0],[0,0,0]]
        """Non-resonant S-factor for singles and coincidences with uncertainties in keVb.
        
        Statistical and systematic uncertainty separately"""
        self.reacSigma =  [[0,0,0],[0,0,0]]
        """Non-resonant cross-setion for singles and coincidences with uncertainties in nbarn.
        
        Statistical and systematic uncertainty separately"""
        self.reacWg =  [[0,0,0],[0,0,0]]
        """OmegaGamma for singles and coincidences with uncertainties."""
        self.RList = [[0,0,0]]
        """Normalisation factors from runs."""
        self.R1List = [[0,0,0]]
        """Normalisation factors from runs using SB1.

        Only used for elastics, not in yield calculation"""
        self.IoverN = [0,0,0]
        """Particles on FC4 over elastics on SB0 at the beginning and end of run with uncertainties.

        Not Rutherford corrected."""
        self.IoverN1 = [0,0,0]
        """Particles on FC4 over elastics on SB1 at the beginning and end of run with uncertainties.

        Not Rutherford corrected."""
        self.IoNList =[[0,0,0]]
        self.IoN1List =[[0,0,0]]
        self.beamInTarget = [0,0,0]
        """Using the beam that made it into the target, i.e. without multiplying with target transmission.
        
        Statistic and systematic uncertainty separately"""
        self.totTrans = [[0,0,0],[0,0,0]]
        """Total separator efficiency singles, coincidences with uncertainties
        
        Statistic and systematic uncertainty separately"""
        self.targetDensity = [0,3]
        """Number density In 1/cm^3

        Second value is uncertainty in % of value.
        This is replaced by the calculated value in calc_target_density"""
        self.stoppingPower = [0,0,0]
        """Stopping power in kev/(ug*cm^2).
        
        Statistical and systematic uncertainty separately"""
        self._upperLimit = [0,0]
        """Upper limit switch: 0 for normal calculation.
        
        1 for cases where the Feldman&Cousins limits will be used,
        2 for cases with 0 events, where only an upper limit (again based on Feldman&Cousins) is used.
        Two values for singles and coincidences"""

    def average_list(self, aList):
        """Returns the average of all entries in aList of values and uncertainties.

        If entries are nan, they are removed from the list"""
        nX = 0
        result = 0
        dResult = 0
        dRes2 = 0
        chi=0
        for x,dx,dx2 in aList[:]:
            #            if (isnan(r) or r==0):
            if (isnan(x)):
                aList.remove([x,dx,dx2])
                continue
            result += x
            dResult += dx**2
            nX += 1
            if (dRes2 == 0):
                dRes2 = dx2/x
        try:
            result /= nX
            dResult = sqrt(dResult)/nX
            dRes2 = result * dRes2
        except ZeroDivisionError:
            result = dResult = dRes2= float('nan')
            return [result,dResult,dRes2]     
        for x,dx,dx2 in aList[:]:
            chi += (x-result)**2/dx**2
        if chi > 1.0:
            dResult *= sqrt(chi)
        return [result,dResult,dRes2]
    
    def calc_averages(self):
        """Calculate average target pressure and HI lt."""
        self.pTarget = 0
        for x in self.pTargetList:
            self.pTarget += x/len(self.pTargetList)
        self.hiLt[0] = 0
        uncert = 0
        for x,dx in self.hiLtList:
            self.hiLt[0] += x/len(self.hiLtList)
            uncert += dx**2
        self.hiLt[1] = sqrt(uncert)/len(self.hiLtList)
        
    def calc_beam_in_target(self):
        """Calculates the beam in the target.

        The target transmission is not corrected for;
        it is assumed that the target transmission is largely determined by the entrance apperture,
        so what doesn't make it through the target doesn't make it into the target and cannot react."""
        self.beamInTarget[0] = self.sb0Tot[0] * (self.eBeam[0]**2)/self.pTarget * self.R[0]
        err1 = ( (self.eBeam[0]**2)/self.pTarget * self.R[0])**2 * self.sb0Tot[1]**2
        #     err2 = (self.sb0Tot[0] * 2 * (self.eBeam[0])/self.pTarget * self.R[0])**2 * self.dEnergy**2
        #err3 = (self.sb0Tot[0] * (self.eBeam[0]**2)/(self.pTarget**2) * self.R[0])**2 * self.dPressure**2
        err4 = (self.sb0Tot[0] * (self.eBeam[0]**2)/(self.pTarget**2))**2 *  self.R[1]**2
        self.beamInTarget[1] = sqrt(err1 + err4)
        systUncert = sqrt( (self.dFc4/100.0)**2 + (self.tarTrans[1]/self.tarTrans[0])**2 + (self.R[2]/self.R[0])**2  + 2*(self.dEnergy/self.eBeam[0])**2 + (self.dPressure/100.0)**2) 
        self.beamInTarget[2] = self.beamInTarget[0]*systUncert

    def calc_beam_in_target_sb1(self):
        """Calculates the beam in the target.

        The target transmission is not corrected for;
        it is assumed that the target transmission is largely determined by the entrance apperture,
        so what doesn't make it through the target doesn't make it into the target and cannot react."""
        self.beamInTarget[0] = self.sb1Tot[0] * (self.eBeam[0]**2)/self.pTarget * self.R1[0]
        err1 = ( (self.eBeam[0]**2)/self.pTarget * self.R1[0])**2 * self.sb1Tot[1]**2
        #     err2 = (self.sb0Tot[0] * 2 * (self.eBeam[0])/self.pTarget * self.R[0])**2 * self.dEnergy**2
        #err3 = (self.sb0Tot[0] * (self.eBeam[0]**2)/(self.pTarget**2) * self.R[0])**2 * self.dPressure**2
        err4 = (self.sb1Tot[0] * (self.eBeam[0]**2)/(self.pTarget**2))**2 *  self.R1[1]**2
        self.beamInTarget[1] = sqrt(err1 + err4)
        systUncert = sqrt( (self.dFc4/100.0)**2 + (self.tarTrans[1]/self.tarTrans[0])**2 + (self.R1[2]/self.R1[0])**2  + 2*(self.dEnergy/self.eBeam[0])**2 + (self.dPressure/100.0)**2) 
        self.beamInTarget[2] = self.beamInTarget[0]*systUncert

    def calc_elastics(self, maxDev=0.4):
        """Calculate elastics.

        Calculates from values for individual runs: normalisation factor R (SB0) and R1 (SB1), FC4/SB0, and FC4/SB1.
        \'maxDev\' is how much an individual value may deviate from the average (in ave * maxDev) before being skipped. """
        self.calc_R(maxDev)
        for aList in [self.R1List, self.IoNList, self.IoN1List]:
            self.remove_zero_from_list(aList)
            self.remove_outliers_from_list(aList,3)
        self.R1 = self.calc_value_from_list(self.R1List,maxDev)
        self.IoverN = self.calc_value_from_list(self.IoNList,maxDev)
        self.IoverN1 = self.calc_value_from_list(self.IoN1List,maxDev)
        
    def calc_R(self,maxDev=0.4):
        """Calculate normalisation factor R from all run Rs.

        \'maxDev\' is how much an individual R may deviate from Rave (in Rave * maxDev) before being skipped."""
        self.remove_zero_from_list(self.RList)
        self.remove_outliers_from_list(self.RList,3)
        self.R = self.calc_value_from_list(self.RList,maxDev)
        
    def calc_R1(self,maxDev=0.4):
        """Calculate normalisation factor R based on SB1 from all run Rs.

        \'maxDev\' is how much an individual R may deviate from Rave (in Rave * maxDev) before being skipped."""
        self.remove_zero_from_list(self.R1List)
        self.remove_outliers_from_list(self.R1List,3)
        self.R1 = self.calc_value_from_list(self.R1List,maxDev)

    def calc_Sfactor(self, projZ, targZ, reducedMass):
        """Calculate the S-factor at the cm energy.

        calc_sigma must have been called first."""
        twoPiEta = 31.29 *projZ*targZ * sqrt(reducedMass/(self.Ecm[0]*1000))
        err1 = - 31.29 * projZ*targZ * sqrt(reducedMass) * 1/2.0 * ( (self.Ecm[0]*1000)**(-1.5)) * self.Ecm[1] * 1000
        for i in range(2):
            self.reacSFactor[i][0] = self.reacSigma[i][0]/1e9 * self.Ecm[0] * 1000 * exp(twoPiEta )
            err2 = ( 1/1e9 * self.Ecm[0] * 1000 * exp(twoPiEta) )**2 * self.reacSigma[i][1]**2
#            err3 = ( self.reacSigma[i][0]/1e9 * 1000 * exp(twoPiEta) )**2 * self.Ecm[1]**2
#            err4 = (self.reacSigma[i][0]/1e9 *1000 * self.Ecm[0] * exp(twoPiEta) )**2 * err1**2
            if isnan(err2): err2 = 0
            self.reacSFactor[i][1] = sqrt( err2  )
            if self.reacSFactor[i][0] ==0:
                systUncert = sqrt( (self.Ecm[1]/self.Ecm[0])**2 + (err1)**2 + (self.reacSigma[i][2]/self.reacSigma[i][1])**2)
                self.reacSFactor[i][2] = self.reacSFactor[i][1] * systUncert
            else:
                systUncert = sqrt( (self.Ecm[1]/self.Ecm[0])**2 + (err1)**2 + (self.reacSigma[i][2]/self.reacSigma[i][0])**2)
                self.reacSFactor[i][2] = self.reacSFactor[i][0] * systUncert

            
#            if self._upperLimit[i] == 2:
#                self.reacSFactor[i][0] += self.reacSFactor[i][1]
#                self.reacSFactor[i][1] = float('nan')
    
    def calc_sigma(self, targetLength=[12.3,0.1]):
        """Calculate non-resonant cross section for singles and coincidences

        \'targetLength\' is in cm"""
        for i in range(2):
            self.reacSigma[i][0] = self.reacYield[i][0]/(self.targetDensity[0] * targetLength[0])*1e9*1e24
            err1 = (1/(self.targetDensity[0] * targetLength[0]))**2 * self.reacYield[i][1]**2
#           err2 = (self.reacYield[i][0]/(self.targetDensity[0]**2 * targetLength[0]))**2 * self.targetDensity[1]**2
#           err3 = (self.reacYield[i][0]/(self.targetDensity[0] * targetLength[0]**2))**2 * targetLength[1]**2
            if isnan(err1): err1 = 0
            self.reacSigma[i][1] = sqrt(err1 ) * 1e9 * 1e24
            if (self.reacYield[i][0] == 0):
                self.reacSigma[i][2] = self.reacSigma[i][1] *sqrt((self.reacYield[i][2]/self.reacYield[i][1])**2 + (targetLength[1]/targetLength[0])**2  + (self.targetDensity[1]/self.targetDensity[0])**2)
            else:
                self.reacSigma[i][2] = self.reacSigma[i][0] * sqrt((self.reacYield[i][2]/self.reacYield[i][0])**2 + (targetLength[1]/targetLength[0])**2  + (self.targetDensity[1]/self.targetDensity[0])**2)                
            # if self._upperLimit[i] == 2:
            #     self.reacSigma[i][0] += self.reacSigma[i][1]
            #     self.reacSigma[i][1] = float('nan')

    def calc_stopping_power(self,projA, targetMassKg, targetLength=[12.3,0.1]):
        """Calculate the stopping power based on beam energy before and after the target.

        Calculates first the target density in ug/cm^2, then the stopping power in keV/(ug/cm^2).
        \'projA\' is the beam mass number, \'targetMassKg\' is the mass in kg, \'targetLength\' is in cm.
        All uncertainties assumed systematic, so stoppingPower[1] is 0, and stoppingPower[2] is the systematic uncertainty."""
        density = self.targetDensity[0] * targetMassKg * 1e9 * targetLength[0]
        err1 = (targetMassKg * 1e9 * targetLength[0])**2 * self.targetDensity[1]**2
        err2 = (self.targetDensity[0] * targetMassKg * 1e9)**2 *  targetLength[1]**2
        dDensity = sqrt(err1 + err2)
#        dDensSys = sqrt( (targetLength[1]/targetLength[0])**2 + (self.targetDensity[1]/self.targetDensity[0])**2)
        try:
            self.stoppingPower[0] = (self.eBeam[0]-self.eBeam[1]) * float(projA)/density
        except ZeroDivisionError:
            print "could not calculate stopping power (target density is 0?)"
            return
        else:
            err1 = ((float(projA)/density)**2 * self.dEnergy**2)
            err2 = ((self.eBeam[0]-self.eBeam[1]) * float(projA)/density**2)**2 * dDensity**2
            self.stoppingPower[1] = 0
#            systUncert = sqrt( (dDensSys/density)**2 + (self.dEnergy/self.eBeam[0])**2 + (self.dEnergy/self.eBeam[1])**2)
#            self.stoppingPower[2] = self.stoppingPower[0] * systUncert
            self.stoppingPower[2] =  sqrt(err1 + err1 + err2)
  
    def calc_target_density(self, nuTarget=1):
        """Calculate target density in 1/cm^3.

        nuTarget is the number of atoms per molecule, i.e. 1 for helium and 2 for hydrogen.
        Uncertainty is systematic."""
        self.targetDensity[0] = 9.66e18 * nuTarget * self.pTarget/300.0
        err1 = (9.66e18 * nuTarget/300.0)**2 * (self.dPressure/100.0*self.pTarget)**2
        self.targetDensity[1] = sqrt(err1)
##        self.targetDensity[1] = self.targetDensity[1]/100.0 * self.targetDensity[0]

    def calc_transmission(self):
        """Calculate the separator efficiency for singles and coincidences.

        \'mcpEff\' is a tuple containing MCP efficiency and uncertainty."""
        self.totTrans[0][0] = self.CSF[0] * self.mcpTrans[0] * self.mcpEff[0] * self.sepTrans[0] * self.hiLt[0] * self.geantTrans[0]
#        err1 = (self.mcpTrans[0] * self.mcpEff[0] * self.sepTrans[0] * self.hiLt[0] * self.geantTrans[0])**2 * self.CSF[1]**2
        err4 = (self.CSF[0] * self.mcpTrans[0] * self.mcpEff[0] * self.sepTrans[0] * self.geantTrans[0])**2 * self.hiLt[1]**2
#        err5 = (self.CSF[0] * self.mcpTrans[0] * self.mcpEff[0] * self.hiLt[0] * self.geantTrans[0])**2 * self.sepTrans[1]**2
#        err6 = (self.CSF[0] * self.mcpTrans[0] * self.mcpEff[0] * self.sepTrans[0] * self.hiLt[0])**2 * self.geantTrans[1]**2
        self.totTrans[0][1] = sqrt( err4 )
        systUncert = sqrt((self.mcpEff[1]/self.mcpEff[0])**2 + (self.mcpTrans[1]/self.mcpTrans[0])**2  + (self.sepTrans[1]/self.sepTrans[0])**2  +(self.geantTrans[1]/self.geantTrans[0])**2  + (self.CSF[1]/self.CSF[0])**2 )
        self.totTrans[0][2] = self.totTrans[0][0] * systUncert
        self.totTrans[1][0] = self.CSF[0] * self.mcpTrans[0] * self.mcpEff[0] * self.sepTrans[0] * self.hiLt[0] * self.geantBGO[0] * self.geantTrans[0]
#        err1 = (self.mcpTrans[0] * self.mcpEff[0] * self.sepTrans[0] * self.hiLt[0] * self.geantBGO[0] * self.geantTrans[0])**2 * self.CSF[1]**2
#        err4 = (self.CSF[0] * self.mcpTrans[0] * self.mcpEff[0] * self.hiLt[0] * self.geantBGO[0] * self.geantTrans[0])**2 * self.sepTrans[1]**2
        err5 = (self.CSF[0] * self.mcpTrans[0] * self.mcpEff[0] * self.sepTrans[0] * self.geantBGO[0] * self.geantTrans[0])**2 * self.hiLt[1]**2
#        err7 = (self.CSF[0] * self.mcpTrans[0] * self.mcpEff[0] * self.sepTrans[0] * self.hiLt[0] * self.geantBGO[0])**2 * self.geantTrans[1]**2
        self.totTrans[1][1] = sqrt(err5)
        systUncert = sqrt((self.mcpEff[1]/self.mcpEff[0])**2 + (self.mcpTrans[1]/self.mcpTrans[0])**2  + (self.sepTrans[1]/self.sepTrans[0])**2  +(self.geantTrans[1]/self.geantTrans[0])**2  + (self.CSF[1]/self.CSF[0])**2  + (self.geantBGO[1]/self.geantBGO[0])**2)
        self.totTrans[1][2] = self.totTrans[1][0] * systUncert

    def calc_value_from_list(self,list,maxDev=0.4):
       """Calculate a value with uncertainty from a list of values and uncertainties

       \'maxDev\' is how much an individual value may deviate from the average (in ave * maxDev) before being skipped. """ 
       ## if self.runString == '':
       ##     self.runString = self.runs_to_string()
       ## print self.runString
       ## print list
       ## print "------------------------------------------------------------"
       result = [0,0,0]
       result = self.average_list(list);
       for x,dx,dx2 in list[:]:
           if abs(x-result[0])>result[0]*maxDev: list.remove([x,dx,dx2])
       result = self.average_list(list)
       ## print list
       ## print "*****************" + str(result) + "**************************"
       return result;
                 
    def calc_wg(self, projA, targetA, stoppingFactor):
        """Calculate resonant omega-gamma.

        \'projA\' and \'targetA\' are the projectile and target mass numbers.
        \'stoppingFactor\' is the conversion factor for the stopping power from  keV/(ug/cm^2) to eV/(10^15/cm^2)."""
        if stoppingFactor == 0: return 1
        reducedMass = (projA*targetA)/float(projA+targetA)
        factor = 1/(4125.5e-9) * reducedMass * targetA/float(projA+targetA)
        for i in range(2):
            self.reacWg[i][0] = factor * self.Ecm[0]*1000 * self.stoppingPower[0] * stoppingFactor * self.reacYield[i][0]
#            err1 = (factor * self.stoppingPower[0] * stoppingFactor * self.reacYield[i][0])**2 * (self.Ecm[1]*1000)**2
            err2 = (factor * self.Ecm[0]*1000 * self.reacYield[i][0])**2 * (stoppingFactor * self.stoppingPower[1])**2 
            err3 = (factor * self.Ecm[0]*1000 * self.stoppingPower[0] * stoppingFactor)**2 * self.reacYield[i][1]**2
            if (self._upperLimit[i] == 2 and isnan(err3)): err3 = 0
            self.reacWg[i][1] = sqrt(err2 + err3)
            if self._upperLimit[i] == 2:
                systUncert = sqrt( (self.Ecm[1]/self.Ecm[0])**2 + (self.reacYield[i][2]/self.reacYield[i][1])**2 + (self.stoppingPower[2]/self.stoppingPower[0])**2)
                self.reacWg[i][2] = self.reacWg[i][1] * systUncert
            else:
                systUncert = sqrt( (self.Ecm[1]/self.Ecm[0])**2 + (self.reacYield[i][2]/self.reacYield[i][0])**2 + (self.stoppingPower[2]/self.stoppingPower[0])**2)
                self.reacWg[i][2] = self.reacWg[i][0] * systUncert
        return 0

    def calc_yield(self):
        """Calculate yield for singles and coincidences"""
        for i in range(2):
            recoilL = []
            value = []
            uncert = []
            if self._upperLimit[i] == 2:
                recoilL.append([lowStatLimits[0][1], 0])
            elif self._upperLimit[i] == 1:
                for j in range(2):
                    recoilL.append([lowStatLimits[self.recoilEvents[i]][j], 0])
            else: recoilL.append([self.recoilEvents[i],self.recoilEvents[i]])
            for [recoils, dRecoils]  in recoilL:
                value.append(recoils/(self.totTrans[i][0]*self.beamInTarget[0]))
                err1 = (1/(self.totTrans[i][0]*self.beamInTarget[0]))**2 * dRecoils
#                err2 = (recoils/((self.totTrans[i][0]**2) * self.beamInTarget[0]))**2 * self.totTrans[i][1]**2 
#                err3 = (recoils/(self.totTrans[i][0]*self.beamInTarget[0]**2))**2 * self.beamInTarget[1]**2
                uncert.append(sqrt(err1 ))
            if self._upperLimit[i] == 2:                   
                self.reacYield[i][0] = 0
                self.reacYield[i][1] = sqrt( value[0]**2 + uncert[0]**2 )
#                self.reacYield[i][1] = float('nan')
            elif self._upperLimit[i] == 1:
                limit = []
                for j in range(len(value)):
                    limit.append(value[j] + (-1)**(j+1) * uncert[j])
                self.reacYield[i][0] = (limit[0] + limit[1])/2
                self.reacYield[i][1] = (limit[1] - limit[0])/2
            else:
                self.reacYield[i][0] = value[0]
                self.reacYield[i][1] = uncert[0]
            uncertSys = sqrt( (self.totTrans[i][2]/self.totTrans[i][0])**2 + (self.beamInTarget[2]/self.beamInTarget[0])**2 )
            if self.reacYield[i][0] ==0:
                self.reacYield[i][2] = uncertSys * self.reacYield[i][1]
            else:
                self.reacYield[i][2] = self.reacYield[i][0]*uncertSys
            

    def get_geant_data(self,line):
        """Get GEANT data from \'line\'.
        
        Format:
        E_in   pTar    transm  +/-    BGOeff  +/- """
        dataList = line.split()
        self.geantTrans = [float(dataList[2]),float(dataList[3])]
        self.geantBGO = [float(dataList[4]),float(dataList[5])]
        
    def get_recoil_data(self,line,upperLimitSwitch=True):
        """Get recoil data from \'line\'.

        Adds events to recoilEvents.
        Format
        Ein   pTar   singles    coincidences
        \'upperLimitSwitch\' determines whether (True) of not (False) the limits by Feldman & Cousins will be used for low numbers"""
        dataList = line.split()
        self.recoilEvents[0] += int(dataList[2])
        self.recoilEvents[1] += int(dataList[3])
        if upperLimitSwitch:
            for i,events in enumerate(self.recoilEvents):
                if events == 0: self._upperLimit[i] = 2
                elif events < len(lowStatLimits): self._upperLimit[i] = 1
                else: self._upperLimit[i] = 0
        
    def get_normalisation(self, format=0):
        """Returns a string representation of the set

        If format=1 the titles are returned, used as the first line in output files."""
        if format == 1:
            return "#run\tEin\tEout\tEcm\t+/-\tpressure\tTargDensity\t+/-\tstopping\t+/-stat\t+/-sys\ttarTrans\t+/-\tsepTrans\t+/-\thiLt\t+/-\tCSF\t+/-\tMCPeff\t+/-\tgeantTransm\t+/-\tgeantBGOeff\t+/-\tSB0tot\t+/-\tSB1tot\t+/-\tR\t+/-stat\t+/-sys\tbeam\t+/-stat\t+/-sys\tR1\t+/-stat\t+/-sys\tFC4/SB0\t+/-stat\t+/-sys\tFC4/SB1\t+/-stat\t+/-sys\tsingles\tcoincidences\n"
        if self.runString == '':
            self.runString = self.runs_to_string()
        result = self.runString + "\t"
        result += "\t".join(map(str, self.eBeam)) + "\t"
        result += "\t".join(map(str, self.Ecm)) + "\t"
        result += str(self.pTarget) + "\t"
        result += "\t".join(map(str, self.targetDensity)) + "\t"
        result += "\t".join(map(str, self.stoppingPower)) + "\t"
        result += "\t".join(map(str, self.tarTrans)) + "\t"
        result += "\t".join(map(str, self.sepTrans)) + "\t"
        result += "\t".join(map(str, self.hiLt)) + "\t"
        result += "\t".join(map(str, self.CSF)) + "\t"
        result += "\t".join(map(str, self.mcpEff)) + "\t"
        result += "\t".join(map(str, self.geantTrans)) + "\t"
        result += "\t".join(map(str, self.geantBGO)) + "\t"
        result += "\t".join(map(str, self.sb0Tot)) + "\t"
        result += "\t".join(map(str, self.sb1Tot)) + "\t"
        result += "\t".join(map(str, self.R)) + "\t"       
        result += "\t".join(map(str, self.beamInTarget)) + "\t"
        result += "\t".join(map(str, self.R1)) + "\t"       
        result += "\t".join(map(str, self.IoverN)) + "\t"       
        result += "\t".join(map(str, self.IoverN1)) + "\t"           
        result += "\t".join(map(str, self.recoilEvents)) + "\n"           
        return result
        
    def make_set_from_run(self,run):
        """Create a new set taking the relevant data from run. """
        self.runNumber.append(run.runNumber)
        self.Ecm = run.Ecm
        self.eBeam = run.eBeam
        self.pTargetList.append(run.pTarget)
        self.pTarget = run.pTarget
        self.tarTrans = run.tarTrans
        self.sepTrans = run.sepTrans
        self.sb0Tot[0] += run.sb0Tot[0]
        self.sb0Tot[1] = sqrt(self.sb0Tot[0])
        self.sb1Tot[0] += run.sb1Tot[0]
        self.sb1Tot[1] = sqrt(self.sb1Tot[0])
        self.hiLtList.append(run.hiLt)
        self.CSF = run.CSF
        self.RList.append(run.R[0])
        self.RList.append(run.R[1])
        self.R1List.append(run.R1[0])
        self.R1List.append(run.R1[1])
        self.IoNList.append(run.IoverN[0])
        self.IoNList.append(run.IoverN[1])
        self.IoN1List.append(run.IoverN1[0])
        self.IoN1List.append(run.IoverN1[1])
        return self

    def make_set_from_set_step2(self,line,upperLimitSwitch=True):
        """Get set data from line.

        This function is used to read in a file with set data written by \'get_normalisation\'. This way, the normalisation can be done up to that point, the output file edited, and then read back  in to continue with the modified sets. 
        Format:
        Runs   Ein   Eout   Ecm  +/-  pressure  TargDensity  +/-   stopping  +/-  tarTrans  +/-   sepTrans   +/-   HI-lt   +/-   CSF   +/-   MCPeff   +/-   geantTransm   +/-   geantBGOeff   +/-   SB0tot   +/-   SB1tot   +/-   R   +/-   beam   +/-  (   R1   +/-     FC4/SB0   +/-     FC4/SB1   +/- ( singles    coincidences)) """
        lList = line.split()
        self.runString = lList[0]
        self.eBeam = [float(lList[1]),float(lList[2])]
        self.Ecm = [float(lList[3]),float(lList[4])]
        self.pTarget = float(lList[5])
        self.targetDensity = [float(lList[6]),float(lList[7])]
        self.stoppingPower = [float(lList[8]),float(lList[9])]
        self.tarTrans = [float(lList[10]),float(lList[11])]
        self.sepTrans = [float(lList[12]),float(lList[13])]
        self.hiLt = [float(lList[14]),float(lList[15])]
        self.CSF = [float(lList[16]),float(lList[17])]
        self.mcpEff = [float(lList[18]),float(lList[19])]
        self.geantTrans = [float(lList[20]),float(lList[21])]
        self.geantBGO = [float(lList[22]),float(lList[23])]
        self.sb0Tot = [float(lList[24]),float(lList[25])]
        self.sb1Tot = [float(lList[26]),float(lList[27])]
        self.R = [float(lList[28]),float(lList[29])]
        self.beamInTarget = [float(lList[30]),float(lList[31])]
        if (len(lList)==38 or len(lList)==40):
            self.R1 =  [float(lList[32]),float(lList[33])]
            self.IoverN =  [float(lList[34]),float(lList[35])]
            self.IoverN1 =  [float(lList[36]),float(lList[37])]
            if (len(lList)==40):
                self.recoilEvents =  [int(lList[38]),int(lList[39])]
                if upperLimitSwitch:
                    for i,events in enumerate(self.recoilEvents):
                        if events == 0: self._upperLimit[i] = 2
                        elif events < len(lowStatLimits): self._upperLimit[i] = 1
                        else: self._upperLimit[i] = 0
        return self
    
    def remove_outliers_from_list(self,list,limit):
        """Removes entries from a list [[0,0],[x,dx]] that are more than limit x other value.

        Requires list of at least 3 elements."""
        if (len(list)>2):
            for i in range(len(list)):
                if (i==0):
                    if (list[i][0]> limit * list[i+1][0] and list[i][0]> limit * list[i+2][0]):
                        list.pop(i)
                        break
                elif i<len(list)-1:
                    if (list[i][0]> limit * list[i+1][0] and list[i][0]> limit * list[i-2][0]):
                        list.pop(i)
                        break
                elif i==len(list)-1:
                    if (list[i][0]> limit * list[i-1][0] and list[i][0]> limit * list[i-2][0]):
                        list.pop(i)
                        break
            else:
                return
            self.remove_outliers_from_list(list,limit)
                    
        
    def remove_zero_from_list(self,list):
        """Removes entries that are 0 from a list [[0,0],[x,dx]]"""
        for x,dx,dx2 in list[:]:
            if (isnan(x) or x==0):
                list.remove([x,dx,dx2])
        
    def runs_to_string(self):
        """Returns a more readable string of the run numbers."""
        result = ""
        i = 0
        if len(self.runNumber) == 1: return str(self.runNumber[0])
        while i < len(self.runNumber)-1:
            sep = ","
            j = i+1
            check = 0
            while j < len(self.runNumber):
                if self.runNumber[j]-self.runNumber[i] == j-i:
                    sep = "-"
                    j += 1
                    check = 1
                else: break
            j = j - check
            result +=  str(self.runNumber[i]) +sep + str(self.runNumber[j])
            i = j + 1
            if i < len(self.runNumber)-1: result += "," 
        return result
    
    def summarise_results(self,format=0):
        """Returns a string representation of the results for the set.

        If format=1 the titles are returned, used as the first line in output files."""
        if format == 1:
            return "#run\tEin\tEout\tEcm\t+/-\tpressure\tR\t+/-stat\t+/-sys\tbeam\t+/-stat\t+/-sys\ttransm(single)\t+/-stat\t+/-sys\ttransm(coinc)\t+/-stat\t+/-sys\tsingles\tcoincidences\tyield(single)\t+/-stat\t+/-sys\tyield(coinc)\t+/-stat\t+/-sys\tsigma(single)\t+/-stat\t+/-sys\tsigma(coinc)\t+/-stat\t+/-sys\twg(singles)\t+/-stat\t+/-sys\twg(coinc)\t+/-stat\t+/-sys\tSfactor(singles)\t+/-stat\t+/-sys\tSfactor(coinc)\t+/-stat\t+/-sys\n"
        if self.runString == '':
            self.runString = self.runs_to_string()
        result = self.runString + "\t"
        result += "\t".join(map(str, self.eBeam)) + "\t"
        result += "\t".join(map(str, self.Ecm)) + "\t"
        result += str(self.pTarget) + "\t"
        result += "\t".join(map(str, self.R)) + "\t"       
        result += "\t".join(map(str, self.beamInTarget)) + "\t"
        result += "\t".join(map(str, self.totTrans[0])) + "\t"       
        result += "\t".join(map(str, self.totTrans[1])) + "\t"       
        result += "\t".join(map(str, self.recoilEvents)) + "\t"       
        result += "\t".join(map(str, self.reacYield[0])) + "\t"
        result += "\t".join(map(str, self.reacYield[1])) + "\t"
        result += "\t".join(map(str, self.reacSigma[0])) + "\t"
        result += "\t".join(map(str, self.reacSigma[1])) + "\t"
        result += "\t".join(map(str, self.reacWg[0])) + "\t"       
        result += "\t".join(map(str, self.reacWg[1])) + "\t"       
        result += "\t".join(map(str, self.reacSFactor[0])) + "\t"       
        result += "\t".join(map(str, self.reacSFactor[1])) + "\n"       
        return result
        
                

            

############# Run  ####################

class Run(RunSettings):
    """Contains data for a single run."""
    def __init__(self):
        RunSettings.__init__(self)
        self.runNumber = 0
        self.fc4 = [0,0]
        """Cup readings before and after the run."""
        self.dFc4 = 3
        """FC4 reading uncertainty, in percent of the value"""
        self.sb0Start = 0
        """SB0/s at the beginning of the run."""
        self.sb0End = 0
        """SB0/s at the end of the run."""
        self.sb1Start = 0
        """SB1/s at the beginning of the run."""
        self.sb1End = 0
        """SB1/s at the end of the run."""
        self.R = [[0,0,0],[0,0,0]]
        """Normalisation factor R (Rutherford corrected) at the beginning and end of run with uncertainties."""
        self.R1 = [[0,0,0],[0,0,0]]
        """Normalisation factor R (Rutherford corrected) at the beginning and end of run with uncertainties.

        Calculated using SB1, this is just for comparison, R (for SB0) is used for the yield calculations."""
        self.IoverN = [[0,0,0],[0,0,0]]
        """Particles on FC4 over elastics on SB0 at the beginning and end of run with uncertainties.

        Not Rutherford corrected."""
        self.IoverN1 = [[0,0,0],[0,0,0]]
        """Particles on FC4 over elastics on SB1 at the beginning and end of run with uncertainties.

        Not Rutherford corrected."""

        
    def calc_IoverN(self):
        """Calculates the individual FC4 current per SB0 for the beginning and the end of the run."""
        try:
            self.IoverN[0][0] = self.fc4[0]/self.sb0Start[0] * self.tarTrans[0]
            err1= err2=err3 = 0
            err1 = (1./self.sb0Start[0] * self.tarTrans[0])**2 * (self.fc4[0]*self.dFc4/100.0)**2
            err2 = (self.fc4[0]/(self.sb0Start[0]**2) * self.tarTrans[0])**2 * self.sb0Start[1]**2
            err3 = (self.fc4[0]/self.sb0Start[0] )**2 * self.tarTrans[1]**2
            self.IoverN[0][1] = sqrt(err1 + err2 + err3) 
        except ZeroDivisionError:
            self.IoverN[0][0] = float('nan')
            self.IoverN[0][1] = float('nan')
        try:
            err1 = err2 = err3 = err4 = 0
            self.IoverN[1][0] = self.fc4[1]/self.sb0End[0] * self.tarTrans[0]
            err1 = (1./self.sb0End[0] * self.tarTrans[0])**2 * (self.fc4[1]*self.dFc4/100.0)**2
            err2 = (self.fc4[1]/(self.sb0End[0]**2) * self.tarTrans[0])**2 * self.sb0End[1]**2
            err3 = (self.fc4[1]/self.sb0End[0] )**2 * self.tarTrans[1]**2
            self.IoverN[1][1] = sqrt(err1 + err2 + err3) 
        except ZeroDivisionError:
            self.IoverN[1][0] = float('nan')
            self.IoverN[1][1] = float('nan')
            
    def calc_IoverN1(self):
        """Calculates the individual FC4 current per SB1 for the beginning and the end of the run."""
        try:
            self.IoverN1[0][0] = self.fc4[0]/self.sb1Start[0] *self.tarTrans[0]
            err1= err2=err3 = 0
            err1 = (1./self.sb1Start[0] * self.tarTrans[0])**2 * (self.fc4[0]*self.dFc4/100.0)**2
            err2 = (self.fc4[0]/(self.sb1Start[0]**2) *self.tarTrans[0])**2 * self.sb1Start[1]**2
            err3 = (self.fc4[0]/self.sb1Start[0])  * self.tarTrans[1]**2
            self.IoverN1[0][1] = sqrt(err1 + err2 + err3) 
        except ZeroDivisionError:
            self.IoverN1[0][0] = float('nan')
            self.IoverN1[0][1] = float('nan')
        try:
            err1 = err2 = err3 = 0
            self.IoverN1[1][0] = self.fc4[1]/self.sb1End[0] *self.tarTrans[0]
            err1 = (1./self.sb1End[0] *self.tarTrans[0])**2 * (self.fc4[1]*self.dFc4/100.0)**2
            err2 = (self.fc4[1]/(self.sb1End[0]**2) *self.tarTrans[0])**2 * self.sb1End[1]**2
            err3 = (self.fc4[1]/self.sb1End[0] )**2 * self.tarTrans[1]**2
            self.IoverN1[1][1] = sqrt(err1 + err2 + err3) 
        except ZeroDivisionError:
            self.IoverN1[1][0] = float('nan')
            self.IoverN1[1][1] = float('nan')
            
    def calc_R(self):
        """Calculates the individual normalisation factor for the beginning and the end of the run.

        Does not contain FC4 uncertainty (systematic!)"""
        try:
            ## print "run", self.runNumber
            ## print "fc4" , self.fc4[0]
            ## print "hiLt" , self.hiLt[0]
            ## print "sb0Start" , self.sb0Start[0]
            ## print "tarTrans" , self.tarTrans[0]
            ## print "eBeam" , self.eBeam[0]
            self.R[0][0] = self.fc4[0]/self.sb0Start[0] *self.tarTrans[0] * self.pTarget/self.eBeam[0]**2
            err1= err2=err3=err4= 0
            err1 = (self.fc4[0]/(self.sb0Start[0]**2)*self.tarTrans[0] * self.pTarget/self.eBeam[0]**2)**2 * self.sb0Start[1]**2
#            err2 = (self.fc4[0]/self.sb0Start[0]  * self.pTarget/self.eBeam[0]**2)**2 * self.tarTrans[1]**2
#            err3 = (self.fc4[0]/self.sb0Start[0] *self.tarTrans[0] * 1/self.eBeam[0]**2)**2 * self.dPressure**2
#            err4 = (self.fc4[0]/self.sb0Start[0] *self.tarTrans[0] * 2 * self.pTarget/(self.eBeam[0]**3))**2 * self.dEnergy**2
            self.R[0][1] = sqrt(err1 ) 
        except ZeroDivisionError:
            self.R[0][0] = float('nan')
            self.R[0][1] = float('nan')
            self.R[0][2] = float('nan')
        try:
            err1 = err2 = err3 = err4 = err5 = 0
            self.R[1][0] = self.fc4[1]/self.sb0End[0] *self.tarTrans[0] * self.pTarget/self.eBeam[0]**2
#            err1 = (1./self.sb0End[0] *self.tarTrans[0] * self.pTarget/self.eBeam[0]**2)**2 * (self.fc4[1]*self.dFc4/100.0)**2
            err2 = (self.fc4[1]/(self.sb0End[0]**2) *self.tarTrans[0] * self.pTarget/self.eBeam[0]**2)**2 * self.sb0End[1]**2
#            err3 = (self.fc4[1]/self.sb0End[0]  * self.pTarget/self.eBeam[0]**2)**2 * self.tarTrans[1]**2
#            err4 = (self.fc4[1]/self.sb0End[0] *self.tarTrans[0] * 1/self.eBeam[0]**2)**2 * self.dPressure**2
#            err5 = (self.fc4[1]/self.sb0End[0] *self.tarTrans[0] * 2* self.pTarget/(self.eBeam[0]**3))**2 * self.dEnergy**2
            self.R[1][1] = sqrt(err2 ) 
        except ZeroDivisionError:
            self.R[1][0] = float('nan')
            self.R[1][1] = float('nan')
        systUncert = sqrt( (self.tarTrans[1]/self.tarTrans[0])**2 + (self.dPressure/100.0)**2 + (self.dEnergy/self.eBeam[0])**2 + (self.dFc4/100.0)**2)
        self.R[0][2] = self.R[0][0] * systUncert 
        self.R[1][2] = self.R[1][0] * systUncert 

            
    def calc_R1(self):
        """Calculates the individual normalisation factor based on SB1 for the beginning and the end of the run."""
        try:
            self.R1[0][0] = self.fc4[0]/self.sb1Start[0] *self.tarTrans[0] * self.pTarget/self.eBeam[0]**2
            err1= err2=err3=err4=err5 = 0
            err1 = (1./self.sb1Start[0] *self.tarTrans[0] * self.pTarget/self.eBeam[0]**2)**2 * (self.fc4[0]*self.dFc4/100.0)**2
            err2 = (self.fc4[0]/(self.sb1Start[0]**2) *self.tarTrans[0] * self.pTarget/self.eBeam[0]**2)**2 * self.sb1Start[1]**2
            err3 = (self.fc4[0]/self.sb1Start[0]  * self.pTarget/self.eBeam[0]**2)**2 * self.tarTrans[1]**2
            err4 = (self.fc4[0]/self.sb1Start[0] *self.tarTrans[0] * 1/self.eBeam[0]**2)**2 * (self.dPressure/100.0*self.pTarget)**2
            err5 = (self.fc4[0]/self.sb1Start[0] *self.tarTrans[0] * 2 * self.pTarget/(self.eBeam[0]**3))**2 * self.dEnergy**2
            self.R1[0][1] = sqrt(err1 + err2 + err3 + err4 + err5) 
        except ZeroDivisionError:
            self.R1[0][0] = float('nan')
            self.R1[0][1] = float('nan')
        try:
            err1 = err2 = err3 = err4 = err5 = 0
            self.R1[1][0] = self.fc4[1]/self.sb1End[0] *self.tarTrans[0] * self.pTarget/self.eBeam[0]**2
            err1 = (1./self.sb1End[0] *self.tarTrans[0] * self.pTarget/self.eBeam[0]**2)**2 * (self.fc4[1]*self.dFc4/100.0)**2
            err2 = (self.fc4[1]/(self.sb1End[0]**2) *self.tarTrans[0] * self.pTarget/self.eBeam[0]**2)**2 * self.sb1End[1]**2
            err3 = (self.fc4[1]/self.sb1End[0]  * self.pTarget/self.eBeam[0]**2)**2 * self.tarTrans[1]**2
            err4 = (self.fc4[1]/self.sb1End[0] *self.tarTrans[0] * 1/self.eBeam[0]**2)**2 * (self.dPressure/100.0*self.pTarget)**2
            err5 = (self.fc4[1]/self.sb1End[0] *self.tarTrans[0] * 2* self.pTarget/(self.eBeam[0]**3))**2 * self.dEnergy**2
            self.R1[1][1] = sqrt(err1 + err2 + err3 + err4 + err5) 
        except ZeroDivisionError:
            self.R1[1][0] = float('nan')
            self.R1[1][1] = float('nan')
            
    def convert_fc4(self,chargeState):
        """Convert cup readings from eA to particles/s.

        Requires charge state of incoming beam."""
        for i in range(len(self.fc4)):
            if self.fc4[i] > 1e-6:
                self.fc4[i] = float('nan')
            self.fc4[i] = self.fc4[i]/(chargeState*1.602176E-019)
            
    def get_normalisation(self,format=0):
        """Return string of run data

        If format=1 the titles are returned, used as the first line in output files."""
        if format == 1:
            return "#run\tEin\tEout\tEcm\t+/-\tpressure\tFC4before\t+/-\tFC4after\t+/-\ttarTrans\t+/-\tsepTrans\t+/-\tHI_lt\t+/-\tCSF\t+/-\tSB0tot\t+/-\tSB0start\t+/-\tSB0end\t+/-\tRstart\t+/-stat\t+/-sys\tRend\t+/-stat\t+/-sys\tFC4/SB0start\t+/-stat\t+/-sys\tFC4/SB0end\t+/-stat\t+/-sys\tSB1tot\t+/-\tSB1start\t+/-\tSB1end\t+/-\tR1start\t+/-stat\t+/-sys\tR1end\t+/-stat\t+/-sys\tFC4/SB1start\t+/-stat\t+/-sys\tFC4/SB1end\t+/-stat\t+/-sys\n"
        result = str(self.runNumber) + "\t"
        result += "\t".join(map(str, self.eBeam)) + "\t"
        result += "\t".join(map(str, self.Ecm)) + "\t"
        result +=  str(self.pTarget) + "\t"
        #result += "\t".join(map(str, self.fc4)) + "\t"
        result += str(self.fc4[0]) + "\t" + str(self.fc4[0]*self.dFc4/100.0) + "\t" + str(self.fc4[1]) + "\t" + str(self.fc4[1]*self.dFc4/100.0) + "\t"
        result += "\t".join(map(str, self.tarTrans)) + "\t"
        result += "\t".join(map(str, self.sepTrans)) + "\t"
        result += "\t".join(map(str, self.hiLt)) + "\t"
        result += "\t".join(map(str, self.CSF)) + "\t"
        result += "\t".join(map(str, self.sb0Tot)) + "\t"
        result += "\t".join(map(str, self.sb0Start)) + "\t"
        result += "\t".join(map(str, self.sb0End)) + "\t"
        result += "\t".join(map(str, self.R[0])) + "\t"
        result += "\t".join(map(str, self.R[1])) + "\t"
        result += "\t".join(map(str, self.IoverN[0])) + "\t"
        result += "\t".join(map(str, self.IoverN[1])) + "\t"
        result += "\t".join(map(str, self.sb1Tot)) + "\t"
        result += "\t".join(map(str, self.sb1Start)) + "\t"
        result += "\t".join(map(str, self.sb1End)) + "\t"
        result += "\t".join(map(str, self.R1[0])) + "\t"
        result += "\t".join(map(str, self.R1[1])) + "\t"        
        result += "\t".join(map(str, self.IoverN1[0])) + "\t"
        result += "\t".join(map(str, self.IoverN1[1])) + "\n"
        return result

    def get_sb_data(self,line):
        """Get SB and lt data from line.

        Format:
        Run    SB0   +/-   SB1   +/-   lt_HI  +/-   lt_BGO   +/-    SB0start/s   +/-     SB0end/s    +/-    SB1start/s   +/-     SB1end/s    +/- """
        dataList = line.split()
        self.sb0Tot = [float(dataList[1]),float(dataList[2])]
        self.sb1Tot = [float(dataList[3]),float(dataList[4])]
        self.hiLt = [float(dataList[5]),float(dataList[6])]
        self.sb0Start = [float(dataList[9]),float(dataList[10])]
        if self.sb0Start[0]==0:
            self.sb0Start = [float('nan'),float('nan')]
        self.sb0End = [float(dataList[11]),float(dataList[12])]
        if self.sb0End[0]==0:
            self.sb0End = [float('nan'),float('nan')]
        self.sb1Start = [float(dataList[13]),float(dataList[14])]
        if self.sb1Start[0]==0:
            self.sb1Start = [float('nan'),float('nan')]
        self.sb1End = [float(dataList[15]),float(dataList[16])]
        if self.sb1End[0]==0:
            self.sb1End = [float('nan'),float('nan')]

    def make_input_from_run(self):
        result = str(self.runNumber) + "\t"
        result += "\t".join(map(str, self.eBeam)) + "\t"
        result +=  str(self.pTarget) + "\t"
        result += "\t".join(map(str, self.fc4)) + "\t"
        result += "\t".join(map(str, self.tarTrans)) + "\t"
        result += "\t".join(map(str, self.sepTrans)) + "\t"
        result += "\t".join(map(str, self.CSF)) + "\t"
        return result

    def make_run_from_file_step1(self,line):
        """Get run data from line.

        This function is used to read in a file with pre-calculated run data written by \'get_normalisation\'. This way, the normalisation can be done up to that point, the output file edited, and then read back  in to continue with the modified runs. 
        Format:
        Run   Ein   Eout   Ecm   +/-   pressure   FC4before   +/-   FC4after   +/-    tarTrans   +/-   sepTrans   +/-   HI_lt   +/-   CSF   +/-   SB0tot   +/-   SB0start   +/-   SB0end   +/-   Rstart   +/-   Rend   +/-   FC4/SB0start   +/-   FC4/SB0end   +/-   SB1tot   +/-   SB1start   +/-   SB1end   +/-   R1start   +/-   R1end   +/-   FC4/SB1start   +/-   FC4/SB1end   +/-"""
        lList = line.split()
        self.runNumber = int(lList[0])
        self.eBeam = [float(lList[1]),float(lList[2])]
        self.Ecm = [float(lList[3]),float(lList[4])]
        self.pTarget = float(lList[5])
        self.fc4 = [float(lList[6]),float(lList[8])]
        self.tarTrans = [float(lList[10]),float(lList[11])]
        self.sepTrans = [float(lList[12]),float(lList[13])]
        self.hiLt = [float(lList[14]),float(lList[15])]
        self.CSF = [float(lList[16]),float(lList[17])]
        self.sb0Tot = [float(lList[18]),float(lList[19])]
        self.sb0Start =[float(lList[20]),float(lList[21])]
        self.sb0End =[float(lList[22]),float(lList[23])]
        self.R = [[float(lList[24]),float(lList[25])],[float(lList[26]),float(lList[27])]]
        self.IoverN = [[float(lList[28]),float(lList[29])],[float(lList[30]),float(lList[31])]]
        self.sb1Tot = [float(lList[32]),float(lList[33])]
        self.sb1Start =[float(lList[34]),float(lList[35])]
        self.sb1End =[float(lList[36]),float(lList[37])]
        self.R1 = [[float(lList[38]),float(lList[39])],[float(lList[40]),float(lList[41])]]
        self.IoverN1 = [[float(lList[42]),float(lList[43])],[float(lList[44]),float(lList[45])]]
        return self
    
    def make_run_from_input(self,line):
        """Get run data from line. 

        Format:
        run    Ein[keV/u]     Eout[keV/u]    pressure      FC4before [eA]    FC4after [eA]    target  +/-  separator    +/-    CSF     +/-.
        Returns new Run object"""
        dataList = line.split()
        self.runNumber = int(dataList[0])
        self.eBeam = [float(dataList[1]),float(dataList[2])]
        self.pTarget = float(dataList[3])
        self.fc4 = [float(dataList[4]),float(dataList[5])]
        self.tarTrans = [float(dataList[6]),float(dataList[7])]
        self.sepTrans = [float(dataList[8]),float(dataList[9])]
        self.CSF = [float(dataList[10]),float(dataList[11])]
        return self
    

        
    
            
