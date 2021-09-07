"""
# =============================================================================
# P-POP PHOTOMETRY
# A photometry tool for P-POP
# =============================================================================
"""


# =============================================================================
# IMPORTS
# =============================================================================

import numpy as np

import SystemReader


# =============================================================================
# PHOTOMETRYCOMPUTER
# =============================================================================

class PhotometryComputer():
    
    def __init__(self,
                 PathPlanetTable,
                 Filters,
                 Sstar,
                 Splanet,
                 Unit,
                 Mission,
                 SummaryPlots,
                 FigDir,
                 block):
        """
        Parameters
        ----------
        PathPlanetTable: str
            Path of the planet table to be read.
        Filters: list
            List of instances of class Filter.
        Sstar: list
            List of modules of type Photometry for computing the host star
            signal.
        Splanet: list
            List of modules of type Photometry for computing the planet
            signal.
        Unit: 'uJy', 'ph'
            Unit in which the photometry should be computed.
        Mission: 'MIR', 'VIS'
            Wavelength range in which the mission is operating.
        SummaryPlots: bool
            If True, makes summary plots after importing a module.
        FigDir: str
            Directory to which summary plots are saved.
        block: bool
            If True, blocks plots when showing.
        """
        
        # Print.
        print('--> Initializing PhotometryComputer')
        
        self.PathPlanetTable = PathPlanetTable
        self.SysRdr = SystemReader.SystemReader(self.PathPlanetTable)
        self.SysRdr.Open()
        
        self.Filters = Filters
        self.Nfilters = len(self.Filters)
        
        self.Sstar = []
        self.Nsstar = len(Sstar)
        for i in range(self.Nsstar):
            self.Sstar += [Sstar[i].Photometry()]
            if (SummaryPlots == True):
                self.Sstar[-1].SummaryPlots(FigDir=FigDir,
                                            block=block)
        
        self.Splanet = []
        self.Nsplanet = len(Splanet)
        for i in range(self.Nsplanet):
            self.Splanet += [Splanet[i].Photometry()]
            if (SummaryPlots == True):
                self.Splanet[-1].SummaryPlots(FigDir=FigDir,
                                              block=block)
        
        if (Unit == 'uJy'):
            self.Unit = Unit
        elif (Unit == 'ph'):
            self.Unit = Unit
        else:
            print('--> WARNING: '+str(Unit)+' is an unknown unit')
            self.Unit = 'uJy'
        print('--> Using unit '+str(self.Unit))
        
        if (Mission == 'MIR'):
            self.Mission = Mission
        elif (Mission == 'VIS'):
            self.Mission = Mission
        else:
            print('--> WARNING: '+str(Mission)+' is an unknown mission')
            self.Mission = 'MIR'
        print('--> Using mission '+str(self.Mission))
        
        pass
    
    def Run(self):
        """
        """
        
        # Go through all filters.
        for i in range(self.Nfilters):
            
            print('--> Filter %.0f of %.0f: ' % (i+1, self.Nfilters)+self.Filters[i].Name)
            
            # Reset the table flag and the line counter.
            self.TableFlag = False
            self.SysRdr.Reset()
            
            # Get the first system. Then compute the signal of the host star
            # and the planet until the end of the planet population table is
            # reached.
            Sys = self.SysRdr.nextSystem()
            while (Sys is not None):
                
                # Compute the signal of the host star.
                Fstar = []
                for j in range(self.Nsstar):
                    Fstar += [self.Sstar[j].Compute(self.Filters[i],
                                                    Sys,
                                                    self.Unit,
                                                    self.Mission)]
                Fstar = np.array(Fstar)
                
                # Compute the signal of the planet.
                Fplanet = []
                for j in range(self.Nsplanet):
                    Fplanet += [self.Splanet[j].Compute(self.Filters[i],
                                                        Sys,
                                                        self.Unit,
                                                        self.Mission)]
                Fplanet = np.array(Fplanet)
                
                # Create a new photometry table (if it hasn't already been
                # created) and write the computed fluxes to it.
                if (self.TableFlag == False):
                    temp = self.Filters[i].Name.rfind('/')+1
                    Name = self.PathPlanetTable[:-4]+'_'+self.Filters[i].Name[temp:]
                    self.write(Name,
                               Fstar,
                               Fplanet)
                    self.TableFlag = True
                else:
                    self.append(Name,
                                Fstar,
                                Fplanet)
                
                # Get the next system.
                Sys = self.SysRdr.nextSystem()
            
            print('')
        
        pass
    
    def write(self,
              Name,
              Fstar,
              Fplanet):
        """
        Parameters
        ----------
        Name: str
            Name of the output planet table.
        """
        
        Table = open(Name+'.txt', 'w')
        
        Header = ''
        for i in range(self.Nsstar):
            name = str(self.Sstar[i])
            temp = name.rfind('Photometry')-1
            Header += name[1:temp]+'\t'
        for i in range(self.Nsplanet):
            name = str(self.Splanet[i])
            temp = name.rfind('Photometry')-1
            Header += name[1:temp]+'\t'
        Header += '\n'
        
        # Old header.
        Table.write('Ftherm_star\tFtherm_planet\tFrefl_planet\t\n')
        
        # New header.
        Table.write(Header)
        
        Table.close()
        
        self.append(Name,
                    Fstar,
                    Fplanet)
        
        pass
    
    def append(self,
               Name,
               Fstar,
               Fplanet):
        """
        Parameters
        ----------
        Name: str
            Name of the output planet table.
        
        """
        
        Table = open(Name+'.txt', 'a')
        
        # Write the computed fluxes to the photometry table.
        for i in range(Fstar.shape[1]):
            temp = ''
            for j in range(Fstar.shape[0]):
                temp += '%018.12f\t' % Fstar[j, i]
            for j in range(Fplanet.shape[0]):
                temp += '%018.12f\t' % Fplanet[j, i]
            temp += '\n'
            
            Table.write(temp)
        
        Table.close()
        
        pass
