# generate NAMD/Gromacs/Colvars config files

class configTemplate:
    ''' generate Colvars config files
        In the Colvars config file, ndx and xyz files are used to indicate the group of atoms
        The non-hydrogen atoms of protein are labeled as 'protein' in complex.ndx
        The non-hydrogen atoms of ligand are labeled as 'ligand' in complex.ndx
        The non-hydrogen atoms of user-defined reference are labeled as 'reference' in complex.ndx '''

    def __init__(self, unit='namd'):
        """initialize template of colvars config

        Args:
            unit (str): 'namd' (A, kcal) or 'gromacs' (nm, kJ). Defaults to 'namd'.
        """
        
        assert(unit == 'namd' or unit == 'gromacs')
        self.unit = unit
        
        # scale factor of energy and length
        # by default, kcal and A are used
        if self.unit == 'namd':
            self.energyFactor = 1.0
            self.lengthFactor = 1.0
        elif self.unit == 'gromacs':
            self.energyFactor = 4.184
            self.lengthFactor = 0.1

    def cvRMSDTemplate(self, setBoundary, lowerBoundary, upperBoundary, refFile):
        """RMSD CV template

        Args:
            setBoundary (bool): whether set boundary (for free-energy calculation)
            lowerBoundary (float): lower boundary of free-energy calculaton
            upperboundary (float): upper boundary of free-energy calculation
            refFile (str): path to the reference file
        
        Returns:
            str: string of RMSD definition
        """

        string = f'\
colvar {{                                    \n\
    name RMSD                                \n'

        if setBoundary:
            string += f'\
    width {0.05 * self.lengthFactor}         \n\
    lowerboundary {lowerBoundary:.1f}        \n\
    upperboundary {upperBoundary:.1f}        \n\
    subtractAppliedForce on                  \n\
    expandboundaries  on                     \n\
    extendedLagrangian on                    \n\
    extendedFluctuation {0.05 * self.lengthFactor}      \n'

        string += f'\
    rmsd {{                                  \n\
        atoms {{                             \n\
            indexGroup  ligand               \n\
        }}                                   \n\
        refpositionsfile  {refFile}          \n\
    }}                                       \n\
}}                                           \n'
        return string
    
    def cvAngleTemplate(self, setBoundary, lowerBoundary, upperBoundary, angle, refFile, oldDefinition = True):
        """Eulaer and polar angle template

        Args:
            setBoundary (bool): whether set boundary (for free-energy calculation)
            lowerBoundary (float): lower boundary of free-energy calculaton
            upperBoundary (float): upper boundary of free-energy calculation
            angle (str): 'eulerTheta', 'eulerPhi', 'eulerPsi', 'polarTheta' or 'polarPhi'
            refFile (str): path to the reference file
            oldDefinition (bool, optional): Whether use old definition of angles
                                            for compatibility. Defaults to True.
        """
        
        assert(
            angle == 'eulerTheta' or angle == 'eulerPhi' or angle == 'eulerPsi' or \
            angle == 'polarTheta' or angle == 'polarPhi'
        )
        
        if angle == 'eulerTheta' or angle == 'eulerPhi' or angle == 'eulerPsi':
            if oldDefinition:
                return self.cvEulerAngleTemplate(setBoundary, lowerBoundary, upperBoundary, angle, refFile)
            else:
                return self.newCvEulerAngleTemplate(setBoundary, lowerBoundary, upperBoundary, angle, refFile)
        elif angle == 'polarTheta' or angle == 'polarPhi':
            if oldDefinition:
                return self.cvPolarAngleTemplate(setBoundary, lowerBoundary, upperBoundary, angle, refFile)
            else:
                return self.newCvPolarAngleTemplate(setBoundary, lowerBoundary, upperBoundary, angle, refFile)

    def cvEulerAngleTemplate(self, setBoundary, lowerBoundary, upperBoundary, angle, refFile):
        """Euler angle template

        Args:
            setBoundary (bool): whether set boundary (for free-energy calculation)
            lowerBoundary (float): lower boundary of free-energy calculaton
            upperboundary (float): upper boundary of free-energy calculation
            angle (str): 'eulerTheta', 'eulerPhi' or 'eulerPsi'
            refFile (str): path to the reference file
            
        Returns:
            string: string of Euler angle definition
        """

        assert(angle == 'eulerTheta' or angle == 'eulerPhi' or angle == 'eulerPsi')
        
        string = f'\
colvar {{                              \n\
    name {angle}                   \n'

        if angle == 'eulerTheta':
            string += f'\
    customFunction asin(2 * (q1*q3-q4*q2)) * 180 / 3.1415926\n'
        elif angle == 'eulerPhi':
            string += f'\
    customFunction atan2(2*(q1*q2+q3*q4), 1-2*(q2*q2+q3*q3)) * 180 / 3.1415926\n'
        elif angle == 'eulerPsi':
            string += f'\
    customFunction atan2(2*(q1*q4+q2*q3), 1-2*(q3*q3+q4*q4)) * 180 / 3.1415926\n'

        if setBoundary:
            string += f'\
    width 1                            \n\
    lowerboundary {lowerBoundary:.1f}      \n\
    upperboundary {upperBoundary:.1f}      \n\
    subtractAppliedForce on            \n\
    expandboundaries  on               \n\
    extendedLagrangian on              \n\
    extendedFluctuation 1              \n'

        string += f'\
    Orientation {{                             \n\
        name  q                                \n\
        atoms {{                               \n\
            indexGroup  ligand                 \n\
            centerReference    on              \n\
            rotateReference    on              \n\
	        enableFitGradients no              \n\
            fittingGroup {{                    \n\
                indexGroup  protein            \n\
            }}                                 \n\
            refpositionsfile  {refFile}        \n\
         }}                                    \n\
         refpositionsfile  {refFile}           \n\
    }}                                         \n\
}}                                             \n'

        return string

    def cvPolarAngleTemplate(self, setBoundary, lowerBoundary, upperBoundary, angle, refFile):
        """Polar angle template

        Args:
            setBoundary (bool): whether set boundary (for free-energy calculation)
            lowerBoundary (float): lower boundary of free-energy calculaton
            upperboundary (float): upper boundary of free-energy calculation
            angle (str): 'polarTheta' or 'polarPhi'
            refFile (str): path to the reference file
            
        Return:
            str: string of polar angle definition
        """

        assert(angle == 'polarTheta' or angle == 'polarPhi')
        
        string = f'\
colvar {{                                   \n\
    name {angle}                            \n'

        if angle == 'polarTheta':
            string += f'\
    customFunction acos(-i2) * 180 / 3.1415926\n'
        elif angle == 'polarPhi':
            string += f'\
    customFunction atan2(i3, i1) * 180 / 3.1415926\n\
    period  360                             \n\
    wrapAround 0.0                          \n'

        if setBoundary:
            string += f'\
    width 1                                 \n\
    lowerboundary {lowerBoundary:.1f}           \n\
    upperboundary {upperBoundary:.1f}           \n\
    subtractAppliedForce on                 \n\
    expandboundaries  on                    \n\
    extendedLagrangian on                   \n\
    extendedFluctuation 1                   \n'

        string += f'\
    distanceDir {{                          \n\
        name  i                             \n\
        group1 {{                           \n\
            indexGroup  reference           \n\
            centerReference    on           \n\
            rotateReference    on           \n\
            enableFitGradients no           \n\
            fittingGroup {{                 \n\
                indexGroup  protein         \n\
            }}                              \n\
            refpositionsfile  {refFile}     \n\
        }}                                  \n\
        group2 {{                           \n\
            indexGroup  ligand              \n\
            centerReference    on           \n\
            rotateReference    on           \n\
            enableFitGradients no           \n\
            fittingGroup {{                 \n\
                indexGroup  protein         \n\
            }}                              \n\
            refpositionsfile  {refFile}     \n\
        }}                                  \n\
    }}                                      \n\
}}                                          \n'
        return string
    
    def newCvEulerAngleTemplate(self, setBoundary, lowerBoundary, upperBoundary, angle, refFile):
        """new definition Euler angle template, probably the pinning down the protein is not required

        Args:
            setBoundary (bool): whether set boundary (for free-energy calculation)
            lowerBoundary (float): lower boundary of free-energy calculaton
            upperboundary (float): upper boundary of free-energy calculation
            angle (str): 'eulerTheta', 'eulerPhi' of 'eulerPsi'
            refFile (str): path to the reference file
            
        Returns:
            string: string of Euler angle definition
        """

        assert(angle == 'eulerTheta' or angle == 'eulerPhi' or angle == 'eulerPsi')
        
        string = f'\
colvar {{                              \n\
    name {angle}                   \n'

        if setBoundary:
            string += f'\
    width 1                            \n\
    lowerboundary {lowerBoundary:.1f}      \n\
    upperboundary {upperBoundary:.1f}      \n\
    subtractAppliedForce on            \n\
    expandboundaries  on               \n\
    extendedLagrangian on              \n\
    extendedFluctuation 1              \n'

        string += f'\
    {angle} {{                             \n\
        atoms {{                               \n\
            indexGroup  ligand                 \n\
            centerReference    on              \n\
            rotateReference    on              \n\
            centerToOrigin     on              \n\
	        enableFitGradients on              \n\
            fittingGroup {{                    \n\
                indexGroup  protein            \n\
            }}                                 \n\
            refpositionsfile  {refFile}        \n\
         }}                                    \n\
         refpositionsfile  {refFile}           \n\
    }}                                         \n\
}}                                             \n'

        return string
    
    def newCvPolarAngleTemplate(self, setBoundary, lowerBoundary, upperBoundary, angle, refFile):
        """new definition of Polar angle template, probably the pinning down the protein is not required

        Args:
            setBoundary (bool): whether set boundary (for free-energy calculation)
            lowerBoundary (float): lower boundary of free-energy calculaton
            upperboundary (float): upper boundary of free-energy calculation
            angle (str): 'polarTheta' or 'polarPhi'
            refFile (str): path to the reference file
            
        Return:
            str: string of polar angle definition
        """

        assert(angle == 'polarTheta' or angle == 'polarPhi')
        
        string = f'\
colvar {{                                   \n\
    name {angle}                            \n'
    
        if angle == 'polarTheta':
            string += f'\
    customFunction acos(-sin(t / 180 * 3.1415926) * sin(p / 180 * 3.1415926)) * 180 / 3.1415926\n'
        elif angle == 'polarPhi':
            string += f'\
    customFunction atan2(cos(t / 180 * 3.1415926), cos(p / 180 * 3.1415926) * sin(t / 180 * 3.1415926)) * 180 / 3.1415926\n\
    period  360                             \n\
    wrapAround 0.0                          \n'

        if setBoundary:
            string += f'\
    width 1                                 \n\
    lowerboundary {lowerBoundary:.1f}           \n\
    upperboundary {upperBoundary:.1f}           \n\
    subtractAppliedForce on                 \n\
    expandboundaries  on                    \n\
    extendedLagrangian on                   \n\
    extendedFluctuation 1                   \n'

        string += f'\
    polarTheta {{                             \n\
        name        t                          \n\
        atoms {{                               \n\
            indexGroup  ligand                 \n\
            centerReference    on              \n\
            rotateReference    on              \n\
            centerToOrigin     on              \n\
            fittingGroup {{                    \n\
                indexGroup  protein            \n\
            }}                                 \n\
            refpositionsfile  {refFile}        \n\
         }}                                    \n\
    }}                                         \n\
    polarPhi {{                             \n\
        name        p                      \n\
        atoms {{                               \n\
            indexGroup  ligand                 \n\
            centerReference    on              \n\
            rotateReference    on              \n\
            centerToOrigin     on              \n\
            fittingGroup {{                    \n\
                indexGroup  protein            \n\
            }}                                 \n\
            refpositionsfile  {refFile}        \n\
         }}                                    \n\
    }}                                         \n\
}}                                             \n'
        return string

    def cvRTemplate(self, setBoundary, lowerBoundary, upperBoundary):
        """r distance template

        Args:
            setBoundary (bool): whether set boundary (for free-energy calculation)
            lowerBoundary (float): lower boundary of free-energy calculaton
            upperboundary (float): upper boundary of free-energy
        
        Returns:
            str: string of distance r definition
        """
        
        string = f'\
colvar {{                            \n\
    name    r                        \n'
        if setBoundary:
            string += f'\
    width {0.1 * self.lengthFactor}      \n\
    lowerboundary {lowerBoundary:.1f}    \n\
    upperboundary {upperBoundary:.1f}    \n\
    subtractAppliedForce on          \n\
    expandboundaries  on             \n\
    extendedLagrangian on            \n\
    extendedFluctuation {0.1 * self.lengthFactor}          \n'

        string += f'\
    distance {{                            \n\
        forceNoPBC       yes               \n\
        group1 {{                          \n\
            indexGroup  reference          \n\
	    }}                                 \n\
        group2 {{                          \n\
            indexGroup  ligand             \n\
        }}                                 \n\
    }}                                     \n\
}}                                         \n'

        return string

    def cvHeadTemplate(self, indexFile):
        """return the head of colvars file

        Args:
            indexFile (str): name of ndx file

        Returns:
            str: head of colvars file
        """
        
        return f'\
colvarsTrajFrequency      5000             \n\
colvarsRestartFrequency   5000            \n\
indexFile                 {indexFile}      \n'

    def cvHarmonicWallsTemplate(self, cv, lowerWall, upperWall):
        ''' template of harmonic wall bias
        
        Args:
            cv (str): name of the colvars
            lowerWall (float): lower wall of the bias
            upperWall (float): upper wall of the bias
                
        Returns:
            str: string of the harmonic wall bias definition '''
            
        string = f'\
harmonicWalls {{                           \n\
    colvars           {cv}                 \n\
    lowerWalls        {lowerWall:.1f}      \n\
    upperWalls        {upperWall:.1f}      \n\
    lowerWallConstant {0.2 * self.energyFactor}      \n\
    upperWallConstant {0.2 * self.energyFactor}      \n\
}}                                         \n'
        return string

    def cvHarmonicTemplate(self, cv, constant, center, tiWindows=0, tiForward=True, targetForceConstant = 0):
        """template for a harmonic restraint

        Args:
            cv (str): name of the colvars
            constant (float): force constant of the restraint (in kcal/mol)
            center (float): center of the restraint
            tiWindows (int): number of windows of the TI simulation (if runs a TI simulation). Defaults to 0.
            tiForward (bool, optional): whether the TI simulation is forward (if runs a TI simulation). Defaults to True.
            targetForceConstant (int, optional): targeted force constant of the restraint in TI simulation
                                                 (if runs a TI simulation)  (in kcal/mol). Defaults to 0.
        
        Returns:
            str: string of the harmonic restraint definition
        """

        string = f'\
harmonic {{                          \n\
    colvars         {cv}             \n\
    forceConstant   {constant * self.energyFactor:.1f}   \n\
    centers         {center:.1f}     \n'
        
        if tiWindows != 0:
            string += f'\
    targetNumSteps      500000                       \n\
    targetEquilSteps    100000                       \n\
    targetForceConstant {targetForceConstant * self.energyFactor}        \n\
    targetForceExponent 4                            \n'

            schedule = ''
            if tiForward:
                schedule += ' '.join([str(float(i) / float(tiWindows)) for i in range(tiWindows+1)])
            else:
                schedule += ' '.join([str(float(i) / float(tiWindows)) for i in range(tiWindows, -1, -1)])
            
            string += f'    lambdaSchedule {schedule}\n'

        string += '}\n'
        return string

    def cvABFTemplate(self, cv):
        ''' template for WTM-eABF bias
        
        Args:
            cv (str): name of the colvars
            
        Returns:
            str: string of the WTM-eABF definition '''
            
        string = f'\
abf {{                            \n\
    colvars        {cv}           \n\
    FullSamples    10000          \n\
    historyfreq    50000          \n\
    writeCZARwindowFile           \n\
}}                                \n\
metadynamics {{                   \n\
    colvars           {cv}        \n\
    hillWidth         3.0         \n\
    hillWeight        {0.05 * self.energyFactor}        \n\
    wellTempered      on          \n\
    biasTemperature   4000        \n\
}}                                \n'
        return string

    def cvProteinTemplate(self, centerCoor, refFile):
        """the template of restraining the protein

        Args:
            centerCoor (np.array, 3): (x,y,z), center of the protein 
            refFile (str): path of the reference file
        
        Returns:
            str: string of the restraining the protein
        """
        
        string = f'\
colvar {{                         \n\
  name translation                \n\
  distance {{                     \n\
    group1 {{                     \n\
      indexGroup  protein         \n\
    }}                            \n\
    group2 {{                     \n\
      dummyAtom ({centerCoor[0] * self.lengthFactor}, {centerCoor[1] * self.lengthFactor}, {centerCoor[2] * self.lengthFactor})    \n\
    }}                            \n\
  }}                              \n\
}}                                \n\
harmonic {{                       \n\
  colvars       translation       \n\
  centers       0.0               \n\
  forceConstant {100.0 * self.energyFactor}    \n\
}}                                \n\
                                  \n\
colvar {{                         \n\
  name orientation                \n\
  orientation {{                  \n\
    atoms {{                      \n\
      indexGroup  protein         \n\
    }}                            \n\
    refPositionsFile   {refFile}  \n\
  }}                              \n\
}}                                \n\
harmonic {{                       \n\
  colvars       orientation       \n\
  centers       (1.0, 0.0, 0.0, 0.0)    \n\
  forceConstant {2000.0 * self.energyFactor}   \n\
}}                                \n'
        return string
