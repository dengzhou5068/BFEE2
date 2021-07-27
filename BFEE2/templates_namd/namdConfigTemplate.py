# generate NAMD/Gromacs/Colvars config files

class configTemplate:
    ''' generate NAMD config files '''

    def __init__(self):
        pass

    def namdConfigTemplate(
                            self,
                            forceFieldType,
                            forceFieldFiles,
                            topFile,
                            coorFile,
                            NAMDRestartCoor,
                            NAMDRestartVel,
                            NAMDRestartXsc,
                            PBCCondition,
                            outputPrefix,
                            temperature,
                            numSteps,
                            cvFile = '',
                            cvDefinitionFile = '',
                            CVRestartFile = '',
                            fepFile = '',
                            fepWindowNum = 20,
                            fepForward = True,
                            fepDoubleWide = False,
                            fepMinBeforeSample = False,
                            membraneProtein = False
                            ):
        """the namd config file template

        Args:
            forceFieldType (str): 'charmm' or 'amber'
            forceFieldFiles (list of str): name of charmm force field files
            topFile (str): name of the topology file (psf, parm) 
            coorFile (str): name of the coordinate file (pdb, rst)
            NAMDRestartCoor (str): name of namd binary restart coor file (if restart from a previous simulation)
            NAMDRestartVel (str): name of namd binary restart vel file
            NAMDRestartXsc (str): name of namd restart xsc file
            PBCCondition (np.array, 2*3): PBC vector, ((lengthX, lengthY, lengthZ),(centerX, centerY, centerZ))
            outputPrefix (str): prefix of output file
            temperature (float): temperature of the simulation
            numSteps (int): number of steps of the simulation
            cvFile (str): name of Colvars file. Defaults to ''.
            cvDefinitionFile (str, optional): name of a TCL file defining new CVs. Defaults to ''.
            CVRestartFile (str, optional): name of Colvars restart file. Defaults to ''.
            fepFile (str, optional): name of fep file, indicating which atoms will be generated/removed 
                                     (if run alchemical simulation). Defaults to ''.
            fepWindowNum (int, optional): number of fep windows. Defaults to 20.
            fepForward (bool, optional): whether this is a forward fep simulation. Defaults to True.
            fepDoubleWide (bool, optional): whether this is a double-wide fep simulation. Defaults to False.
            fepMinBeforeSample (bool, optional): whether do minimization before sampling in each FEP window.
                                                 Defaults to False.
            membraneProtein (bool, optional): whether simulating a membrame protein. Defaults to False.
            
        Returns:
            str: a NAMD config string if succeed, and empty string otherwise
        """

        assert(forceFieldType == 'charmm' or forceFieldType == 'amber')

        configString = f'\
coordinates    {coorFile}                   \n'

        # force field files
        if forceFieldType == 'charmm':
            configString += f'\
structure      {topFile}                \n\
paraTypeCharmm    on                    \n'
            for ff in forceFieldFiles:
                configString += f'\
parameters    {ff}                      \n'
        elif forceFieldType == 'amber':
            configString += f'\
parmFile      {topFile}                 \n\
amber    on                             \n'
        else:
            # error
            return ''

        # structure
        if forceFieldType == 'charmm':
            configString += f'\
exclude    scaled1-4                    \n\
1-4scaling    1.0                       \n\
switching            on                 \n\
switchdist           10.0               \n\
cutoff               12.0               \n\
pairlistdist         14.0               \n'
        elif forceFieldType == 'amber':
            configString += f'\
exclude    scaled1-4                    \n\
1-4scaling    0.83333333                \n\
switching            on                 \n\
switchdist           8.0                \n\
cutoff               9.0                \n\
pairlistdist         11.0               \n'
        else:
            # error
            return ''

        if NAMDRestartCoor == '' and NAMDRestartVel == '' and NAMDRestartXsc == '' and PBCCondition != '':
            # set temperature
            configString += f'\
temperature    {temperature}                      \n\
cellBasisVector1    {PBCCondition[0][0]} 0 0         \n\
cellBasisVector2    0 {PBCCondition[0][1]} 0         \n\
cellBasisVector3    0 0 {PBCCondition[0][2]}         \n\
cellOrigin    {PBCCondition[1][0]} {PBCCondition[1][1]} {PBCCondition[1][2]}  \n'
        elif NAMDRestartCoor != '' and NAMDRestartVel != '' and NAMDRestartXsc != '' and PBCCondition == '':
            # restart from files
            configString += f'\
bincoordinates    {NAMDRestartCoor}                          \n\
binvelocities    {NAMDRestartVel}                            \n\
ExtendedSystem    {NAMDRestartXsc}                           \n'
        else:
            # error
            return ''

        # other parameters
        configString += f'\
binaryoutput         yes                        \n\
binaryrestart        yes                        \n\
outputname           {outputPrefix}             \n\
dcdUnitCell          yes                        \n\
outputenergies       5000                       \n\
outputtiming         5000                       \n\
outputpressure       5000                       \n\
restartfreq          5000                       \n\
XSTFreq              5000                       \n\
dcdFreq              5000                       \n\
hgroupcutoff         2.8                        \n\
wrapAll              off                        \n\
wrapWater            on                         \n\
langevin             on                         \n\
langevinDamping      1                          \n\
langevinTemp         {temperature}              \n\
langevinHydrogen     no                         \n\
langevinpiston       on                         \n\
langevinpistontarget 1.01325                    \n\
langevinpistonperiod 200                        \n\
langevinpistondecay  100                        \n\
langevinpistontemp   {temperature}              \n\
usegrouppressure     yes                        \n\
PME                  yes                        \n\
PMETolerance         10e-6                      \n\
PMEInterpOrder       4                          \n\
PMEGridSpacing       1.0                        \n\
timestep             2.0                        \n\
fullelectfrequency   2                          \n\
nonbondedfreq        1                          \n\
rigidbonds           all                        \n\
rigidtolerance       0.00001                    \n\
rigiditerations      400                        \n\
stepspercycle        10                         \n\
splitpatch           hydrogen                   \n\
margin               2                          \n'

        # membrane protein
        if membraneProtein:
            configString += f'\
useflexiblecell      yes                        \n\
useConstantRatio     yes                        \n'
        else:
            configString += f'\
useflexiblecell      no                         \n\
useConstantRatio     no                         \n'

        # colvars definition
        if cvFile != '':
            configString += f'\
colvars    on                                   \n\
colvarsConfig    {cvFile}                       \n'

            if CVRestartFile != '':
                configString += f'\
colvarsInput     {CVRestartFile}                \n'


        if cvDefinitionFile != '':
            configString += f'\
source     {cvDefinitionFile}                   \n'

        # fep
        if fepFile == '':
            if NAMDRestartCoor == '' and NAMDRestartVel == '' and NAMDRestartXsc == '':
                configString += f'\
minimize    500                                 \n\
reinitvels    {temperature}                     \n'
            configString += f'\
run    {numSteps}                               \n'
        else:
            # currently the alchemical route is somewhat hard-coded
            # this will be improved in the future
            configString += f'\
source ../fep.tcl                                  \n\
alch on                                         \n\
alchType FEP                                    \n\
alchFile {fepFile}                              \n\
alchCol B                                       \n\
alchOutFile {outputPrefix}.fepout               \n\
alchOutFreq 50                                  \n\
alchVdwLambdaEnd 0.7                            \n\
alchElecLambdaStart 0.5                         \n\
alchEquilSteps 100000                           \n'

            if fepForward:
                if not fepDoubleWide:
                    if fepMinBeforeSample:
                        # minimize before sampling
                        configString += f'\
runFEPmin 0.0 1.0 {1.0/fepWindowNum} 500000 1000 {temperature}\n'
                    else:
                        configString += f'\
runFEP 0.0 1.0 {1.0/fepWindowNum} 500000\n'

                else:
                    # double wide simulation
                    configString += f'\
runFEP 0.0 1.0 {1.0/fepWindowNum} 500000 true\n'

            else:
                # backward
                if not fepDoubleWide:
                    if fepMinBeforeSample:
                        # minimize before sampling
                        configString += f'\
runFEPmin 1.0 0.0 {-1.0/fepWindowNum} 500000 1000 {temperature}\n'
                    else:
                        configString += f'\
runFEP 1.0 0.0 {-1.0/fepWindowNum} 500000\n'

                else:
                    # double wide simulation
                    configString += f'\
runFEP 1.0 0.0 {-1.0/fepWindowNum} 500000 true\n'

        return configString