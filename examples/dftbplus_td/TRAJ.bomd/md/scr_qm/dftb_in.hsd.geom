Geometry = GenFormat{
  <<< 'geometry.gen'
}
Hamiltonian = DFTB{
  SCC = Yes
  SCCTolerance = 1e-06
  MaxSCCIterations = 300
  Mixer = Broyden{}
  Charge = 0.0
  Filling = Fermi{
    Temperature[K] = 0.0
  }
  MaxAngularMomentum = {
    C = 'p'
    H = 's'
  }
  SlaterKosterFiles = Type2FileNames{
    Prefix = '/home/islee/program/dftb/slko/3ob-3-1/'
    Separator = '-'
    Suffix = '.skf'
    LowerCaseTypeName = No
  }
}
Analysis = {
  CalculateForces = Yes
  WriteBandOut = Yes
  WriteEigenvectors = Yes
  MullikenAnalysis = Yes
}
Options = {
  WriteDetailedXml = Yes
  WriteDetailedOut = Yes
  TimingVerbosity = -1
}
ExcitedState = Casida{
  NrOfExcitations = 6
  StateOfInterest = 1
  Symmetry = singlet
  WriteTransitions = Yes
  WriteMulliken = Yes
  WriteXplusY = No
  ExcitedStateForces = No
}
ParserOptions = {
  ParserVersion = 7
}
