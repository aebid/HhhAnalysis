import os, ROOT
from math import hypot, pi

# ========= UTILITIES =======================


def deltaPhi(phi1, phi2):
    # Catch if being called with two objects
    if type(phi1) != float and type(phi1) != int:
        phi1 = phi1.phi
    if type(phi2) != float and type(phi2) != int:
        phi2 = phi2.phi
    # Otherwise
    dphi = (phi1 - phi2)
    while dphi > pi:
        dphi -= 2 * pi
    while dphi < -pi:
        dphi += 2 * pi
    return dphi


def deltaR(eta1, phi1, eta2=None, phi2=None):
    # catch if called with objects
    if eta2 == None:
        return deltaR(eta1.eta, eta1.phi, phi1.eta, phi1.phi)
    # otherwise
    return hypot(eta1 - eta2, deltaPhi(phi1, phi2))



def ensureTFile(filename,option='READ',verbose=False):
  """Open TFile, checking if the file in the given path exists."""
  if not os.path.isfile(filename):
    raise IOError("File in path '%s' does not exist!"%(filename))
  if verbose:
    print("Opening '%s'..."%(filename))
  file = ROOT.TFile.Open(filename,option)
  if not file or file.IsZombie():
    raise IOError("Could not open file by name '%s'"%(filename))
  return file

def extractTH1(file,histname,setdir=True):
  """Get histogram by name from a given file."""
  close = False
  if isinstance(file,str):
    file  = ensureTFile(file,'READ')
    close = True
  if not file or file.IsZombie():
    raise IOError("Could not open file for histogram '%s'!"%(histname))
  hist = file.Get(histname)
  if not hist:
    raise IOError("Did not find histogram '%s' in file '%s'!"%(histname,file.GetName()))
  if setdir and isinstance(hist,ROOT.TH1):
    hist.SetDirectory(0)
    if close:
      file.Close()
  return hist
