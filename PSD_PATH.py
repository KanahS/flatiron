import glob, os

#repertory="path"


def PSD_PATH_IN_REPERTORY(repertory):
    os.chdir(repertory)
    path_tab=glob.glob("kplr_*.fits")
    return path_tab

def BACK_PATH_IN_REPERTORY(repertory):
    os.chdir(repertory)
    path_tab=glob.glob("*/back_*.fits")
    return path_tab

def BACK_CORRECT_PATH_IN_REPERTORY(repertory):
    os.chdir(repertory)
    path_tab=glob.glob("*/back_*_corrected.fits")
    return path_tab

def L0_PATH_IN_REPERTORY(repertory):
    os.chdir(repertory)
    path_tab=glob.glob("*/*.pkb")
    return path_tab

if __name__ == '__main__' :
  print(PSD_PATH_IN_REPERTORY(repertory))
