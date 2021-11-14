import subprocess, sys, datetime
from functools import reduce

scratchDir = '/tmp/'
myOrca     = 'orca'

def main() -> None:
    if len(sys.argv) > 0:
        print(" * Input file is ", sys.argv[0])
    else:
        raise ValueError("Give the name of your input file! Fowlloing are the instructions. \n 1) Input file name should not include any dot except the one at the end e.g., ch2-rhf.inp is OK but ch2.rhf.inp is NOT OK!\n 2) You should call RunOrca in the same directory as the input file e.g., RunOrca ./o2.inp (or simply o2.inp) is OK but RunOrca ../o2.inp is NOT OK!")
    myLabel : str = reduce(lambda x,y: x + '.' + y, (sys.argv[0]).split('.')[:-1])

    # Get current time to make foot note
    time = datetime.datetime.now()
    myScratchDir : str = scratchDir + myLabel + '-' + str(time) + '/'

    print(" * Scratch directory : ", myScratchDir)

    # Create the scratch directory
