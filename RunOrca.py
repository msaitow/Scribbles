import subprocess, sys, datetime
from functools import reduce

scratchDir = '/tmp/'
myOrca     = 'orca'

def main() -> None:
    if len(sys.argv) > 1:
        print(" * Input file is ", sys.argv[1])
    else:
        raise ValueError("Give the name of your input file! Fowlloing are the instructions. \n 1) Input file name should not include any dot except the one at the end e.g., ch2-rhf.inp is OK but ch2.rhf.inp is NOT OK!\n 2) You should call RunOrca in the same directory as the input file e.g., RunOrca ./o2.inp (or simply o2.inp) is OK but RunOrca ../o2.inp is NOT OK!")
    myLabel : str = reduce(lambda x,y: x + '.' + y, (sys.argv[1]).split('.')[:-1])

    # Get current time to make foot note
    time = datetime.datetime.now()
    myScratchDir : str = scratchDir + myLabel + '-' + str(time).replace(' ','-') + '/'

    print(" * Scratch directory : ", myScratchDir)

    # Create the scratch directory
    subprocess.call('mkdir ' + myScratchDir, shell=True)

    # Copy the input file into the scratch
    if len(myLabel) == 0: raise ValueError('Input file name should end by \".inp\"!')
    else: subprocess.call('cp ./' + myLabel + '.inp ' + myScratchDir, shell=True)

    # Submit a calculation
    subprocess.call(myOrca + ' ' + myScratchDir + myLabel + '.inp | tee ./' + myLabel + '.out', shell=True)
    print('Calculation finished!')

    # Convert the gbw file into molden file
    subprocess.call('mv ' + myScratchDir + myLabel + '.gbw' + ' ./', shell=True)
    subprocess.call(myOrca + '_2mkl ' + myLabel + ' -molden', shell=True)
    subprocess.call('mv ' + myLabel + '.molden.input' + ' ./' + myLabel + '.molden', shell=True)
    subprocess.call('rm ' + myLabel + '.gbw', shell=True)

    print('Conversion of GBW files into molden file!')

    
if __name__ == '__main__':
    main()
    
