"""
This file contains the methods for processing of a trajectory

"""
from steinhardt import *
import numpy as np
import io
import gzip

def traj_to_systems(filename,natoms,snaps=False,npy=False,zipped=False, usecols=[]):
    """
    Reads in LAMMPS trajectory file and converts it to systems.

    Parameters
    ----------
    
    filename : string
        name of input file
    natoms : int
        number of input atoms
    snaps : bool
        splits the input file into smaller files of one time slice each and
        return a list of filenames if True.
        If False, returns an array of systems with the atom and box coordinates
        read in.
    npy : bool
        If true saves the systems as a npy file.
    usecols : array like
        If custom column numbers are to be used, they can be specified here.
        for example if usecols = [3,4,5,8,9] - the first three columns are
        used to read in x,y, and z values are the rest are read in as custom
        properties. They can be accessed by atom.gcustom(). The minimum length
        of usecols has to be three(for positions). The column number starts 
        from 0. Atom ID always has to be column 0.
        Restrictions - * ONLY used if snaps is set to False *

    Returns
    -------
    
    systems : array like
        array of System. If snaps=True, the atom coordinates are not read in.
        Otherwise, atom and box coordinates are already set.
    
    files : array like
        list of files names of one time slice each. Only returned if
        snaps=True

    """
    if zipped:
        gz = gzip.open(filename,'rb')
        infile = io.BufferedReader(gz)
    else:
        gz = open(filename,'r')
        infile = gz        

    if snaps:
        nblock = natoms+9
        startblock = 0
        count=1
        files = []


        for line in infile:
            if(count==1):
                ff = 'snap.'+str(startblock)+'.dat'
                files.append(ff)
                fout = open(ff,'w')
                fout.write(line)
                #count+=1
            elif(count<nblock):
                fout.write(line)
                #count+=1

            else:
                fout.write(line)
                fout.close()
                count=0
                startblock+=1
            count+=1
        infile.close()

        #qtraj = getQTrajectory(files)
        systems=[]
        for file in files:
            sys = System()
            sys.set_inputfile(file)
            systems.append(sys)
        return systems,files

    else:
        atoms = []
        systems=[]

        custom = False
        if len(usecols) > 0:
            if len(usecols)>=3:
                xid = usecols[0]
                yid = usecols[1]
                zid = usecols[2]
                customcols = usecols[3:]
                custom = True
        else:
            xid = 3
            yid = 4
            zid = 5

        iterval = infile
        count = 1
        timestep = 0
        for line in iterval:
            if count==6:
                box = []
                raw = line.strip().split()
                dimx = float(raw[1])-float(raw[0])          
            
            elif count==7: 
                raw = line.strip().split()
                dimy = float(raw[1])-float(raw[0])
                
            elif count==8:
                raw = line.strip().split()
                dimz = float(raw[1])-float(raw[0])
                boxd = [dimx,dimy,dimz]
                
            elif count>9:
                line = line.strip().split()
                idd = int(line[0]) 
                x = float(line[xid])
                y = float(line[yid])
                z = float(line[zid])
                cavals = []
                if custom:
                    for cindex in customcols:
                        cavals.append(float(line[cindex]))

               
                xx = [x,y,z]
                a = Atom()
                a.sx(xx)
                a.sid(idd)
                a.scustom(cavals)
                #print a.gcustom()  
                atoms.append(a)
                if count==natoms+9:
                    count=0
                    #do stuff here
                    sys = System()
                    sys.assign_particles(atoms,boxd)
                    systems.append(sys)
                    del atoms[:]
                    del boxd[:]
                    atoms = []
                    boxd = []
            count+=1
        infile.close()
        
        if npy:
            np.save(".".join([filename,"npy"]),systems)
        return systems        
