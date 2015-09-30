#!/usr/bin/env python
# -*- coding: latin-1 -*-

#
# fs_shapeDNA
#
# script to compute ShapeDNA of FreeSurfer structures
#
# Original Author: Martin Reuter
# Date: Dec-16-2014
#

import warnings
warnings.filterwarnings('ignore', '.*negative int.*')
import os
import sys
import shlex
import optparse
import logging
import subprocess
import tempfile
import shutil
import time
import math
import stat
import uuid
import fs_shapeDNA

HELPTEXT = """

fs_brainprint.py V2.0
Author: Martin Reuter, 2015

SUMMARY

Computes the BrainPrint for a FreeSurfer subject.

The BrainPrint consists of the shape descriptors (Shape-DNA) [1]
of a selection of both cortical and subcortical structures [2].

Here is a list of structures and FreeSurfer aseg label ids:

CorpusCallosum                  [251, 252, 253, 254, 255]
Cerebellum                      [7, 8, 16, 46, 47]
Ventricles                      [4, 5, 14, 24, 31, 43, 44, 63]
3rd-Ventricle                   [14, 24]
4th-Ventricle                   15
Brain-Stem                      16
Left-Striatum                   [11, 12, 26]
Left-Lateral-Ventricle          [4, 5, 31]
Left-Cerebellum-White-Matter    7
Left-Cerebellum-Cortex          8
Left-Thalamus-Proper            10
Left-Caudate                    11
Left-Putamen                    12
Left-Pallidum                   13
Left-Hippocampus                17
Left-Amygdala                   18
Left-Accumbens-area             26
Left-VentralDC                  28
Right-Striatum                  [50, 51, 58]
Right-Lateral-Ventricle         [43, 44, 63]
Right-Cerebellum-White-Matter   46
Right-Cerebellum-Cortex         47
Right-Thalamus-Proper           49
Right-Caudate                   50
Right-Putamen                   51
Right-Pallidum                  52
Right-Hippocampus               53
Right-Amygdala                  54
Right-Accumbens-area            58
Right-VentralDC                 60

And the following cortical structures:
lh-white-2d    (left white matter surface triangles)
lh-white-3d    (left white matter volume tetrahedra)
lh-pial-2d     (left pial surface triangles)
lh-pial-3d     (left pial+white volume tetrahedra)
rh-white-2d    (same for right hemisphere ...)
rh-white-3d
rh-pial-2d
rh-pial-3d

Note, the 3d (tetrahedra) descriptors are only computed if the 
necessary software (meshfix, gmsh and shapeDNA-tetra) are
available in the $SHAPEDNA_HOME path. As a minimum the file
shapeDNA-tria and the key.txt file need to exist and the 
environment variable $SHAPEDNA_HOME needs to point to that
directory. The key file and shapeDNA-tria can be obtained
from http://reuter.mit.edu/software/shapedna/

If used for a publication, please cite both [1] for the shape
descriptor method and [2] for the application to brain MRI and
definiton of the BrainPrint.


REFERENCES
==========

[1] M. Reuter, F.-E. Wolter and N. Peinecke.
Laplace-Beltrami spectra as "Shape-DNA" of surfaces and solids.
Computer-Aided Design 38 (4), pp.342-366, 2006.
http://dx.doi.org/10.1016/j.cad.2005.10.011

[2] C. Wachinger, P. Golland, W. Kremen, B. Fischl, M. Reuter.
BrainPrint: A discriminative characterization of brain morphology.
NeuroImage Volume 109, pp.232-248, 2015.
http://dx.doi.org/10.1016/j.neuroimage.2015.01.032

"""

def options_parse():
    """
    Command Line Options Parser:
    initiate the option parser and return the parsed object
    """
    parser = optparse.OptionParser(version='$Id: fs_brainPrint,v 1.21 2013/05/17 15:14:08 mreuter Exp $', usage=HELPTEXT)
    
    # help text
    h_sid        = '(REQUIRED) subject ID (FS processed directory inside the subjects directory)'
    h_sdir       = 'FS subjects directory (or set environment $SUBJECTS_DIR)'
    h_num        = 'Number of eigenvalues/vectors to compute (default: 50)'
    h_outdir     = 'Output directory (default: <sdir>/<sid>/surf )'
    h_brainprint = 'Output BrainPrint file (default: <outdir>/<sid>.brainprint_<num>.csv )'
    h_keeptmp    = 'Keep intermediate surface, tet-mesh and ev files'
    h_gsmooth    = 'Geometry smoothing iterations (for surfaces), default 0'
    h_tsmooth    = 'Tangential smoothing iterations (for surface mesh improvement), default 3'
    h_skip3d     = 'Skip 3D tet-meshing and computation'

    parser.add_option('--sid',        dest='sid',        help=h_sid)
    parser.add_option('--sdir',       dest='sdir',       help=h_sdir)
    parser.add_option('--num' ,       dest='num',        help=h_num,     default=50, type='int')
    parser.add_option('--outdir',     dest='outdir',     help=h_outdir)
    parser.add_option('--brainprint', dest='brainprint', help=h_brainprint)
    parser.add_option('--keeptmp',    dest='keeptmp',    help=h_keeptmp, default=False, action='store_true')
    parser.add_option('--tsmooth',    dest='tsmooth',    help=h_tsmooth, default=3, type='int')   
    parser.add_option('--gsmooth',    dest='gsmooth',    help=h_gsmooth, default=0, type='int')   
    parser.add_option('--skip3d',     dest='skip3d',     help=h_skip3d,  default=False, action='store_true')

    (options, args) = parser.parse_args()
    #if len(args) == 0:
    #    parser.print_help()
    #    print '\n'
    #    sys.exit(1)

    # WITHOUT FREESURFER DO NOTHING
    fshome = os.getenv('FREESURFER_HOME')
    if fshome is None:
        parser.print_help()
        print '\nERROR: Environment variable FREESURFER_HOME not set.'
        print '        You need to source FreeSurfer 5.3 or newer.\n'
        sys.exit(1)

    sdnahome = os.getenv('SHAPEDNA_HOME')
    if sdnahome is None:
        parser.print_help()
        print '\nERROR: Environment variable SHAPEDNA_HOME not set.'
        print '       Set that variable to point to the directory containing'
        print '       shapeDNA-tria, e.g.'
        print '       setenv SHAPEDNA_HOME /user/me/shapedna/ (cshell)'
        print '       export SHAPEDNA_HOME=/user/me/shapedna/ (bash)\n'
        sys.exit(1)

    sdna = os.path.join(sdnahome,"shapeDNA-tria")
    if not os.path.exists(sdna):
        parser.print_help()
        print '\nERROR: Cannot find shapeDNA-tria in $SHAPEDNA_HOME\n'
        print '       Set that variable to point to the directory containing'
        print '       shapeDNA-tria, e.g.'
        print '       setenv SHAPEDNA_HOME /user/me/shapedna/ (cshell)'
        print '       export SHAPEDNA_HOME=/user/me/shapedna/ (bash)\n'
        sys.exit(1)
        
    if fs_shapeDNA.which("shapeDNA-tetra") is None and not options.skip3d:
        print '\nWARNING: Cannot find shapeDNA-tetra in $SHAPEDNA_HOME'
        print '       Make sure that binary is at $SHAPEDNA_HOME'
        print '       Switching 3D computation OFF for now (--skip3d)! \n'
        options.skip3d = True
        
    if options.sdir is None:
        options.sdir = os.getenv('SUBJECTS_DIR')

    if options.sdir is None:
        parser.print_help()
        print '\nERROR: specify subjects directory via --sdir or $SUBJECTS_DIR\n'
        sys.exit(1)
        
    if options.sid is None:
        parser.print_help()
        print '\nERROR: Specify --sid\n'
        sys.exit(1)
            
    subjdir = os.path.join(options.sdir,options.sid)
    if not os.path.exists(subjdir):
        parser.print_help()
        print '\nERROR: cannot find sid in subjects directory\n'
        sys.exit(1)
    
    if options.outdir is None:
        options.outdir = os.path.join(subjdir,'brainprint')
    try:
        os.mkdir(options.outdir)
    except OSError, e:
        if e.errno != os.errno.EEXIST:
            raise e
        pass
        
    if options.brainprint is None:
        options.brainprint = os.path.join(options.outdir,options.sid+'.brainprint_'+str(options.num)+'.csv')
    
    return options



def run_cmd(cmd,err_msg):
    """
    execute the comand
    """
    print cmd+'\n'
    args = shlex.split(cmd)
    retcode = subprocess.call(args)
    if retcode != 0 :
        print 'ERROR: '+err_msg
        sys.exit(1)
    print '\n'


def get_evals(evfile):
    """
    returns string list of area, volume and evals
    """
    if not os.path.isfile(evfile) :
        return []
    area = ''
    volume = ''
    evals = []
    with open(evfile, 'r') as inF:
        for line in inF:
            if 'Area:' in line:
                strlst = line.split()
                area = strlst[1]
            if 'Volume:' in line:
                strlst = line.split()
                volume = strlst[1]
            if 'Eigenvalues:' in line:
                evline = inF.next()
                evstr = ''
                while (evline is not None) and (not '}' in evline):
                    evstr = evstr+evline
                    evline = inF.next()
                evstr = evstr+evline
                evstr = evstr.translate(None,'{} \n')
                evals = evstr.split(';')
                evals.insert(0,volume)
                evals.insert(0,area)
                return evals
    return []


def compute_shapeDNAs(options):

    # combined and individual aseg labels
    # combined:
    # Left  Striatum: left  Caudate+Putamen+Accumbens
    # Right Striatum: right Caudate+Putamen+Accumbens
    # CorpusCallosum: (5 sub regions combinded)
    # Cerebellum: brainstem+ (left+right) cerebellum WM and GM
    # Ventricles: (left+right) lat.vent+inf.lat.vent+choroidplexus +3rdVent+CSF
    # Lateral-Ventricle: lat.vent+inf.lat.vent+choroidplexus
    # 3rd-Ventricle: 3rd-Ventricle + CSF
    
    
    ## Original BrainPrint definition, but not everything of this was used in Wachinger, 2015
    ## some structures such as choroid-plexus were excluded in analysis scripts
    ## names for table output:
    #structures = ['Left-Striatum','Right-Striatum','CorpusCallosum','Cerebellum','Ventricles',
    #              'Left-Lateral-Ventricle','Left-Inf-Lat-Vent','Left-Cerebellum-White-Matter',
    #              'Left-Cerebellum-Cortex','Left-Thalamus-Proper','Left-Caudate','Left-Putamen',
    #              'Left-Pallidum','3rd-Ventricle','4th-Ventricle','Brain-Stem','Left-Hippocampus',
    #              'Left-Amygdala','CSF','Left-Accumbens-area','Left-VentralDC','Left-choroid-plexus',
    #              'Right-Lateral-Ventricle','Right-Inf-Lat-Vent','Right-Cerebellum-White-Matter',
    #              'Right-Cerebellum-Cortex','Right-Thalamus-Proper','Right-Caudate','Right-Putamen',
    #              'Right-Pallidum','Right-Hippocampus','Right-Amygdala','Right-Accumbens-area',
    #              'Right-VentralDC','Right-choroid-plexus',
    #              'lh-white-2d','lh-white-3d','lh-pial-2d','lh-pial-3d',
    #              'rh-white-2d','rh-white-3d','rh-pial-2d','rh-pial-3d']
    ## label ids for aseg structures
    #labels = [[11, 12, 26], [50, 51, 58], [251, 252, 253, 254, 255], [7, 8, 16, 46, 47], [4, 5, 14, 24, 31, 43, 44, 63], 4, 5, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 24, 26, 28, 31, 43, 44, 46, 47, 49, 50, 51, 52, 53, 54, 58, 60, 63]

    # new definition for better ventricle treatment
    # (i)  choroid plexus is partially in Lat.Vent and Inf.Lat.Vent, so we merge those three (call em Lateral-Ventricle)
    # (ii) merge 3rd Ventricle and CSF (separation not clear), call 3rd-Ventricle
    # now we therefore drop Inf.Lat.Vent,ChoroidPlexus,CSF
    # maybe exclude brainstem and ventralDC ?
    # names for table output:
    structures = ['CorpusCallosum','Cerebellum','Ventricles',
                  '3rd-Ventricle','4th-Ventricle','Brain-Stem',
                  'Left-Striatum','Left-Lateral-Ventricle',
                  'Left-Cerebellum-White-Matter','Left-Cerebellum-Cortex',
                  'Left-Thalamus-Proper','Left-Caudate','Left-Putamen',
                  'Left-Pallidum','Left-Hippocampus','Left-Amygdala',
                  'Left-Accumbens-area','Left-VentralDC',
                  'Right-Striatum','Right-Lateral-Ventricle',
                  'Right-Cerebellum-White-Matter','Right-Cerebellum-Cortex',
                  'Right-Thalamus-Proper','Right-Caudate','Right-Putamen',
                  'Right-Pallidum','Right-Hippocampus','Right-Amygdala',
                  'Right-Accumbens-area','Right-VentralDC',
                  'lh-white-2d','lh-white-3d','lh-pial-2d','lh-pial-3d',
                  'rh-white-2d','rh-white-3d','rh-pial-2d','rh-pial-3d']
    # label ids for aseg structures
    labels = [[251, 252, 253, 254, 255], [7, 8, 16, 46, 47], [4, 5, 14, 24, 31, 43, 44, 63],
              [14, 24], 15, 16,
              [11, 12, 26], [4, 5, 31],
               7, 8,
               10, 11, 12,
               13, 17, 18,
               26, 28,
              [50, 51, 58],[43, 44, 63],
               46, 47,
               49, 50, 51,
               52, 53, 54,
               58, 60]
    


    class sdnaopt(object):
        num      = options.num
        degree   = 1
        refmin   = 3000
        ignorelq = True
        bcond    = 1 #Neumann (for tet mesh)
        tsmooth  = options.tsmooth
        gsmooth  = options.gsmooth
        evec     = False
        param2d  = None
 
    evmat = []
 
    for lab in labels:
        if type(lab) == list:
            astring  = '_'.join(str(x) for x in lab)
        else:
            astring = str(lab)
        surfnamei = 'aseg.init.'+astring+'.vtk'
        asegsurfi  = os.path.join(options.outdir,surfnamei)
        surfnameo = 'aseg.final.'+astring+'.vtk'
        asegsurfo  = os.path.join(options.outdir,surfnameo)
        outev    = asegsurfo+'.ev'
        failed = False
        evs = []
        try:
            fs_shapeDNA.get_aseg_surf(options.sdir,options.sid,astring.split('_'),asegsurfi)
            fs_shapeDNA.run_shapeDNAtria(asegsurfi,outev,asegsurfo,sdnaopt)
            evs = get_evals(outev)
        except subprocess.CalledProcessError as e:
            print 'Error occured, skipping label '+astring
            failed = True
            
        if not evs or failed:
            evs = ['NaN'] * (sdnaopt.num+2)
        evmat.append(evs)
        if not options.keeptmp and not failed:
            cmd ='rm '+asegsurfi
            run_cmd(cmd,'rm temp asegsurfi failed?')
            cmd ='rm '+asegsurfo
            run_cmd(cmd,'rm temp asegsurfo failed?')
            cmd ='rm '+outev
            run_cmd(cmd,'rm temp outev failed?')
       
        
        #lstring  = ','.join(lab)
        #cmd = 'fs_shapeDNA.py --sid '+options.sid+' --sdir '+options.sdir+' --asegid '+lstring+' --num '+options.num
        #if options.outdir is not None:
        #    cmd = cmd+' --outdir '+options.outdir
        #run_cmd(cmd,'fs_shapeDNA.py '+lstring+' failed?')

    # Surfaces (both 2D and 3D tet):
    for hem in ['lh','rh']:
        for typeSurf in ['white', 'pial']:
            surfname = hem+'.'+typeSurf
            insurf   = os.path.join(options.sdir,options.sid,'surf',surfname)
            outsurf  = os.path.join(options.outdir,surfname+'.final.vtk')
            outev2d  = os.path.join(options.outdir,surfname+'.ev')
            outtet   = os.path.join(options.outdir,surfname+'.msh')
            outev3d  = os.path.join(options.outdir,surfname+'.msh.ev')

            failed = False
            try:
                fs_shapeDNA.run_shapeDNAtria(insurf,outev2d,outsurf,sdnaopt)
                evs = get_evals(outev2d)
            except subprocess.CalledProcessError as e:
                print 'Error occured, skipping 2D surface '+surfname
                failed = True

            if not evs or failed:
                evs = ['NaN'] * (sdnaopt.num+2)
            evmat.append(evs)
            if not options.keeptmp and not failed:
                cmd ='rm '+outev2d
                run_cmd(cmd,'rm temp outev2d failed?')
                cmd ='rm '+outsurf
                run_cmd(cmd,'rm temp outsurf failed?')
            
            if options.skip3d:
                print 'Skipping 3D meshing and computation'
                continue

            if fs_shapeDNA.which('meshfix') is None or fs_shapeDNA.which('gmsh') is None or fs_shapeDNA.which('shapeDNA-tetra') is None:
                print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                print '!!! Skipping 3D computations due to missing executables (meshfix, gmsh or shapeDNA-tetra) !!!'
                print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
                continue
            # try with 3 mesh fix iterations:    
            failed = False
            try:
                fixiter = 3
                fs_shapeDNA.get_tetmesh(insurf,outtet,fixiter)
                fs_shapeDNA.run_shapeDNAtetra(outtet,outev3d,sdnaopt)
                evs = get_evals(outev3d)
            except subprocess.CalledProcessError as e:
                print 'Error occured, skipping 3D surface '+surfname
                failed = True

            # if failed, try with 4 mesh fix iterations:    
            if not evs or failed:
                failed = False
                try:
                    fixiter = 4
                    fs_shapeDNA.get_tetmesh(insurf,outtet,fixiter)
                    fs_shapeDNA.run_shapeDNAtetra(outtet,outev3d,sdnaopt)
                    evs = get_evals(outev3d)
                except subprocess.CalledProcessError as e:
                    print 'Error occured, skipping 3D surface '+surfname
                    failed = True

            if not evs or failed:
                evs = ['NaN'] * (sdnaopt.num+2)
            evmat.append(evs)
            if not options.keeptmp and not failed:
                cmd ='rm '+outev3d
                run_cmd(cmd,'rm temp outev3d failed?')
                cmd ='rm '+outtet
                run_cmd(cmd,'rm temp outtet failed?')
                    
            #cmd = 'fs_shapeDNA.py --sid '+options.sid+' --sdir '+options.sdir+' --surf '+surfname+' --num '+options.num
            #if options.outdir is not None:
            #    cmd = cmd+' --outdir '+options.outdir
            #run_cmd(cmd,'fs_shapeDNA.py --surf '+sstring+' failed?')
            #sstring=hem+'.'+typeSurf
            #cmd = 'fs_shapeDNA.py --sid '+options.sid+' --sdir '+options.sdir+' --surf '+surfname+' --num '+options.num+' --dotet'
            #if options.outdir is not None:
            #    cmd = cmd+' --outdir '+options.outdir
            #run_cmd(cmd,'fs_shapeDNA.py --dotet --surf '+sstring+' failed?')

    return (structures, evmat)


def write_evs(outfile,structures,evmat):
    # write final csv
    text_file = open(outfile, "w")
    text_file.write((','.join(structures))+'\n')
    evstrans= zip(*evmat)
    for item in evstrans:
        text_file.write("%s\n" % ','.join(item))
    text_file.close()
    

if __name__=="__main__":
    # Command Line options and error checking done here
    options = options_parse()

    (structures , evmat) = compute_shapeDNAs(options)
    write_evs(options.brainprint,structures,evmat)

    # always exit with 0 exit code
    sys.exit(0)
