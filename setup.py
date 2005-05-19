"""Setup settings for distutils.

You can get debug build with command:
python2.2 setup.py debug build
"""
import sys





##
##   $Log$
##   Revision 1.18.2.1  2005/05/09 07:12:52  kpalin
##   Under development.
##
##   Revision 1.18  2005/03/02 13:31:43  kpalin
##   One more working version.
##
##   Revision 1.17  2005/01/14 13:45:21  kpalin
##   Better installation and added msvcp71.dll and msvcr71.dll for windows
##   installation.
##
##   Revision 1.16  2005/01/14 12:51:48  kpalin
##   Fixes for TCL/TIX gui in windows
##
##   Revision 1.15  2005/01/07 13:41:25  kpalin
##   Works with py2exe. (windows executables)
##
##   Revision 1.14  2005/01/05 09:05:08  kpalin
##   Fixed a one line description to satisfy bdist_rpm
##
##   Revision 1.13  2004/12/22 11:14:16  kpalin
##   Some fixes for better distributability
##
##   Revision 1.12  2004/12/14 13:07:52  kpalin
##
##   Name change from MABS to EEL (Enhancer Element Locator / Monty Python pun
##   "My hovercraft is full of EELs" )
##
##   Revision 1.11  2004/07/30 12:22:10  kpalin
##   Exact multiple alignment and alignedCols.alnColumn
##
##   Revision 1.10  2004/07/14 12:04:46  kpalin
##   New version for genome wide alignments.
##
##   Revision 1.9  2004/04/14 07:54:10  kpalin
##   Checking in new version code for distribution.
##
##



print sys.argv

common_compile_args=["-Wall"]
if len(sys.argv)>1:
    if sys.argv[1]=='debug':
        common_compile_args=["-O0","-fno-inline","-Wall","-g","-UNDEBUG","-DEXTRADEBUG","-DDEBUG_OUTPUT","-DSUBOPTDEBUG"]
        print "Using debug settings"
        del sys.argv[1]
    elif sys.argv[1]=='profile':
        common_compile_args=["-pg"]
        print "Using profiling settings"
        del sys.argv[1]


from distutils.core import setup, Extension


#alignCompileArgs = ["-Wall","-g","-UNDEBUG","-DEXTRADEBUG"]
#alignCompileArgs = ["-fno-inline","-Wall","-g","-UNDEBUG"]
#alignCompileArgs = ["-g","-Wall","-O3"]
alignCompileArgs = []

#Save memory if the matrix would be larger than 512MB
alignCompileArgs = ["-DSAVE_MEM","-DSAVE_MEM_LIMIT=536870912"]
alignLibs=[]


commonLibs=[]
alignLibDirs=[]
compileLibs=[]



extra_data=[]

from glob import glob
if sys.platform=='win32':
    try:
        import py2exe
        SDKdir=r"C:\Program Files\Microsoft.NET\SDK\v1.1\Bin"
        extra_data.append((r"tcl\tix8.1",[x for x in glob(r"C:\Python*\tcl\tix8.1\*") if x[-7:] not in ("bitmaps",".1\\pref")]))
        extra_data.append((r"tcl\tix8.1\bitmaps",glob(r"C:\Python*\tcl\tix8.1\bitmaps\*")))
        extra_data.append((r"tcl\tix8.1\pref",glob(r"C:\Python*\tcl\tix8.1\pref\*")))
        extra_data.append((".",glob(r"C:\Python*\DLLs\tix*.dll")+[SDKdir+r"\msvcr71.dll",SDKdir+r"\msvcp71.dll"]))
    except ImportError:
        pass
    common_compile_args.extend([r"/O2", r"/IC:\Program Files\Microsoft Platform SDK for Windows XP SP2\Include",\
                                r"/IC:\Program Files\Microsoft Visual C++ Toolkit 2003\include"])
    alignCompileArgs.extend([r"/DHAS_NO_TRUNC"])
    pass
else:
    commonLibs.extend(["stdc++"])
    #common_compile_args.extend([r"-O3", r"-Wall"])
    try:
        import gzip
        print "Looks like you have zlib! (It's a bad thing if you don't)"
        alignCompileArgs.extend(["-Isrc/gzstream","-DHAVE_GZSTREAM=1"])
        alignLibDirs.append("src/gzstream")
        alignLibs.extend(["gzstream","z","m"])

        libgzSrcs= ["src/gzstream/gzstream.C"]
        compileLibs.append(('gzstream',
                         {'sources':libgzSrcs,
                          'include_dirs':['./src/gzstream/'
                                          ],
                          'macros':[]},
                         ))

    except ImportError:
        print "Looks like you don't have zlib! You will not be able to use gzip:ed files in alignment"
        pass





modMatrix = Extension('eellib._c_matrix',
                    libraries = commonLibs,
                    sources = ['src/_c_matrix.cc'],
                    extra_compile_args = common_compile_args+["-DSEQ_BUFFER_SIZE=5000000"],
                    extra_link_args = [])




modAlign = Extension('eellib.align',
                     library_dirs = alignLibDirs,
                     libraries = commonLibs+alignLibs,
                     sources = ['src/align.cc'],
                     extra_compile_args=alignCompileArgs+common_compile_args,
                     extra_link_args = [])

modAlignedCols=Extension('eellib.alignedCols',
                         library_dirs = alignLibDirs,
                         libraries = commonLibs+alignLibs,
                         sources = ['src/alignedCols.cc' ],
                         extra_compile_args=alignCompileArgs+common_compile_args,
                         extra_link_args = [])
                         

modDist = Extension("eellib.editdist",
                    sources = ["src/editdist.c"],
                    extra_compile_args = common_compile_args
                    )


modMultiAlign = Extension('eellib.multiAlign',
                          library_dirs = alignLibDirs,
                          libraries = commonLibs+alignLibs,
                          sources = ['src/multiAlign.cc'],
                          extra_compile_args=alignCompileArgs+common_compile_args,
                          extra_link_args = [])


ext_modList= [modMatrix,modAlignedCols,modAlign,modDist]
#ext_modList= [modMatrix,modAlignedCols,modAlign,modMultiAlign,modDist]

setup (name = 'EEL',
       version = '1.0',
       url = "http://www.cs.helsinki.fi/u/kpalin/",
       author = "Kimmo Palin, Matthias Berg",
       author_email = "kimmo.palin@helsinki.fi",
       maintainer = "Kimmo Palin",
       maintainer_email = "kimmo.palin@helsinki.fi",
       license = "GPL (see file COPYING)",
       description = 'A software tool for locating enhancer elements from genomic sequence',
       ext_modules = ext_modList,
       libraries= compileLibs,
       packages = ["eellib"],
       scripts = [ "eel"],
       console = [ "eel"],
       data_files=extra_data)
