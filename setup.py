"""Setup settings for distutils.

You can get debug build with command:
python2.2 setup.py debug build
"""
import sys





##
##   $Log$
##   Revision 1.10  2004/07/14 12:04:46  kpalin
##   New version for genome wide alignments.
##
##   Revision 1.9  2004/04/14 07:54:10  kpalin
##   Checking in new version code for distribution.
##
##



print sys.argv

common_compile_args=["-Wall","-O3"]
if len(sys.argv)>1 and sys.argv[1]=='debug':
    common_compile_args=["-O0","-fno-inline","-Wall","-g","-UNDEBUG","-DEXTRADEBUG","-DDEBUG_OUTPUT"]
    print "Using debug settings"
    del sys.argv[1]


from distutils.core import setup, Extension



modMatrix = Extension('matrix',
                    libraries = ["stdc++"],
                    sources = ['src/matrix.cc'],
                    extra_compile_args = common_compile_args+["-DSEQ_BUFFER_SIZE=5000000"],
                    extra_link_args = [])



#alignCompileArgs = ["-Wall","-g","-UNDEBUG","-DEXTRADEBUG"]
#alignCompileArgs = ["-fno-inline","-Wall","-g","-UNDEBUG"]
#alignCompileArgs = ["-g","-Wall","-O3"]
alignCompileArgs = []

#Save memory if the matrix would be larger than 512MB
alignCompileArgs = ["-DSAVE_MEM","-DSAVE_MEM_LIMIT=536870912"]

alignLibDirs=[]
alignLibs=["stdc++"]
try:
    import gzip
    print "Looks like you have zlib! (It's a bad thing if you don't)"
    alignCompileArgs.extend(["-Isrc/gzstream","-DHAVE_GZSTREAM=1"])
    alignLibDirs.append("src/gzstream")
    alignLibs.extend(["gzstream","z","m"])
except ImportError:
    print "Looks like you don't have zlib! You will not be able to use gzip:ed files in alignment"
    pass


modAlign = Extension('align',
                     library_dirs = alignLibDirs,
                     libraries = alignLibs,
                     sources = ['src/align.cc'],
                     extra_compile_args=alignCompileArgs+common_compile_args,
                     extra_link_args = [])

modAlignedCols=Extension('alignedCols',
                         library_dirs = alignLibDirs,
                         libraries = alignLibs,
                         sources = ['src/alignedCols.cc' ],
                         extra_compile_args=alignCompileArgs+common_compile_args,
                         extra_link_args = [])
                         

modDist = Extension("editdist",
                    sources = ["src/editdist.c"],
                    extra_compile_args = common_compile_args
                    )


modMultiAlign = Extension('multiAlign',
                     library_dirs = alignLibDirs,
                     libraries = alignLibs,
                     sources = ['src/multiAlign.cc'],
                     extra_compile_args=alignCompileArgs+common_compile_args,
                     extra_link_args = [])


ext_modList= [modMatrix,modAlignedCols,modAlign,modMultiAlign,modDist]
#ext_modList= [modAlignedCols,modMultiAlign,modAlign]

setup (name = 'mabs',
       version = '1_beta15',
       url = "http://www.cs.helsinki.fi/u/kpalin/",
       author = "Kimmo Palin, Matthias Berg",
       author_email = "kimmo.palin@helsinki.fi",
       maintainer = "Kimmo Palin",
       maintainer_email = "kimmo.palin@helsinki.fi",
       license = "GPL (see file COPYING)",
       description = 'c++ extension modules:\na binding site matrix\nan alignment function',
       ext_modules = ext_modList,
       packages = ["mabslib"],
       scripts = [ "mabs"] )
