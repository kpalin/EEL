from distutils.core import setup, Extension


modMatrix = Extension('matrix',
                    libraries = ["stdc++"],
                    sources = ['src/matrix.cc'],
                    extra_compile_args = ["-g","-Wall","-O3","-DSEQ_BUFFER_SIZE=5000000"],
                    #extra_compile_args = ["-O0","-fno-inline","-Wall","-g","-UNDEBUG","-DSEQ_BUFFER_SIZE=5000000" ],
                    extra_link_args = [])



alignCompileArgs = ["-Wall","-g","-UNDEBUG","-DEXTRADEBUG"]
#alignCompileArgs = ["-fno-inline","-Wall","-g","-UNDEBUG"]
#alignCompileArgs = ["-g","-Wall","-O3","-DHAVE_GZSTREAM=1"]

alignLibDirs=[]
alignLibs=["stdc++"]
try:
    import gzip
    print "Looks like you have zlib! (It's a bad thing if you don't)"
    alignCompileArgs.extend(["-I./gzstream","-DHAVE_GZSTREAM=1"])
    alignLibDirs.append("src/gzstream")
    alignLibs.extend(["gzstream","z"])
except ImportError:
    print "Looks like you don't have zlib! You will not be able to use gzip:ed files in alignment"
    pass


modAlign = Extension('align',
                     library_dirs = alignLibDirs,
                     libraries = alignLibs,
                     sources = ['src/align.cc'],
                     extra_compile_args=alignCompileArgs,
                     extra_link_args = [])

modDist = Extension("editdist",
                    sources = ["src/editdist.c"],
                    extra_compile_args = ["-Wall","-O3"]
                    )



setup (name = 'mabs',
       version = '1_beta10',
       url = "http://www.cs.helsinki.fi/u/kpalin/",
       author = "Kimmo Palin, Matthias Berg",
       author_email = "kimmo.palin@helsinki.fi",
       maintainer = "Kimmo Palin",
       maintainer_email = "kimmo.palin@helsinki.fi",
       license = "GPL (see file COPYING)",
       description = 'c++ extension modules:\na binding site matrix\nan alignment function',
       ext_modules = [modMatrix, modAlign,modDist],
       packages = ["mabslib"],
       scripts = [ "mabs"] )
