from distutils.core import setup, Extension

module1 = Extension('matrix',
                    libraries = ["stdc++"],
                    sources = ['src/matrix.cc'],
                    extra_compile_args = ["-g","-Wall","-O3","-DSEQ_BUFFER_SIZE=5000000"],
                    #extra_compile_args = ["-O0","-fno-inline","-Wall","-g","-UNDEBUG","-DSEQ_BUFFER_SIZE=5000000" ],
                    extra_link_args = [])

module2 = Extension('align',
                    libraries = ["stdc++"],
                    sources = ['src/align.cc'],
                    #extra_compile_args = ["-Wall","-g","-UNDEBUG","-DEXTRADEBUG"],
                    extra_compile_args = ["-fno-inline","-Wall","-g","-UNDEBUG"],
                    #extra_compile_args = ["-g","-Wall","-O3"],
                    extra_link_args = [])
module3 = Extension("editdist",
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
       ext_modules = [module1, module2,module3],
       packages = ["mabslib"],
       scripts = [ "mabs"] )
