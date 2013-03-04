
DOXYGEN documentation
=====================

Just install doxygen and run:
$ doxygen doxyfile

There will be a new folder html/index.html

If you want callgraphs, install the tool dot.

You can specify source folders in doxyfile line 613
INPUT                  = src/wwmII


Note: Remove all *genmod*.f90 file from the source folder befor you run doxygen.
      The Intel Fortran Compiler has a command line parameter -module modules to put
      all the trash files into a modules/ folder.
Note: Only the wwm_petsc* files a commented with doxygen tags at the moment.
