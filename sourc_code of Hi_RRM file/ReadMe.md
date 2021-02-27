We only modify the mvlmm.cpp file of GEMMA software. 

After confirming the installation of GSL, if there is a problem with the dynamic library, you could try to put my dynamic library file into "/usr/lib64" or "/usr/lib", or make a copy of "libgsl.so.23" then rename it to "libgsl.so.25".

You also could download GEMMA source code(https://github.com/genetics-statistics/GEMMA) and replace this file "./GEMMA/src/mvlmm.cpp" by "mvlmm.cpp", then recompile the Hi_RRM file.
