% a script to get the compiled C++ code for spgL1, useful after a fresh git pull

p = dir('*getmex.m');
t = cd;
addpath(genpath(p.folder))
cd(p.folder)
cd('spgl1-2.1')
spgsetup
cd(t)

