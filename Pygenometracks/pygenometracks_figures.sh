#!/bin/bash

#this code is used for the pygenometracks figures
#CD1A
pyGenomeTracks --tracks ./tracks.ini --region chr1:158000000-158600000 -o ../CD1A.pdf

#MEF2C
pyGenomeTracks --tracks ./tracks.ini --region chr5:88300000-89200000 -o ../MEF2C.pdf

#CD56 (NCAM1)
pyGenomeTracks --tracks ./tracks.ini --region chr11:112800000-113800000 -o ../CD56.pdf