# NOTICE

This is an _unofficial_ repository for the trunk branch of the GENIE Neutrino 
Event MC generator. Please see the [official webpage](http://genie.hepforge.org)
for code that is appropriate for publications-level physics analysis. 

This branch is current with trunk as of SVN commit 7219.

## Versioning and Releases

The complete list of available versions of GENIE can be seen on the official
repository using the Subversion command:

    svn list http://genie.hepforge.org/svn/generator/branches

_This_ repository contains different versions of GENIE stored on branches.
You can see the list of available branches with

    git remote show origin

You can check out a branch directly with, for example, 

    git checkout trunk

to get the `trunk` branch.

## About GENIE

GENIE (Generates Events for Neutrino Interaction Experiments) is a universal 
object-oriented neutrino MC generator supported and developed by an [international 
collaboration](http://genie.hepforge.org/collaboration.html) of scientists whose 
expertise covers a very broad range of neutrino  physics aspects, both 
phenomenological and experimental. GENIE is currently being 
used by T2K, NOvA, MINERvA, MicroBooNE, ArgoNEUT, LAGUNA-LBNO, DUNE, INO, IceCUBE, 
NESSiE and others. 

GENIE uses several third-party HEP codes:

* [ROOT](http://root.cern.ch)
* [LHAPDF](https://lhapdf.hepforge.org)
* [Pythia](http://home.thep.lu.se/~torbjorn/Pythia.html)

GENIE also uses [log4pp](http://log4cpp.sourceforge.net) and the [Gnu 
Scientific Library](http://www.gnu.org/software/gsl/)
