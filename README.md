#Introduction
COMBAT is a comprehensive data pipeline for converting raw diffuse scattering images into 3D reciprocal space maps. The program is built on top of the diffuse imaging software LUNUS (http://lunus.sourceforge.net/), as well as the Computational Crystallography Toolbox (http://cctbx.sourceforge.net/).



#Install
In order to run these scripts, you'll need a working copy of LUNUS along with the most recent PHENIX software distribution.



#Contents

combat.params:                COMBAT parameter file

combat.py:                    Wrapper script called on command line

combat_aniso.py:              Alternative wrapper script with reversed anisotropic map calculation ordering

genlat_labelit.py:            AVB version of map calculation script (called by combat.py)

genlat_labelit_parallel.py:   MEW version of map calculation script called by combat.py

proc.makeref:                 Shell script for creating reference diffuse frame for processed image scaling

proc.scale:                   Shell script for scaling non-reference processed diffuse images

two_run:                      SGE submission script for UCSF QB3 computational cluster



#Instructions
Because COMBAT is still a work in progress (and unpublished), we don't encourage general usage of the pipeline just yet. However, we're working hard on getting there!



#Authors
Andrew Van Benschoten (andrew.vanbenschoten@ucsf.edu)

Michael Wall          (mewall@lanl.gov)

Nicholas Sauter       (nksauter@lbl.gov)