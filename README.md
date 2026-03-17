This pipeline requires the installation of astrosource as Python environment, with the scripts photometry1_sep.py, photometry2_sep.py and spiders_variability.py ran as command lines from the terminal (e.g. "> spiders_variability.py") within the astrosource environment (#!/home/gudrun/marcotu/python-envs/astrosource/bin/python). The astrosource scripts must be also replaced with the ones uploaded in this gitlab repository at astrosource/ and bin/ directories. The pipeline also requires sep, numpy and astropy libraries (check possible errors for further dependencies). Each code listed below must be place in the "$HOME/bin" folder and ran as command line (e.g. "> sortandmedian.py") in the following order:
#If images coordinates are not aligned with WCS coordinates:

fix_wcs_browser.py

#If astrometry.net is installed on the machine with proper index files:

fix_wcs_local.py


sortandmedian.py

photometry1_sep.py

spiders_variability.py

photometry2_sep.py

spiders_variability.py

fieldsofview.py

variable_selection.py

spiders_lightcurves.py

final_plots_gri_rmethod.py

#select the promising variables and put their indeces inside a text file titled as "final_periodic_can_1.txt"

final_plots_match.py

final_plots_spidershunt2_paper.py
