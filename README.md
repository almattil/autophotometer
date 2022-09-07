# Autophotometer

Autophotometer is a semi-automatic Python based pipeline for photometric measurements.

It uses Source Extractor and SCAMP to find the detections from the user's FITS images, 
and finds matches for them from an online database. 
Autophotometer's objective is to provide a quick and 
easy way to study variable and transient stars, 
but can be used for other types of objects as well.

### Prerequisites
Autophotometer requires [SExtractor](https://github.com/astromatic/sextractor) version 2.25.0 or newer,
[SCAMP](https://github.com/astromatic/scamp) version 2.10.0 or newer, and their software requirements. 
Additionally it uses the following python packages: 
  - Astropy
  - Numpy
  - Matplotlib
  - Scipy
  - Configparser
  - photutils
  - CCDProc
  - Plotille
  - Requests

Autophotometer was developed using Python 3.8.

### Installing - Using pip

Install the program using command

    pip install .
    
in the directory of the downloaded files. 
The use of flag --user is recommended, so the program is installed for only the current user. 
It should be noted that when installing in a virtual environment, --user is not necessary.

## Usage

Autophotometer currently has two modes: single image and dual image operations. In single mode the program processes only one image, and in dual mode it runs two images of the field in different filter bands. Dual image mode compares the images to each other and outputs values for the detections in both images.

The dual mode is more constricted, as it only uses detections present in both images, which can cause some issues. The dual image functionality may be removed in future, as the single image mode gives more reliable results.

Running the program is done using the terminal. It can be called from any directory using the command:

    AUTOPHOTOMETER [images]

The argument "images" can contain up to two FITS files, depending on which mode you want to use. The program finds the used filter bands from the fits headers, and will give an error if a filter is unsupported. The fits header must contain keywords that describe the coordinate system, coordinate system reference pixels, and other keywords for SCAMP to be able to process them.

Autophotometer then asks if the user wants to make modifications to the fits files. It can create a backup of the original file, remove cosmic rays, and do the coordinate correction.



After the modifications, the program shows two histograms. The first histogram shows all the FWHM values of the detections in arcseconds. In the second histogram only the detections within a threshold are kept. This is done to exclude detections with too big or small FWHM, which might be galaxies, clusters, cosmic rays, or other non-star objects. The detections that are left are most likely stars.

Then the program shows a figure of the colour term calculations. In the figure the X-axis shows the difference in magnitudes between two adjacent filters and the Y-axis shows the uncorrected instrumental magnitudes. A linear fit to the data is used to find the extinction coefficient and zero point magnitude. The figure can be used to evaluate whether the fit is good or not.

Finally the program draws the fits image, with markers depicting the different detection types and if the detections are present in the database. Note that in the figure, some detections are left unmarked. The detection threshold can be adjusted in the conf\_autophotometer.ini file. 

When a marker is clicked, the program outputs the coordinates of the detected source in degrees, magnitudes, database magnitudes, and errors in the terminal. The program also saves the information of the selected detections and writes it into a file called "Selected\_Mags\_[image name].dat," which is generated in the directory where the program is executed. 
