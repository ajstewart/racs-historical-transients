# askap-image-diagnostic
A module to perform diagnostic analysis on ASKAP images using the SUMSS catalogue. Creates a crossmatch catalogue for SUMSS -> ASKAP sources and also produces diagnostic plots. Also includes the ability to create postage stamps of each crossmatch and to search for transients from either catalogue.

## Installation
I recommend to install the module to a new python environment using, for example, conda or virtualenv.

**Note** This is written for python 2.7, it is not currently compatible with python 3.

To install using pip:

```
pip install git+https://github.com/ajstewart/askap-image-diagnostic.git
```

Or you can clone the git repository and install using ```python setup.py install``` or ```pip install .```.

### Creating a Database
This script was intended to be run on the ada machine which has a installation of postgresql available. To create a database run:

```createdb <db name>``` e.g. ```createdb RACS``` (if you get a denied message contact the system administrator)

This will create an empty database with the chosen name. Make sure to note down the database settings (port, user, name) for use with the pipeline options. If the pipeline is run without first initilising the tables then the tables will be newly created. The easiest way to initilise the tables is by setting up the [website](#Installation of the Website).

### Installation of the Website
Included in the repository is `web_server` which is a basic website built using Django to allow the user to explore the results in a convienient way and for other users to give feedback on the crossmatching.

To install, copy the `web_server` directory to a location where you wish to host the website from and `cd` into the web_server directory.

From here rename the `web_server/settings.py.template` to `web_server/settings.py` and edit the file with the correct database information as above.

Now run the migrations as so, this will essentially create the tables in the database:

```
python manage.py makemigrations
python manage.py migrate
```

Now make a `media` direcoty in the `static` directory:

```
mkdir static/media
```

This is where all the postage stamp and other plots are stored by the pipeline. Also make a note of this directory to use in the pipeline (option `--website-media-dir`).

Now the server can be launched (in the example below port 8005 is used):

```
python manage.py runserver 0.0.0.0:8005
```

## Usage
The built pipeline script, available from the command line, is `processASKAPimage.py`.

By default, which means no askap or sumss csv files are provided, `aegean` will be run on the ASKAP image to extract a source catalogue and the SUMSS catalogue will be automatically fetched from Vizier. The SUMSS catalogue will be trimmed to only those sources that fall within the image area.

More than one image can be passed through the processing script at once - however currently the manual csv inputs do not support multiple entires. Hence let the script automatically do source finding and SUMSS fetching if you want to run more than one image through.

A range of options exist to influence processing:

```
usage: processASKAPimage.py [-h] [-c FILE] [--output-tag OUTPUT_TAG]
                            [--log-level {WARNING,INFO,DEBUG}] [--nice NICE]
                            [--clobber CLOBBER] [--weight-crop WEIGHT_CROP]
                            [--weight-crop-value WEIGHT_CROP_VALUE]
                            [--weight-crop-image WEIGHT_CROP_IMAGE]
                            [--convolve CONVOLVE]
                            [--convolved-image CONVOLVED_IMAGE]
                            [--convolved-non-conv-askap-csv CONVOLVED_NON_CONV_ASKAP_CSV]
                            [--sourcefinder {aegean,pybdsf,selavy}]
                            [--frequency FREQUENCY] [--askap-csv ASKAP_CSV]
                            [--sumss-csv SUMSS_CSV] [--nvss-csv NVSS_CSV]
                            [--askap-csv-format {aegean}]
                            [--remove-extended REMOVE_EXTENDED]
                            [--askap-ext-thresh ASKAP_EXT_THRESH]
                            [--sumss-ext-thresh SUMSS_EXT_THRESH]
                            [--nvss-ext-thresh NVSS_EXT_THRESH]
                            [--use-all-fits USE_ALL_FITS]
                            [--write-ann WRITE_ANN]
                            [--boundary-value {nan,zero}]
                            [--diagnostic-max-separation DIAGNOSTIC_MAX_SEPARATION]
                            [--transient-max-separation TRANSIENT_MAX_SEPARATION]
                            [--postage-stamps POSTAGE_STAMPS]
                            [--postage-stamp-selection {all,good,bad,transients}]
                            [--sumss-mosaic-dir SUMSS_MOSAIC_DIR]
                            [--nvss-mosaic-dir NVSS_MOSAIC_DIR]
                            [--aegean-settings-config AEGEAN_SETTINGS_CONFIG]
                            [--pybdsf-settings-config PYBDSF_SETTINGS_CONFIG]
                            [--selavy-settings-config SELAVY_SETTINGS_CONFIG]
                            [--transients TRANSIENTS]
                            [--transients-askap-sumss-snr-thresh TRANSIENTS_ASKAP_SUMSS_SNR_THRESH]
                            [--transients-large-flux-ratio-thresh TRANSIENTS_LARGE_FLUX_RATIO_THRESH]
                            [--db-engine DB_ENGINE]
                            [--db-username DB_USERNAME] [--db-host DB_HOST]
                            [--db-port DB_PORT] [--db-database DB_DATABASE]
                            [--db-tag DB_TAG]
                            [--website-media-dir WEBSITE_MEDIA_DIR]
                            images [images ...]

positional arguments:
  images                Define the images to process

optional arguments:
  -h, --help            show this help message and exit
  -c FILE, --conf_file FILE
                        Specify config file (default: None)
  --output-tag OUTPUT_TAG
                        Add a tag to the output name. (default: )
  --log-level {WARNING,INFO,DEBUG}
                        Set the logging level. (default: INFO)
  --nice NICE           Set the 'nice' level of processes. (default: 10)
  --clobber CLOBBER     Overwrite output if already exists. (default: False)
  --weight-crop WEIGHT_CROP
                        Crop image using the weights image. (default: False)
  --weight-crop-value WEIGHT_CROP_VALUE
                        Define the minimum normalised value from the weights
                        image to crop to. (default: 0.04)
  --weight-crop-image WEIGHT_CROP_IMAGE
                        Define the weights image to use. (default:
                        weights.fits)
  --convolve CONVOLVE   Convolve the image using CASA to SUMSS resolution for
                        crossmatching. (default: False)
  --convolved-image CONVOLVED_IMAGE
                        Define a convolved image that has already been
                        produced. (default: convolved.fits)
  --convolved-non-conv-askap-csv CONVOLVED_NON_CONV_ASKAP_CSV
                        Define the unconvolved catalogue to use when using
                        convolved mode, otherwise it will be generated
                        automatically (if aegaen or pybdsf) (default: None)
  --sourcefinder {aegean,pybdsf,selavy}
                        Select which sourcefinder to use (default: aegean)
  --frequency FREQUENCY
                        Provide the frequency of the image in Hz. Use if
                        'RESTFRQ' is not in the header (default: 864000000.0)
  --askap-csv ASKAP_CSV
                        Manually define a aegean csv file containing the
                        extracted sources to use for the ASKAP image.
                        (default: None)
  --sumss-csv SUMSS_CSV
                        Manually provide the SUMSS catalog csv. (default:
                        None)
  --nvss-csv NVSS_CSV   Manually provide the NVSS catalog csv. (default: None)
  --askap-csv-format {aegean}
                        Define which source finder provided the ASKAP catalog
                        (currently only supports aegean). (default: aegean)
  --remove-extended REMOVE_EXTENDED
                        Remove perceived extended sources from the catalogues.
                        Uses the following arguments 'askap-ext-thresh' and
                        'sumss-ext-thresh' to set the threshold. (default:
                        False)
  --askap-ext-thresh ASKAP_EXT_THRESH
                        Define the maximum scaling threshold of the size of
                        the ASKAP source compared to the PSF. Used to exclude
                        extended sources. Only 1 axis has to exceed. (default:
                        1.2)
  --sumss-ext-thresh SUMSS_EXT_THRESH
                        Define the maximum scaling threshold of the size of
                        the SUMSS source compared to the PSF. Use to exclude
                        extended sources. Only 1 axis has to exceed. (default:
                        1.2)
  --nvss-ext-thresh NVSS_EXT_THRESH
                        Define the maximum scaling threshold of the size of
                        the NVSS source compared to the PSF. Use to exclude
                        extended sources. Only 1 axis has to exceed. (default:
                        1.2)
  --use-all-fits USE_ALL_FITS
                        Use all the fits from Aegean ignoring all flags.
                        Default only those with flag '0' are used. (default:
                        False)
  --write-ann WRITE_ANN
                        Create kvis annotation files of the catalogues.
                        (default: False)
  --boundary-value {nan,zero}
                        Define whether the out-of-bounds value in the ASKAP
                        FITS is 'nan' or 'zero'. (default: nan)
  --diagnostic-max-separation DIAGNOSTIC_MAX_SEPARATION
                        Maximum crossmatch distance (in arcsec) to be
                        consdiered when creating the diagnostic plots.
                        (default: 5.0)
  --transient-max-separation TRANSIENT_MAX_SEPARATION
                        Maximum crossmatch distance (in arcsec) to be
                        consdiered when searching for transients. (default:
                        45.0)
  --postage-stamps POSTAGE_STAMPS
                        Produce postage stamp plots of the cross matched
                        sources within the max separation. (default: False)
  --postage-stamp-selection {all,good,bad,transients}
                        Select which postage stamps to create. (default:
                        ['all'])
  --sumss-mosaic-dir SUMSS_MOSAIC_DIR
                        Directory containing the SUMSS survey mosaic image
                        files. (default: None)
  --nvss-mosaic-dir NVSS_MOSAIC_DIR
                        Directory containing the NVSS survey mosaic image
                        files. (default: None)
  --aegean-settings-config AEGEAN_SETTINGS_CONFIG
                        Select a config file containing the Aegean settings to
                        be used (instead of defaults if none provided).
                        (default: None)
  --pybdsf-settings-config PYBDSF_SETTINGS_CONFIG
                        Select a config file containing the PyBDSF settings to
                        be used (instead of defaults if none provided).
                        (default: None)
  --selavy-settings-config SELAVY_SETTINGS_CONFIG
                        Select a config file containing the Selavy settings to
                        be used (instead of defaults if none provided).
                        (default: None)
  --transients TRANSIENTS
                        Perform a transient search analysis using the
                        crossmatch data. Requires '--max-separation' to be
                        defined. (default: False)
  --transients-askap-sumss-snr-thresh TRANSIENTS_ASKAP_SUMSS_SNR_THRESH
                        Define the threshold for which ASKAP sources are
                        considered to not have a SUMSS match baseed upon the
                        estimated SUMSS SNR if the source was placed in the
                        SUMSS image. (default: 5.0)
  --transients-large-flux-ratio-thresh TRANSIENTS_LARGE_FLUX_RATIO_THRESH
                        Define the threshold for which sources are considered
                        to have a large flux ratio. Median value +/- threshold
                        x std. (default: 3.0)
  --db-engine DB_ENGINE
                        Define the database engine. (default: postgresql)
  --db-username DB_USERNAME
                        Define the username to use for the database (default:
                        postgres)
  --db-host DB_HOST     Define the host for the databse. (default: localhost)
  --db-port DB_PORT     Define the port for the databse. (default: 5432)
  --db-database DB_DATABASE
                        Define the name of the database. (default: postgres)
  --db-tag DB_TAG       The description field in the databased attached to the
                        image. (default: RACS Analysis)
  --website-media-dir WEBSITE_MEDIA_DIR
                        Copy the image directory directly to the static media
                        directory of the website. (default: none)

````

These options can be entered using a ConfigParser configuration file:

```
[GENERAL]
output_tag=askap_racs_analysis
log_level=INFO
nice=10
clobber=True

[ANALYSIS]
frequency=864e6
weight_crop=True
weight_crop_value=0.04
weight_crop_image=../path/to/weight_cropped_image.fits
convolve=True
convolved_image=../path/to/convolved_image.fits
convolved_non_conv_askap_csv=../path/to/preconvolved_catalog.csv_
sourcefinder=aegean
# aegean_settings_config=None
# pybdsf_settings_config=None
# selavy_settings_config=None
boundary_value=nan

[CATALOGUES]
askap_csv=/path/to/askap_catalog.csv
# sumss_csv=None
# nvss_csv=None
askap_csv_format=aegean
write_ann=True

[CROSSMATCHING]
diagnostic_max_separation=5.0
transient_max_separation=45.0
remove_extended=True
askap_ext_thresh=1.3
sumss_ext_thresh=1.5
nvss_ext_thresh=1.2
use_all_fits=False
postage_stamps=True
# postage_stamp_selection=[all,]
sumss_mosaic_dir=/directory/where/sumss/mosaics/are/kept
nvss_mosaic_dir=/directory/where/sumss/mosaics/are/kept

[TRANSIENTS]
transients=True
transients_askap_sumss_snr_thresh=5.0
transients_large_flux_ratio_thresh=2.0

[DATABASE]
db_engine=postgresql
db_username=user
db_host=localhost
db_port=5432
db_database=RACS
db_tag=Tag to add to pipeline run
website_media_dir=/path/to/the/website/media/dir
```

## Weight Cropping
If the weights.XX.fits file is available and supplied then this can be used to trim the image to remove the edges, leaving just the cleaner part of the image. By default this is set to cut to a value of 0.04 of the maximum weight value.

## Image Convolving
The pipeline is able to convolve the supplied ASKAP image to that of the SUMSS or NVSS resolution. If the image has already been convolved it can supplied to the pipeline using the `convolved_image` argument (note make sure to have 'convolved=True'). In the case of convolving to SUMSS the target beam size follows the 45 x 45 cosec |dec| convention. Cross matching is then done against the convolved image.

If convolving is used then the non-convolved image will also be analysed. The source will be extracted by the pipeline using Aegean or a catalogue can be supplied using the `convolved_non_conv_askap_csv` argument. This is used when searching for transient sources.

## Image Diagnostic Plots
These plots are produced by only using sources that have a crossmatch distance <= the user defined max separation, i.e. good matches. Also if `--remove-extended` is enabled then extended sources are also removed from the list of crossmatches used to produce the diagnostic plots.

**Note** the source numbers plot does not make these exclusions.

## Transient Searching
The pipeline works by matching each SUMSS source in the image with the nearest ASKAP source extracted.

Good matches are deemed those that are <= the max separation defined by the user. Above this is considered to have no match. Using this approach the results are categorised into four categories:

* **No ASKAP Match to SUMSS** - This defines a SUMSS source that has no ASKAP source matched to it within the max separation limit.
* **No SUMSS Match to ASKAP** - This defines an ASKAP source that has not been matched to a SUMSS source within the max separation limit **AND** has an integrated flux density such that it would be at least a 5 sigma detection in the SUMSS image (this does not yet account for spectral index).
* **Large Ratio** - These are sources that are good matches but have a integrated ASKAP/SUMSS flux ratio that is >/< the median flux ratio +/- standard deviation.
* **Good Matches** - The sources that are defined as being a good match (including the large ratio sources).

## Output
In the top level directory will be:

* log files.
* png files of diagnostic plots.
* csv files of askap catalogue, sumss sources and the complete crossmatching result.

Two directories may also be present:

* postage-stamps - this will store all the postage stamp images which will be sorted into good and bad matches.
* transients - this will store all the csv files of categorised transient candidates and in sub-directories will also be copies of the postage stamps if these have been created (renamed for the respective transient source and category).

## Example
**Input**: An ASKAP image called `image.askap.mosaic.restored.fits`. The pixels outside of the image area are NaNs. Our database is on the localhost with the username of `user123`, on the default port and is called `racstest`. Our website media directory is `/my/website/static/media`

**Want**: To crossmatch the ASKAP image with SUMSS and use only matches that are <= 20 arcsec to perform the analysis. Allow the script to automatically build the catalogues and remove extended sources when creating the diagnostic plots. In this case, these are defined as sources that have one axis that is 1.4 X larger than the associated beam size axis. Also want to produce postage stamp images of the crossmatches along with producing kvis annotation files, and finally perform a transient search. We will mark the image in the database with `first test`.

**Command**:
```
processASKAPimage.py image.askap.mosaic.restored.fits --remove-extended --askap-ext-thresh 1.4 --sumss-ext-thresh 1.4 --max-separation 20.0 --postage-stamps --sumss-mosaic-dir /path/to/sumss_mosaics_dir --write-ann --transients --db-username user123 --db-name racstest --db-tag "first test" --website-media-dir /my/website/static/media
```

Or we can use configure the parset file to run:

```
[GENERAL]
output_tag=example
log_level=INFO
nice=10
clobber=True

[ANALYSIS]
frequency=864e6
weight_crop=False
weight_crop_value=0.04
weight_crop_image=../path/to/weight_cropped_image.fits
convolve=False
convolved_image=../path/to/convolved_image.fits
convolved_non_conv_askap_csv=../path/to/preconvolved_catalog.csv_
sourcefinder=aegean
# aegean_settings_config=None
# pybdsf_settings_config=None
# selavy_settings_config=None
boundary_value=nan

[CATALOGUES]
askap_csv=/path/to/askap_catalog.csv
# sumss_csv=None
# nvss_csv=None
askap_csv_format=aegean
write_ann=True

[CROSSMATCHING]
diagnostic_max_separation=5.0
transient_max_separation=20.0
remove_extended=True
askap_ext_thresh=1.4
sumss_ext_thresh=1.4
nvss_ext_thresh=1.2
use_all_fits=False
postage_stamps=True
# postage_stamp_selection=[all,]
sumss_mosaic_dir=/directory/where/sumss/mosaics/are/kept
nvss_mosaic_dir=/directory/where/sumss/mosaics/are/kept

[TRANSIENTS]
transients=True
transients_askap_sumss_snr_thresh=5.0
transients_large_flux_ratio_thresh=2.0

[DATABASE]
db_engine=postgresql
db_username=user123
db_host=localhost
db_port=5432
db_database=racstest
db_tag=first test
website_media_dir=/path/to/the/website/media/dir
```

And then run the pipeline like so:

```
processASKAPimage.py -c myparset.in image.askap.mosaic.restored.fits
```

**Output**: The results will be placed in `image.askap.mosaic.restored_results`.

## Aegean settings
The default aegean settings are:
```
cores=1
maxsummits=5
seedclip=5
floodclip=4
nocov=True
````

These can be changed by providing a config file and supplying it to the argument `--aegean-settings-config`. There should be a standard ConfigParser header `[aegean]`. E.g.:
```
[aegean]
cores=12
maxsummits=5
seedclip=6
floodclip=4
autoload=True
````
To deactivate a setting remove it from the config file.

## PyBDSF settings
Similar to above PyBDSF settings can also be supplied in the same format. These should use the header `[pybdsf]`.
