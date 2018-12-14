# askap-image-diagnostic
A module to perform diagnostic analysis on ASKAP images using the SUMSS catalogue. Creates a crossmatch catalogue for SUMSS -> ASKAP sources and also produces diagnostic plots. Also includes the ability to create postage stamps of each crossmatch.

## Installation
I recommend to install the module to a new python environment using, for example, conda or virtualenv.

**Note** This is written for python 2.7, it is not currently compatible with python 3.

To install using pip:

`pip install git+https://github.com/ajstewart/askap-image-diagnostic.git`

Or you can clone the git repository and install using `python setup.py install` or `pip install .`.

## Usage
The built in processing script for the module, available from the command line, is `processASKAPimage.py`.

By default, which means no askap or sumss csv files are provided, `aegean` will be run on the ASKAP image to extract a source catalogue and the SUMSS catalogue will be automatically fetched from Vizier. The SUMSS catalogue will be trimmed to only those sources that fall within the image area.

More than one image can be passed through the processing script at once - however currently the manual csv inputs do not support multiple entires. Hence let the script automatically do source finding and SUMSS fetching.

A range of options exist to influence processing:

```
usage: processASKAPimage.py [-h] [--output-tag OUTPUT_TAG]
                            [--log-level {WARNING,INFO,DEBUG}] [--clobber]
                            [--askap-csv ASKAP_CSV] [--sumss-csv SUMSS_CSV]
                            [--askap-csv-format {aegean}] [--remove-extended]
                            [--askap-ext-thresh ASKAP_EXT_THRESH]
                            [--sumss-ext-thresh SUMSS_EXT_THRESH]
                            [--use-all-fits] [--write-ann]
                            [--boundary-value {nan,zero}]
                            [--crossmatch-base CROSSMATCH_BASE]
                            [--max-separation MAX_SEPARATION]
                            [--create-postage-stamps]
                            [--sumss-mosaic-dir SUMSS_MOSAIC_DIR]
                            [--aegean-settings-config AEGEAN_SETTINGS_CONFIG]
                            images [images ...]

positional arguments:
  images                Define the images to process

optional arguments:
  -h, --help            show this help message and exit
  --output-tag OUTPUT_TAG
                        Add a tag to the output name. (default: )
  --log-level {WARNING,INFO,DEBUG}
                        Set the logging level. (default: INFO)
  --clobber             Overwrite output if already exists. (default: False)
  --askap-csv ASKAP_CSV
                        Manually define a aegean csv file containing the
                        extracted sources to use for the ASKAP image.
                        (default: None)
  --sumss-csv SUMSS_CSV
                        Manually provide the SUMSS catalog csv. (default:
                        None)
  --askap-csv-format {aegean}
                        Define which source finder provided the ASKAP catalog
                        (currently only supports aegean). (default: aegean)
  --remove-extended     Remove perceived extended sources from the catalogues.
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
  --use-all-fits        Use all the fits from Aegean ignoring all flags.
                        Default only those with flag '0' are used. (default:
                        False)
  --write-ann           Create kvis annotation files of the catalogues.
                        (default: False)
  --boundary-value {nan,zero}
                        Define whether the out-of-bounds value in the ASKAP
                        FITS is 'nan' or 'zero'. (default: nan)
  --crossmatch-base CROSSMATCH_BASE
                        Define the base catalogue in the cross matching
                        (currently not supported). (default: sumss)
  --max-separation MAX_SEPARATION
                        Maximum crossmatch distance (in arcsec) to be
                        consdiered when creating plots. (default: None)
  --create-postage-stamps
                        Produce postage stamp plots of the cross matched
                        sources within the max separation. (default: False)
  --sumss-mosaic-dir SUMSS_MOSAIC_DIR
                        Directory containing the SUMSS survey mosaic image
                        files. (default: None)
  --aegean-settings-config AEGEAN_SETTINGS_CONFIG
                        Select a config file containing the Aegean settings to
                        be used (instead of defaults if none provided).
                        (default: None)

````

## Example
**Input**: An ASKAP image called `image.askap.mosaic.restored.fits`. The pixels outside of the image area are NaNs.

**Want**: To crossmatch the ASKAP image with SUMSS and use only matches that are <= 20 arcsec to perform the analysis. Allow the script to automatically build the catalogues but extended sources need to be removed. In this case, these are defined as sources that have one axis that is 1.4 X larger than the associated beam size axis. Also want to produce postage stamp images of the crossmatches along with producing kvis annotation files.

**Command**:
```
processASKAPimage.py image.askap.mosaic.restored.fits --remove-extended --askap-ext-thresh 1.4 --sumss-ext-thresh 1.4 --max-separation 20.0 --create-postage-stamps --sumss-mosaic-dir /path/to/sumss_mosaics_dir --write-ann
````

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
