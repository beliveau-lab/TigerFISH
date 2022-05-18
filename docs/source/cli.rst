Command Line Interface
######################

`Tigerfish` was written by Robin Aguilar while in the Department of Genome Sciences at the University of Washington.

Official code website:

More documentation and examples: 

The `Tigerfish` workflow is managed by Snakemake which is used to scale and automatically deploy pipeline jobs in parallel as described in our getting started page.

Essentially, each script in the Tigerfish workflow has a set of parameters that may be modified by users. Each script in the main workflow implements these parameters which are called from the pipeline's config.yml file. As a user, one would only need to modify arguments from the config.yml file in order to get Tigerfish working. But to be thorough, named arguments for each script are defined in detail below. 

Main workflow
#############

Here is a simplified DAG of how Snakemake implements scripts in the Tigerfish workflow:


Repeat_ID
---------

The scripts below are presented in the order that they are executed by `Tigerfish`.
 
.. code-block:: bash

    usage: repeat_ID.py [-h] -j JF_COUNT -i INDEX_FILE -chr CHR_NAME -st START
                    [-w WINDOW_LENGTH] [-t THRESHOLD] [-c COMPOSITION_SCORE]
                    -o_b BED_FILE -m MER_LENGTH


