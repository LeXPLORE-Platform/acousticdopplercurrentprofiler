# LéXPLORE Acoustic Doppler Current Profilers

## Project Information

[LéXPLORE](https://lexplore.info) is a state-of-the-art research platform situated on Lake Geneva, dedicated to a wide range of limnological studies. This innovative platform is the result of a collaborative effort between five leading institutions: [Eawag](https://www.eawag.ch/en/), [EPFL](https://www.epfl.ch/en/), [INRAE](https://www6.lyon-grenoble.inrae.fr/carrtel/), [UNIGE](https://unige.ch), [UNIL](https://www.unil.ch/index.html). Since February 2019, the LéXPLORE platform is anchored at 110 m depth off the lake's north-shore (46°30'00.8"N 6°39'39.0"E).

The data presented here is part of the core dataset maintained by the technical team of LéXPLORE.
The data is used and displayed on the [Datalakes website](https://www.datalakes-eawag.ch/) where other related data or products can be visualised and downloaded.

The specific dataset contains raw and processed data recorded by two acoustic Doppler current profilers (ADCPs): an upward-looking RDI 600 kHz samples near-surface waters and a downward-looking RDI 300 kHz samples deep waters. Data includes backscattering and 3D water velocities (see [data section](#data)).

**References**:

Wüest, A., Bouffard, D., Guillard, J., Ibelings, B. W., Lavanchy, S., Perga, M. ‐E., & Pasche, N. (2021). LéXPLORE: a floating laboratory on Lake Geneva offering unique lake research opportunities. Wiley Interdisciplinary Reviews: Water, 8(5), e1544 (15 pp.). https://doi.org/10.1002/wat2.1544

See also the [360° virtual tour](https://www.eawag.ch/repository/lexplore/index.htm).

## Citation

<span style="color:red">*To be added here once there is a DOI number.*</span>

## Keywords
LéXPLORE, Lake Geneva, Water velocity, Backscattering, ADCP, Datalakes

## Sensors

An acoustic Doppler current profiler (ADCP) is a hydroacoustic current meter similar to a sonar, used to measure water current velocities over a depth range using the Doppler effect of sound waves scattered back from particles within the water column. ADCPs provide simultaneous measurements of backscattering from particles and water velocity. 

The two different ADCPs used to cover the water column near LéXPLORE are cable linked to the platform and record live data. Characteristics and standard setup of each ADCP are given below. Note that the setup can be modified for specific studies. Please get in touch with the science contact if you are interested in a specific configuration.

### Upward-looking RDI 600 kHz (near-surface velocities and backscattering)
- **Brand, Model & SN**:    Teledyne RD Instruments, Workhorse Sentinel 600 kHz, SN 20515 (before 12.01.2023 and after 19.04.2023), SN 1178 (between 19.01.2023 and 19.04.2023)
- **System integration**: VmDAS
- **Accuracy**:             ​±0.3% of measured velocity ± 0.3 cm/s
- **Setup**:                Bin size: 0.25 m; Averaging ensemble interval: 10 minutes; Number of pings per ensemble varies with the deployment period
- **Deployment depth**: varies with the deployment period (8-28 m), see `scripts/parameters.json` for more information

### Downward-looking RDI 300 kHz (deep velocities and backscattering)
- **Brand, Model & SN**:    Teledyne RD Instruments, Workhorse Sentinel 300 kHz, SN 20360 (before 12.01.2023), SN 20847 (after 19.01.2023)
- **System integration**: VmDAS
- **Accuracy**:             ±1.0% of measured velocity ± 0.5 cm/s
- **Setup**:                Bin size: 1 m; Averaging ensemble interval: 10 minutes; Number of pings per ensemble varies with the deployment period
- **Deployment depth**: varies with the deployment period (8-29 m), see `scripts/parameters.json` for more information

More information about Workhorse Sentinel ADCPs can be found in the manuals in the foler `notes`.

## Geospatial Information

The two ADCPs are deployed 35 m north-east of LéXPLORE: 46°30′02.057″N, 6°39′39.741″E.

## Temporal coverage 
- start: July 2020
- stop: live dataset

## Scripts

## Code

[![License: MIT][mit-by-shield]][mit-by]

:warning You need to have [git](https://git-scm.com/downloads) and [git-lfs](https://git-lfs.github.com/) installed in order to successfully clone the repository.

- Clone the repository to your local machine using the command: 

 `git clone  https://gitlab.renkulab.io/lexplore/acousticdopplercurrentprofiler`

 Note that the repository will be copied to your current working directory.

- Use Python 3 and install the requirements with:

 `pip install -r requirements.txt`

 The python version can be checked by running the command `python --version`. In case python is not installed or only an older version of it, it is recommend to install python through the anaconda distribution which can be downloaded [here](https://www.anaconda.com/products/individual). 


### Process new data

In order to process new data locally on your machine the file path needs to be adapted to your local file system. The following steps are therefore necessary: 

- Edit the `scripts/input_batch.bat` file. Change all the directory paths to match your local file system. This file contains all the file paths necessary to launch the batch scripts `runfile.bat`.

- Edit the `scripts/input_python.py` file. Change all the directory paths to match your local file system. This file contains all the directories where the python script outputs data to.

To process new data, place the data in the input directory which you specified in the `scripts/input_batch.bat` file. Update `scripts/parameters.json` if the ADCPs configuration has changed compared to the previous deployment period (e.g., deployment depth, orientation). Double-clicking on the `runfile.bat` file will automatically 
process all the data in the input directory and store the output in the directories specified in the `scripts/input_python.py` file. 

### Adapt/Extend data processing pipeline

The python script `scripts/main.py` defines the different processing steps while the python script `scripts/adcp.py` contains the python class ADCP with all the corresponding 
class methods to process the data. To add a new processing or visualization step, a new class method can be created in the `adcp.py` file and the step can be added in `main.py` file.
Both above mentioned python scripts are independent of the local file system.

In addition, `scripts/functions.py` and `scripts/general_functions.py` contain ADCP-specific and more general functions, respectively. ADCP-specific quality checks are defined as functions in `scripts/quality_checks_adcp.py`. The script `scripts/quality_assurance.py` runs advanced quality checks based on `quality_assurance.json`. The functions `scripts/download_data.py` and `scripts/upload_data.py` are used to download and upload data, respectively, between the local repository and the cloud (see `data/README.md` for more information). 

## Data

### License

[![CC BY 4.0][cc-by-shield]][cc-by] 

This data is released under the Creative Commons license - Attribution - CC BY (https://creativecommons.org/licenses/by/4.0/). This license states that consumers ("Data Users" herein) may distribute, adapt, reuse, remix, and build upon this work, as long as they give appropriate credit, provide a link to the license, and indicate if changes were made.
 
The Data User has an ethical obligation to cite the data source (see the DOI number) in any publication or product that results from its use. Communication, collaboration, or co-authorship (as appropriate) with the creators of this data package is encouraged. 
 
Extensive efforts are made to ensure that online data are accurate and up to date, but the authors will not take responsibility for any errors that may exist in data provided online. Furthermore, the Data User assumes all responsibility for errors in analysis or judgment resulting from use of the data. The Data User is urged to contact the authors of the data if any questions about methodology or results occur. 


### Data Structure

- **Level 0**: Raw data collected from the different sensors.

- **Level 1**: Raw data stored to NetCDF file where attributes (such as sensors used, units, description of data, etc.) are added to the data. Column with quality flags are added to the Level 1A data. Quality flag >1 indicates that the data point did not pass one or several 
quality checks and further investigation is needed (see section [quality assurance](#quality-assurance)), quality flag "0" indicates that no further investigation is needed. Masked data can be found in the L1 product with the extention '_qual'.

- **Level 2**: Smooth data with moving average filter and additional parameters such as mU, mdir and Sv.


**Netcdf file info (Level 2):**
* Coordinates: 
    * depth [m]
    * time [UTC] 
* Data variables:
    * eastward velocity, *u* (depth,time) [m s^(-1)]
    * westward velocity, *v* (depth,time) [m s^(-1)]
    * upward velocity, *w* (depth,time) [m s^(-1)]
    * temperature, *temp* (time) [°C]
    * velocity error, *vel_err* (depth,time) [m s^(-1)]
    * echo amplitude from beam #1, *echo1* (depth,time) [counts]
    * echo amplitude from beam #2, *echo2* (depth,time) [counts]
    * echo amplitude from beam #3, *echo3* (depth,time) [counts]
    * echo amplitude from beam #4, *echo4* (depth,time) [counts]
    * correlation from beam #1, *corr1* (depth,time) [counts]
    * correlation from beam #2, *corr2* (depth,time) [counts]
    * correlation from beam #3, *corr3* (depth,time) [counts]
    * correlation from beam #4, *corr4* (depth,time) [counts]
    * percentage good #1 (percentage of data obtained with 3 beams due to one beam with bad data), *prcnt_gd1* (depth,time) [%]
    * percentage good #2 (percentage of data with lower velocity error than threshold), *prcnt_gd2* (depth,time) [%]
    * percentage good #3 (percentage of data rejected because more than one beam with bad data), *prcnt_gd3* (depth,time) [%]
    * percentage good #4 (percentage of data obtained with 4 beams), *prcnt_gd4* (depth,time) [%]
    * battery level, *battery* (time) [counts]
    * pitch angle (rotation around x-axis going from beam #1 to beam #2), *pitch* (time) [°]
    * roll angle (rotation around y-axis going from beam #3 to beam #4), *roll* (time) [°]
    * heading angle (rotation around z-axis, 0° if beam #3 is facing north), *heading* (time) [°]
    * depth-averaged velocity magnitude, *mU* (time) [m s^(-1)]
    * direction of depth-averaged velocity (anticlockwise from east), *mdir* (time) [°]
    * absolute backscatter, *Sv* (depth,time) [dB]

An example of visualisation of the Netcdf files is provided as a Jupyter Notebook (`notebooks/vis_adcp.ipynb`)

**Reading from NetCDF**.
There are a number of resources that give detailed information on how to read and interact with NetCDF files. Linked below are some suggested resources.

- Python		https://unidata.github.io/netcdf4-python/
- R		https://cran.r-project.org/web/packages/ncdf4/ncdf4.pdf


**Reading Time.** Datetime is in Unix time format (seconds since 01 Sept. 1970 [UTC]). *Warning*: the display in [Datalakes website](https://www.datalakes-eawag.ch/) is in local time but the downloaded data are always in UTC. Most languages have a function for parsing this format to a datetime object.

- Python		
```
from datetime import datetime 
dt = datetime.utcfromtimestamp(unixdatetime)
```
- R	
```
library(anytime)
dt <- anytime(unixdatetime)
```

## Quality assurance

### Basic checks
Quality checks include but are not limited to range validation, data type checking and flagging missing data. The basic quality checks are defined in `quality_assurance.json` and include the following tests (flag index is "1" if the test is not passed, "0" otherwise):
* "time": numeric, bounds: [1514764800, "now"],
* "u": numeric, bounds: [-5, 5],
* "v": numeric, bounds: [-5, 5],
* "w": numeric, bounds: [-5, 5].

### ADCP-specific tests

In addition, the following ADCP-specific tests are performed to increase the flag index of "u", "v" and "w" following a base 2 format (see function `scripts/quality_checks_adcp.py` for the tests):
* interface detection (flag index: 2): data above the surface (upward ADCP) or below the lake bottom (downward ADCP) is flagged based on the transducer depth and bottom depth specified in `scripts/parameters.json`.
* minimum correlation (flag index: 4): data with at least one beam below a correlation threshold is flagged (default threshold: 64 counts).
* minimum good data percentage (flag index: 8): data with PG1+PG4 < threshold is flagged (default threshold: 25 %).
* maximum bad data percentage (flag index: 16): data with PG3 > threshold is flagged (default threshold: 25 %).
* maximum velocity error (flag index: 32): data with a velocity error above a threshold is flagged (default threshold: 0.05 m/s).
* maximum tilt (flag index: 64): data with either pitch or roll angle above a threshold is flagged (default threshold: 15°).
* maximum 4-beams correlation difference (flag index: 128): data with differences in correlation ratio (0-1) between the 4 beams above a threshold is flagged (default standard deviation threshold: 0.01).
* minimum vertical echo difference (flag index: 256): data with at least one beam with vertical echo difference between two consecutive bins above a threshold is flagged (default threshold: 30 counts).

The flag indices can be found in the data variables "u_qual", "v_qual" and "w_qual" of the L1 and L2 products.

### Advanced quality assurance

Advanced quality assurance can be run using the `scripts/quality_assurance.py` function. 

###  Events 

Maintenance dates, interesting or surprising events, non identified by the basic quality check are listed in `notes/events.csv`.
Check also `notes/sensor_history.csv`.

## Collaborators

- **Concept, finances, project management** Damien Bouffard, Jean Guillard, Bas Ibelings, Natacha Pasche, Marie-Elodie Perga, Alfred Wüest   
- **Installation, maintenance, data collection** Guillaume Cunillera, Matteo Gios, Roxane Fillon, Jeremy Keller, Sébastien Lavanchy, Floreana Miesen, Michael Plüss, Philippe Quetin
- **Data pipeline** Damien Bouffard, Tomy Doda, James Runnalls
- **Data review** Damien Bouffard, Tomy Doda

## Contact
- **Contact science** [Damien Bouffard](mailto:damien.bouffard@eawag.ch), [Tomy Doda](mailto:tomy.doda@unil.ch)
- **Contact software** [James Runnalls](mailto:james.runnalls@eawag.ch)
- **Contact tech** [Guillaume Cunillera](mailto:guillaume.cunillera@epfl.ch)


[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-g.svg?label=Data%20License
[mit-by]: https://opensource.org/licenses/MIT
[mit-by-shield]: https://img.shields.io/badge/License-MIT-g.svg?label=Code%20License
