# LAOS Fitting
This program fits an oscillatory shear rheology dataset containing small and large amplitudes to the constitutive model developed by Ruud van der Sman. The research paper is available [here](https://doi.org/10.1016/j.foodhyd.2023.109586).  


## Data Reading
- data_reading.ipynb : the script to read the specified data file in the ./experimental data
- ./data_reading : is teh directory where teh data read by the data_reading.ipynb is stored.
- data_reading_inputs.py : input file for data_reading.ipynb, it is automatically read in the data_reading.ipynb.
 
## Data Fitting
fitting.ipynb :  the fitting script, whivh takes fitting_inputs.py as an input. 
