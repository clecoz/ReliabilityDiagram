# Documentation  

Project description --> describe what the project does

_References_:

[1] 

[2] 

[3] 

## _Installation:_

```sh
pip install PackageName
```

## _Parameters:_
observation,forecast,climatology,event_lbound,event_ubound,closed_ends='both',nbins=5,weights=None

**observation**: numpy.ndarray
Timeseries of observation/truth (1D).

**forecast**: numpy.ndarray
Forecast array (2D) with first dimension equal to the dimension of observation, and second dimension equal to the ensemble size.
    
**climatology**: numpy.ndarray
Climatological array (2D) with first dimension equal to the dimension of observation, and second dimension equal to the number of years considered in climatology.

**event_lbound**: int or float
The lower bound for the event considered. 
E.g. 1: If the event considered is lower tercile, then event_lbound value is 0. 
E.g. 2: If the event considered is upper tercile, then event_lbound value is 2/3. 
_NOTE_: The value should be in the range of 0 to 1.

**event_ubound**: int or float
The upper bound for the event considered.
E.g. 1: If the event considered is lower tercile, then event_ubound value is 1/3.
E.g. 2: If the event considered is upper tercile, then event_ubound value is 1. 
_NOTE_: The value should be in the range of 0 to 1.
            
**closed_ends**: str, optional
The bounds (upper and/or lower) to include in the event formulation. Options: 'left', 'right', 'none', 'both'. 
E.g. 1: If closed_ends = 'left', then the event becomes event_lbound <= event < event_ubound. 
E.g. 2: If closed_ends = 'both', then the event becomes event_lbound <= event <= event_ubound.
_NOTE_: The default value is 'both'. The value is case sensitive.
            
**nbins**: int, optional
Number of bins to stratify the forecasts into. The default is 5.
_NOTE_: The bins are of equal width.
            
**weights**: numpy.ndarray, optional
The weights for the forecast data. It should have the same shape as the forecast data.


## _Method(s):_

**contingency_table()**:

_Returns_:


**confidence_intervals()**:

_Returns_:


**forecast_attributes()**:

_Returns_:


**plot_diagram()**:

_Returns_:


## _Demonstration:_

```sh
import numpy as np
import PackageName as something
```

Example - 1:
```sh
In [1]: 
Out[1]: 
```

Example - 2:
```sh
In [2]: 
In [3]: 
Out[3]: 
In [4]: 
Out[4]: 
In [5]: 
Out[5]: 
```

