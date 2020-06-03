# Replication Code for "Fallout of Lead over Paris from the 2019 Notre-Dame Cathedral Fire"

## Data:
``Soil2.csv``: the location and values fo the collected soil samples. 

## Map of the sampled locations:
[google map](https://www.google.com/maps/d/u/0/viewer?ll=48.85633192823713%2C2.339671533651302&z=14&mid=12dQJIcE-QTZzrH5eCsxET9-NBrvDFmJk)


## Code:
``soilPb.R``;
``type.stan``:
We fit a Gaussian process regression model that takes into account both the spatial distribution and the type of soils.