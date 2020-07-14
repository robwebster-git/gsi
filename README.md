# Trace Species Identifier

A python tool written to address _*GSI 654*_

The tool currently has two modes with different input data.

## Method 1 - Using percentage adjusted species rasters as inputs

*inputs* : pc_adj (percentage adjusted) species rasters.  Trace species are identified as those which contribute an overall percentage in the valid (ie not masked) pixels in the AOI of less than a threshold percentage (default 1%).

*masks* : FNF raster and/or shapefile with property boundaries or similar

*outputs* : a text file which includes the filenames of those input files pertaining to trace species.  Second output, if flag set, is a new raster with pixel values equal to the summed contributions of each species identified overall as trace in the AOI.

## Method 2 - Using a dominant species raster as input

*inputs* : a dominant species raster.  Trace species are identified as those which are dominant in less than a threshold percentage (default 1%) of valid (ie not masked) pixels.

*masks* : FNF raster and/or shapefile with property boundaries or similar

*output* : currently just a textfile with details on the species and their contributions, along with which are identified as trace species.  These species are in the form of FIA codes, so depending on what needs to be done with this information, they may well have to be converted to USFS common name “codes” or similar.

## Examples of Usage

Calculate using method 1, providing a GeoJSON file and a threshold of 2%

`python trace_species.py --pc_adj input_files_directory/ --shapefile area_of_interest.geojson --write_trace_raster output.tif --threshold 2`

Calculate using method 2, with a threshold of 5%

`python trace_species.py --dom dominant.tif --threshold 5`

Calculate using method 2, specifying a Forest/Non-Forest (FNF) mask and an area file with which to select the relevant area and mask out invalid pixels properly.  The AOI is specified as "LV_KC", though this is also the default.  A threshold of 4% is used.

`python trace_species.py --dom trace_test_files/dominant_10_nf.tif --aoi LV_KC --fnf trace_test_files/nonforest_all_roads_LV_KC.tif --shapefile trace_test_files/top_of_maine.geojson --threshold 4`

Calculate using method 1, supplying a FNF mask, a shapefile, and an output trace species raster filename.  The output raster will be cropped to the extent of the shapefile features, and will have nodata value of -9999, including for masked areas.  

`python trace_species.py --pc_adj trace_test_files/ --fnf trace_test_files/nonforest_all_roads_LV_KC.tif --shapefile trace_test_files/top_of_maine.geojson --write_trace_raster trace_species_output3.tif`

## Help

```bash
python trace_species.py --help             Usage: trace_species.py [OPTIONS]

Options:
  --pc_adj PATH              Path to directory containing percentage adjusted
                             species rasters.  Use EITHER this method OR the
                             --dom method below.

  --dom PATH                 Specifies that you want to calculate using a
                             dominant species raster.  Provide such a raster
                             file here

  --out PATH                 Specify a directory path for the output files.
                             The default is the current directory

  --aoi TEXT                 The string code for the area of interest.
                             Default is LV_KC

  --threshold INTEGER        Threshold for calculating trace species.  Default
                             is 1%

  --fnf TEXT                 Add a forest/non-forest mask - specify the
                             filepath

  --shapefile TEXT           Calculate trace species within this area -
                             provide path to shapefile

  --write_trace_raster TEXT  Write output trace species raster - provide
                             output filename

  --help                     Show this message and exit.
```
