import rasterio as rio
import rasterio.mask
import rasterio.plot
import numpy as np
import numpy.ma as ma
from glob import glob
import os
import click
import sys
import fiona

"""
Trace Species ID, version 1.1
Rob Webster - July 2020
Ticket Number - GSI 654
"""


def find_trace_by_pc_adj(files, fnf_exists, fnf_array, threshold, shapes):
    """
    Takes in a list of species rasters (percentage adjusted), an optional FNF (forest/non-forest mask) and an optional
    shapefile if you only want to analyse the area within certain polygons.  The threshold is the percentage contribution of a particular
    species, below which the species is said to be "trace".
    
    Returns:
    
    trace_species : list - filenames of trace species
    trace_data : a list of masked arrays of the data from the species identified as trace
    profile : return metadata for writing the trace species data to a new raster
    process_info : list - textual info about the process
    """
    trace_species = []
    trace_data = []
    process_info = []

    profile = rio.open(files[0]).profile.copy()

    # Loop through each adjusted species raster
    for file in files:
        
        name = file.split('_pc_adj')[0]

        with rio.open(file) as f:

            #  If there's a shapefile supplied, crop each image to the shapefile and update metadata

            if shapes:
                cropped_image, cropped_transform = rasterio.mask.mask(f, shapes, crop=True, nodata=-9999, filled=True)
                cropped_meta = f.meta
                profile['height'] = cropped_image.shape[1]
                profile['width'] = cropped_image.shape[2]
                profile['transform'] = cropped_transform
                
                # If a FNF mask is supplied, apply this mask (which has already been cropped to the shapefile extent) to the cropped 
                # species file

                if fnf_exists:
                    data = np.ma.masked_array(cropped_image, mask=fnf_array)

            #  If there is no shapefile supplied, but there is a FNF mask supplied, apply the mask to each species raster
            elif fnf_exists: 
                data_unmasked = f.read()
                data = np.ma.masked_array(data_unmasked, mask=fnf_array)

            #  If no shapefile and no FNF are supplied, simply read the species file data:
            else:   
                data = f.read(1, masked=True)

            #  Find the number of masked pixels, for adjusting the calculations later so that they don't include masked pixels
            total_masked_pixels = ma.count_masked(data)
            
            #  Find the mean pixel value (which is a percentage) of all the unmasked pixels in the dataset 
            m = data.mean()
            
            #  Apply the threshold, selecting only species whose mean pixel value is below it
            if m < threshold:
                trace = True
                trace_species.append(file)
                trace_data.append(data)

            #  Otherwise, the species is not a trace species
            else:
                trace = False

            print(f'Species {name} : Mean: {m} - Trace : {trace}')
            
            #  Update the log / process_info list with a line for each species
            process_info.append(f'Species {name} : Mean: {m} - Trace : {trace}')

    return trace_species, trace_data, profile, process_info


def find_trace_by_dominance(domfile, fnf_exists, fnf_array, dom_threshold, shapes):
    """
    Takes in a dominant species raster, an optional FNF (forest/non-forest mask) and optional
    shapefile if you only want to analyse the pixels within certain areas.  The percentage of
    valid (ie not masked) pixels in which each species is dominant is calculated.  the threshold is then
    applied, and species contributing at a level below the threshold are set as "trace" species.
    
    Returns:
    
    trace_species : list - info on trace species including FIA codes
    process_info : list - textual info about the process
    """

    trace_species = []
    process_info = []

    #  Create connection to the dominant species raster GeoTIFF by opening with rasterio
    with rio.open(domfile) as dom:
        
        #  Calculate total pixels in the image for statistics later on
        total_pixels = dom.shape[0] * dom.shape[1]

        #  Display some metadata to help ensure things are working correctly 
        print(f'Shape of original dominant species raster: {dom.shape}\n')
        print(f'Total pixels : {total_pixels}\n')

        #  Append the same information to the "process_info" list which is effectively a log file
        process_info.append(f'Shape of original dominant species raster: {dom.shape}')
        process_info.append(f'Total pixels : {total_pixels}')


        if shapes and fnf_exists:

        #  If both a shapefile and a FNF mask are provided, first crop the dominant species raster using the features in the
        #  shapefile.  Then apply the FNF mask (which has already been cropped before the function was called).
        #  Also, update the values of "total_pixels" and "total_masked_pixels" in light of cropping the data.

            data, cropped_dom_transform = rasterio.mask.mask(dom, shapes, crop=True, nodata=-9999, filled=True)
            print(f"Cropped dominant species raster shape : {data.shape}\n")
            process_info.append(f"Cropped dominant species raster shape : {data.shape}")
            total_pixels = data.shape[1] * data.shape[2] 
            data = np.ma.masked_array(data, mask=fnf_array)
            total_masked_pixels = ma.count_masked(data)

        elif shapes:
        
        #  If only a shapefile is provided, crop the dominant species raster using the features in "shapes",
        #  then update the "total_pixels" and "total_masked_pixels"

            data, cropped_dom_transform = rasterio.mask.mask(dom, shapes, crop=True, nodata=-9999, filled=True)
            print(f"Cropped dominant species raster shape : {data.shape}\n")
            process_info.append(f"Cropped dominant species raster shape : {data.shape}")
            total_pixels = data.shape[1] * data.shape[2]
            total_masked_pixels = ma.count_masked(data)

        elif fnf_exists:

        #  If only a FNF mask is provided, use it to mask the dominant species raster,
        #  and update the number of masked pixels

            data = dom.read()
            data_fnf_masked = np.ma.masked_array(data, mask=fnf_array)
            total_masked_pixels = ma.count_masked(data_fnf_masked)
            data = data_fnf_masked
        
        else:

        #  Neither a shapefile nor a FNF mask is provided, so simply read the data and count the masked pixels

            data = dom.read(1, masked=True)
            total_masked_pixels = ma.count_masked(data)
        
        #  Calculate valid / data pixels as the difference between all pixels and all masked pixels

        data_pixels = total_pixels - total_masked_pixels

        #  Inform user of pixel info

        print(f'Number of masked pixels: {total_masked_pixels}')
        print(f'Number of data pixels: {data_pixels}')

        #  Update log file

        process_info.append('\nAfter mask and shapefile processing:\n')
        process_info.append(f'Number of masked pixels: {total_masked_pixels}')
        process_info.append(f'Number of data pixels: {data_pixels}')


    #  Use numpy.unique() function to count the number of times each species occurs in the dominant species raster.
    #  Each data pixel value is a species code, representing the species which is dominant in that pixel.  By counting
    #  all the pixels in which a particular species is dominant, a percentage can be found which tells the user what proportion of
    #  pixels a given species is dominant in.  This number can be compared to the threshold value to identify trace species.
    #  For example, if a species is dominant in fewer than 1% of valid pixels, it may be deemed "trace".
    
    unique, counts = np.unique(data, return_counts=True)  # calculate the number of times each unique pixel value occurs
    info = zip(unique, counts)  # zip unique raster value and their associated counts into an easily iterable object

    #  Iterate over the results, calculating percentage and applying the threshold for each unique species code
    #  u, c are just shorthand for unique value (ie species code), and count (ie number of times it occurs).

    for u, c in info:
        pc_dominant = (c / data_pixels) * 100
        if u != 0 or u != 0.0 or u != -9999:  # Filter out zeros or nodata (-9999) values
            
            #  Apply threshold.  If lower than the threshold, add the info string to the "trace_species" list

            if pc_dominant < dom_threshold:
                print(f'Species {u} occurs {c} times, and is dominant in {pc_dominant:.2f} % of pixels (TRACE)')
                trace_species.append(f'Species {u} occurs {c} times, and is dominant in {pc_dominant:.2f} % of pixels (TRACE)')

            #  If equal to or above threshold, this is not a trace species.  Inform the user with an info string.

            else:
                print(f'Species {u} occurs {c} times, and is dominant in {pc_dominant:.2f} % of pixels')

    return trace_species, process_info


def write_results_txt(trace_species, process_info, filename):
    
    """  Takes in some info in the form of lists of strings, and write them to text files  """

    with open(filename, 'w') as outfile:
        for info in process_info:
            outfile.write(f'{info}\n')
        outfile.write('\n')
        for trace in trace_species:
            outfile.write(f'{trace}\n')
    return

def write_trace_species_raster(data, profile, outname):

    """
    This function takes a list of numpy arrays (each representing a trace species)
    along with a profile dictionary and an output filename.
    It calculates the sum of the values of each layer of the stacked array - ie of the pixels "lying on top of each other"
    / or along axis=0.  It then writes this aggregate 2D array as a new raster.
    Effectively this represents the cumulative contributions of all species masked as "trace" in each pixel.
    """

    #  Stack the input data 
    all_layers = ma.stack(data)
    
    #  Add up all the pixels lying "on top of" each other in the stack
    layers_sum = all_layers.sum(axis=0)

    print(f'\n Writing output raster with summed trace species contributions: {outname}\n')

    #  Write the data
    with rio.open(outname, 'w', **profile) as dest:
        dest.write(layers_sum.filled(-9999))

    return

@click.command()
@click.option('--pc_adj', 'pc_adj_filepath', type=click.Path(exists=True), help="Path to directory containing percentage adjusted species rasters.  Use EITHER this method OR the --dom method below.")
@click.option('--dom','dominant_species_raster', type=click.Path(exists=True), help="Specifies that you want to calculate using a dominant species raster.  Provide such a raster file here")
@click.option('--out', 'outpath', default='.', type=click.Path(exists=True), help="Specify a directory path for the output files.  The default is the current directory")
@click.option('--aoi', default='LV_KC', help='The string code for the area of interest.  Default is LV_KC')
@click.option('--threshold', default=1, help="Threshold for calculating trace species.  Default is 1%")
@click.option('--fnf', default=None, help='Add a forest/non-forest mask - specify the filepath')
@click.option('--shapefile', default=None, help='Calculate trace species within this area - provide path to shapefile')
@click.option('--write_trace_raster', help='Write output trace species raster - provide output filename')
def main(pc_adj_filepath, dominant_species_raster, outpath, aoi, threshold, fnf, shapefile, write_trace_raster):
    
    fnf_exists = False
    fnf_array = None
    shapes = None

    if pc_adj_filepath and dominant_species_raster:
        print('Please choose either calculation by dominance or by total prevelance')
        sys.exit()

    if fnf and shapefile:
        #  If both FNF mask and shapefile are provided, start by extracting the shapes from the shapefile
        #  and cropping the FNF mask to keep only those areas within the shapefile polygons

        with fiona.open(shapefile, "r") as shp:
            shapes = [feature["geometry"] for feature in shp]
            with rio.open(fnf) as src:
                fnf_array, cropped_fnf_transform = rasterio.mask.mask(src, shapes, crop=True, nodata=-9999, filled=True)
                fnf_exists = True
        
    elif fnf:

        #  If FNF is provided (and no shapefile), read the FNF mask file and store the data in fnf_array
        
        process_info.append("\n FNF mask supplied...")
        
        with rio.open(fnf) as src:
            fnf_array = src.read()
            fnf_exists = True

    elif shapefile:

        #  If only a shapefile is supplied, extract the shapefile features, and store in a list called "shapes"
        
        process_info.append("\n Shapefile supplied...")

        with fiona.open(shapefile, "r") as shp:
            shapes = [feature["geometry"] for feature in shp]

    else:
        process_info.append("\n No FNF mask or shapefile supplied...\n")


    if pc_adj_filepath:
        #  If the adjusted percentage method is used (ie a path a pc_adj.{AOI}.tif is provided), then collect up the 
        #  necessary filenames into a list using glob, run the function "find_trace_by_pc_adj"
        #  then write the outputs
        pc_files = glob(os.path.join(pc_adj_filepath, f'*_pc_adj.{aoi}.tif'))
        
        # Print list of files found at the specified path:
        print(pc_files)
        
        #  Print threshold
        print(f'\nUsing threshold {threshold}%\n')
        
        #  Call the function
        trace_species, trace_data, profile, process_info = find_trace_by_pc_adj(pc_files, fnf_exists, fnf_array, threshold, shapes)
        
        #  Set up filenames for log outputs
        if fnf:
            filename = f'trace_results_by_pc_adj_fnf_{aoi}.txt'
        else:
            filename = f'trace_results_by_pc_adj_{aoi}.txt'   

        #  Write the resulting text files    
        write_results_txt(trace_species, process_info, filename)

        #  Write a trace species raster, if the flag "--write_trace_raster" is set
        if write_trace_raster:
            write_trace_species_raster(trace_data, profile, write_trace_raster)



    #  The other method of calculating trace species, based on a dominant species raster.
    if dominant_species_raster:

        dom_file = os.path.join(dominant_species_raster)
        
        #  Print out threshold and filename to confirm to user correct settings being used
        print(f'Dominant species raster : {dom_file}')
        print(f'Using threshold {threshold}%')

        #  Call the function
        trace_species, process_info = find_trace_by_dominance(dom_file, fnf_exists, fnf_array, threshold, shapes)
        
        #  Set up filenames for log outputs
        if fnf:
            filename = f'trace_results_by_dominance_fnf_{aoi}.txt'
        else:
            filename = f'trace_results_by_dominance_{aoi}.txt' 

        #  Write the resulting text files
        write_results_txt(trace_species, process_info, filename)


if __name__ == "__main__":
    main()