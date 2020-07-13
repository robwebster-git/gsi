#!/usr/bin/env python
# coding: utf-8

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
from pprint import pprint as p


def find_trace_by_pc_adj(files, fnf_exists, fnf_array, threshold, shapes):
    """
    Takes in a list of species rasters (percentage adjusted), an optional FNF (forest/non-forest mask) and and optional
    shapefile if you only want to analyse the area within certain areas.  The threshold is the percentage contirubution of a particular
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
        with rio.open(file) as f:
            #  If there's a shapefile supplied, crop each image to the shapefile and update metadata
            if shapes:
                cropped_image, cropped_transform = rasterio.mask.mask(f, shapes, crop=True)
                cropped_meta = f.meta
                #rasterio.plot.show(cropped_image, transform=cropped_transform)
                
                # If a FNF mask is supplied, apply this mask (which has already been cropped to the shapefile extent) to the cropped 
                # species file
                if fnf_exists:
                    data = np.ma.masked_array(cropped_image, mask=fnf_array, fill=-9999)
                    np.ma.set_fill_value(data, -9999)
                    p(data)
                    p(cropped_meta)
                    #rasterio.plot.show(data, transform=cropped_transform)
        
            ####   SHAPEFILE CROPPING NOT FULLY DEBUGGED YET!!!
            
            # get a string of the species name, for readability of outputs and for showing progress 
            name = file.split('_pc_adj')[0][2:]

            print(f'Processing {name}')

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

            #print(f'Species {name} : Mean: {m} - Trace : {trace}')
            
            #  Update the log / process_info list with a line for each species
            process_info.append(f'Species {name} : Mean: {m} - Trace : {trace}')

    return trace_species, trace_data, profile, process_info


def find_trace_by_dominance(domfile, fnf_exists, fnf_array, dom_threshold):

    trace_species = []
    process_info = []

    with rio.open(domfile) as dom:
        print(f'Shape of dominant species raster: {dom.shape}')
        process_info.append(f'Shape of dominant species raster: {dom.shape}')

        total_pixels = dom.shape[0] * dom.shape[1]
        print(f'Total pixels : {total_pixels}')
        process_info.append(f'Total pixels : {total_pixels}')

        if fnf_exists:
            data = dom.read()
            data_fnf_masked = np.ma.masked_array(data, mask=fnf_array)
            total_masked_pixels = ma.count_masked(data_fnf_masked)
        else:    
            data = dom.read(1, masked=True)
            total_masked_pixels = ma.count_masked(data)
        
        print(f'Number of masked pixels: {total_masked_pixels}')
        process_info.append(f'Number of masked pixels: {total_masked_pixels}')
        
        data_pixels = total_pixels - total_masked_pixels
        
        print(f'Number of data pixels: {data_pixels}')
        process_info.append(f'Number of data pixels: {data_pixels}')
        
    unique, counts = np.unique(data, return_counts=True)
    info = zip(unique, counts)

    for u, c in info:
        pc_dominant = (c / data_pixels) * 100
        if u != 0 or u != 0.0 or u != -9999:
            if pc_dominant < dom_threshold:
                trace_species.append(f'Species {u} occurs {c} times, and is dominant in {pc_dominant:.2f} % of pixels (TRACE)')
            else:
                print(f'Species {u} occurs {c} times, and is dominant in {pc_dominant:.2f} % of pixels')

    return trace_species, process_info


def write_results_txt(trace_species, process_info, filename):
    with open(filename, 'w') as outfile:
        for info in process_info:
            outfile.write(f'{info}\n')
        outfile.write('\n')
        for trace in trace_species:
            outfile.write(f'{trace}\n')


def write_trace_species_raster(data, profile, outname):

    all_layers = np.stack(data)
    
    layers_sum = all_layers.sum(axis=0)
    print(layers_sum.shape)
    print(data[0].shape)

    with rio.open(outname, 'w', **profile) as dest:
        dest.write(layers_sum.filled(-9999), 1)

    return

@click.command()
@click.option('--pc_adj', 'pc_adj_filepath', type=click.Path(exists=True))
@click.option('--dom','dominant_species_raster', type=click.Path(exists=True))
@click.option('--out', 'outpath', default='.', type=click.Path(exists=True))
@click.option('--pc_threshold', default=1)
@click.option('--dom_threshold', default=1)
@click.option('--fnf', default=None, help='Add a forest/non-forest mask')
@click.option('--shapefile', default=None, help='Calculate trace species within this area')
@click.option('--write_trace_raster', is_flag=True, default=False, help='Write out trace species raster')
def main(pc_adj_filepath, dominant_species_raster, outpath, pc_threshold, dom_threshold, fnf, shapefile, write_trace_raster):
    
    fnf_exists = False
    fnf_array = None
    shapes = None

    if pc_adj_filepath and dominant_species_raster:
        print('Please choose either calculation by dominance or by total prevelance')
        sys.exit()

    if fnf and shapefile:
        with fiona.open(shapefile, "r") as shp:
            shapes = [feature["geometry"] for feature in shp]
            with rio.open(fnf) as src:
                fnf_array, cropped_fnf_transform = rasterio.mask.mask(src, shapes, crop=True)
                fnf_exists = True
        
    elif fnf:
        with rio.open(fnf) as src:
            fnf_array = src.read()
            fnf_exists = True

    elif shapefile:
        with fiona.open(shapefile, "r") as shp:
            shapes = [feature["geometry"] for feature in shp]

    if pc_adj_filepath:
        pc_files = glob(os.path.join(pc_adj_filepath, '*_pc_adj.LV_KC.tif'))
        print(pc_files)
        print(f'Using threshold {pc_threshold}%')
        trace_species, trace_data, profile, process_info = find_trace_by_pc_adj(pc_files, fnf_exists, fnf_array, pc_threshold, shapes)
        if fnf:
            filename = 'trace_results_by_pc_adj_fnf.txt'
        else:
            filename = 'trace_results_by_pc_adj.txt'   
        write_results_txt(trace_species, process_info, filename)

        if write_trace_raster:
            write_trace_species_raster(trace_data, profile, 'trace_species2.tif')

    if dominant_species_raster:
        dom_file = os.path.join(dominant_species_raster)
        print(f'Dominant species raster : {dom_file}')
        print(f'Using threshold {dom_threshold}%')
        trace_species, process_info = find_trace_by_dominance(dom_file, fnf_exists, fnf_array, dom_threshold, shapes)
        if fnf:
            filename = 'trace_results_by_dominance_fnf.txt'
        else:
            filename = 'trace_results_by_dominance.txt'   
        write_results_txt(trace_species, process_info, filename)


    #tifs = list(map(rio.open, files))
    #data = np.stack(list(map(lambda x: x.read(1), tifs)))
    #profiles = list(map(lambda x: x.profile, tifs))

if __name__ == "__main__":
    main()