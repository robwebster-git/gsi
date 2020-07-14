import numpy as np
import fiona
import rasterio as rio
from rasterio import features
from shapely.geometry import shape
import time
import click

@click.command()
@click.option('--domfile', default='trace_test_files/dominant_maine.tif', type=click.Path(exists=True), help="Path to dominant species raster")
@click.option('--threshold', default=1000, help='Minimum area in m2.  1000m2 = 10 contigous pixels')
def main(domfile, threshold):
    
    with rio.open(domfile) as src:
        dom = src.read(1, masked=True)
        profile = src.profile

    shapes = features.shapes(dom, transform=src.transform)

    shape_list = list(shapes)

    valid_shapes = []

    # temporarily hard-coded but can be got from the trace_species.py script.  This example is from Top of Maine area.
    trace_species = [71, 94, 95, 125, 129, 261, 315, 316, 319, 531, 543, 741, 761, 951]

    for shp in shape_list:
        a = shape(shp[0])
        v = int(shp[1])
        if a.area > threshold and v in trace_species:
            valid_shapes.append((v,a.area))


    valid_shapes_array = np.array(valid_shapes, dtype=np.int64)
    values, counts = np.unique(valid_shapes_array[:,0], return_counts=True)
    info = zip(values, counts)

    for v, c in info:
        print(f'Species {v} occurs {c} times')



if __name__ == "__main__":
    start = time.time()
    main()
    end = time.time()
    print(f'Execution took {end - start} seconds')