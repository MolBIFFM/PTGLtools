import os
import urllib.request

def download_pdb(pdbcode, datadir="./", downloadurl="https://files.rcsb.org/download/"):
    """
    Taken from stackoverflow (Andras Aszodi) and modified.
    Downloads a PDB file from the Internet and saves it in a data directory.
    :param pdbcode: The standard PDB ID e.g. '3ICB' or '3icb'
    :param datadir: The directory where the downloaded file will be saved
    :param downloadurl: The base PDB download URL, cf.
        `https://www.rcsb.org/pages/download/http#structures` for details
    :return: the full path to the downloaded PDB file or None if something went wrong
    """
    pdbfn = pdbcode + ".cif"
    url = downloadurl + pdbfn
    outfnm = os.path.join(datadir, pdbfn)
    try:
        urllib.request.urlretrieve(url, outfnm)
        return outfnm
    except Exception as err:
        print(str(err), file=sys.stderr)
        return None

import numpy as np
import pandas as pd

def setViewToTransformationMatrix(setViewString: str):
    """ Use the set_view string of pymol obtained by the get_view command and prepare the transformation matrix """
    array: list = []
    row: list = []
    columns: int = 0
    rows: int = 0
    for item in setViewString.split():
        if columns >= 3:
            rows+=1
            columns = 0
            array.append(row)
            row: list = []
        if rows < 3:
            try:
                item = float(item.rstrip(","))
                columns+=1
                row.append(item)
            except:
                pass                   
        else:
            return np.array(array)

def createLayout(setViewString: str, coord_orig, alpha=False):
    """ Creates a new layout for the vertices based on the geometrical centers of the chains (taken from the pseudoatoms of the PyMOL representation. If needed the z coordinates can be written out and used as alpha values for depth view. """
    # The transformation matrix is prepared to modify the coordinates of the pseudoatoms to match the desired view 
    transform_matrix = setViewToTransformationMatrix(setViewString)
    # The coordinates are transformed and the new layout is created
    coord_trans = coord_orig.dot(transform_matrix)
    new_layout: list = []
    if alpha:
        alphaValues: list = []
    maximum_y_value: float = coord_trans[0][1]
    for line in coord_trans:
        if line[1] > maximum_y_value:
            maximum_y_value = line[1]
    if alpha:
        for line in coord_trans:
            new_layout.append([round(line[0],3),round(-line[1]+maximum_y_value,3)])
            alphaValues.append(round(line[2],3))
        xmin = min(alphaValues) 
        xmax=max(alphaValues)
        for i, x in enumerate(alphaValues):
            #alphaValues[i] = 0.25 + ((x-xmin) / (1.5*(xmax-xmin))) # Normalization to the range of 0.25-1
            alphaValues[i] = 0.5 + ((x-xmin) / (2*(xmax-xmin))) # Normalization to the range of 0.5-1
        return new_layout, alphaValues
    else:
        for line in coord_trans:
            new_layout.append([round(line[0],3),round(-line[1]+maximum_y_value,3)])
        return new_layout

def modifyLayoutPosition(vertexName: str, layout, graph, coordinates:tuple = [0.0, 0.0]):
    """ Modify individual positions of the vertices in the representation. """
    for index, vertex in enumerate(graph.vs['label']):
        if vertex == vertexName:
            layout[index] = coordinates
    return layout

