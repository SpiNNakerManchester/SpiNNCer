from colorama import Fore, Style, init as color_init
from spinncer.utilities.constants import (CELL_NAME_FOR_ID,CELL_IO_STATUS,IO_Status)
import numpy as np

def population_reporting(positions):
    """
    Helper function reporting on various aspects of the Populations
    :param positions: X Y Z positions of individual cells
    :type positions: np.ndarray
    :return: None
    :rtype: None
    """
    color_init(strip=False)
    unique_ids = np.unique(positions[:, 1]).astype(int)
    print("=" * 60)
    print("The file contains information about",
          Fore.GREEN, unique_ids.size, Style.RESET_ALL, "populations")
    print("=" * 60)
    print("Mapping of IDs:")
    print("-" * 60)
    for ui in unique_ids:
        _cell_name = CELL_NAME_FOR_ID[ui]
        _status = CELL_IO_STATUS[_cell_name]
        if _status == IO_Status.INPUT:
            _color = Fore.GREEN
        elif _status == IO_Status.OUTPUT:
            _color = Fore.RED
        else:
            _color = ''
        print("\t{:2d} -> {:10} ".format(ui, _cell_name),
              _color, "[{:16}]".format(_status), Style.RESET_ALL)
    print("=" * 60)
    print("Number of neurons in each population")
    print("-" * 60)
    for ui in unique_ids:
        _cell_name = CELL_NAME_FOR_ID[ui]
        _no_cells = positions[positions[:, 1] == ui, :].shape[0]
        print("\t{:10} -> {:10} ".format(_cell_name, _no_cells))