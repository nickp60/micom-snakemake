from micom.workflows import grow, save_results, build, tradeoff, load_results
import os
import sys
import glob
import pandas as pd

from zipfile import BadZipFile

def merge_gr(gra, grb):
    from copy import deepcopy
    new = deepcopy(gra)
    _exchanges = pd.concat([gra.exchanges, grb.exchanges], ignore_index=True)
    _growth_rates =  pd.concat([gra.growth_rates, grb.growth_rates], ignore_index=True)
    _annotations =  pd.concat([gra.annotations, grb.annotations], ignore_index=True).drop_duplicates().reset_index(drop=True)
    new.exchanges = _exchanges
    new.growth_rates = _growth_rates
    new.annotations = _annotations
    return(new)

if __name__ == "__main__":
    infiles = glob.glob(sys.argv[1] + "*.zip")
    print(infiles)
    growth = load_results(infiles[0])
    for i in range(2, len(infiles)):
        print(i)
        try:
            growthb = load_results(infiles[i])
            growth = merge_gr(growth, growthb)
            print(growth.growth_rates.shape)
        except BadZipFile:
            print(f"coulnt load {infiles[i]}")
    save_results(growth, sys.argv[2])
