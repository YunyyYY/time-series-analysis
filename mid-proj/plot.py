import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


avo = pd.read_csv("avocado.csv")
_df = pd.read_csv("merge.csv")


def plot_trend(region):
    plt.clf()
    avo[(avo.region == region) & (avo.type == "conventional")].AvgPrice.plot.line(label="conv")
    avo[(avo.region == region) & (avo.type == "organic")].AvgPrice.plot.line(label="org")
    plt.legend()
    plt.savefig(f'img/{region}.png')
    
    
def plot_diff(region):
    plt.clf()
    _df[_df.region == region]['diff'].plot.line()
    plt.savefig(f'img/{region}_diff.png')
    

def main():
    for region in avo.region.unique():
        plot_diff(region)
        # plot_trend(region)
    
    
if __name__ == '__main__':
    main()