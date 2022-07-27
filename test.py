import ReliabilityDiagram as rd
import numpy as np


####################################################################
# Load data
folder_data = "C:/Users/cmllecoz/PycharmProjects/reliability_diagram/data"

fc_name = "ecmwf"
for k in range(6):
    print("Week {}".format(k+1))
    fcast = np.load("{}/forecasts_ecmwf_week{}.npy".format(folder_data,k+1))
    clima = np.load("{}/climatology_week{}.npy".format(folder_data,k+1))
    obs = np.load("{}/observations_week{}.npy".format(folder_data,k+1))

    data = rd.ReliabilityDiagram(obs,fcast,clima,0,1/3,closed_ends='both',nbins=10)

    tab_ls = data.contingency_table()
    print(tab_ls)
    print(np.sum(tab_ls,axis=1))
    exit()
    #print(tab_us)
    #rd.plot_diagram(tab_us, tab_ls, tab_ms, prob_bins=np.arange(0,1,0.1),save_file="{}/{}_week{}.png".format(folder_results,fc_name,k+1))