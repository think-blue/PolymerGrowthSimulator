import numpy as np
import pandas as pd
import simulation as sim
import matplotlib.pyplot as plt


class Difference:
    def __init__(self, file_name):
        self.exp_df = pd.read_excel(file_name)
        self.exp_molmass = self.exp_df[self.exp_df.columns[0]].values
        self.exp_values = self.exp_df[self.exp_df.columns[1]].values
        self.exp_chainlen = (self.exp_molmass / 99.13).astype(int)
        self.exp_cl_min = self.exp_chainlen.min()
        self.exp_chainlen = np.concatenate((np.arange(1, self.exp_cl_min), self.exp_chainlen))
        self.exp_cl_val = dict(zip(self.exp_chainlen, self.exp_values))
        self.exp_cl_val.update(zip(np.arange(1, self.exp_cl_min + 1), np.zeros(self.exp_cl_min)))
        self.exp_val = np.array(list(self.exp_cl_val.values()))
        self.exp_val = np.concatenate((np.zeros(self.exp_cl_min - 1), self.exp_val))

    def get_difference(self, sim_data):
        # sim_data = sim.polymer(1000, 100000, 0.3, 0.02, 0.3, 0, 0, 1, -10000000, 0.23, 0.23, 0.5)
        dead, living, coupled = sim_data
        sim_data = np.concatenate((dead, living, coupled))
        sim_cl_max = sim_data.max()
        sim_val, sim_bins = np.histogram(sim_data, bins=np.arange(sim_cl_max + 1))

        diff = int(sim_cl_max - self.exp_val.shape[0])

        if diff > 0:
            exp_val = np.concatenate((self.exp_val, np.zeros(abs(diff))))
        elif sim_cl_max < self.exp_val.shape:
            sim_val = np.concatenate((sim_val, np.zeros(abs(diff))))

        exp_val_max = self.exp_val.max()
        exp_norm = (self.exp_val) / exp_val_max

        sim_val_max = sim_val.max()
        sim_norm = (sim_val) / sim_val_max

        # fig, axes = plt.subplots(ncols=2)
        #
        # axes[0].bar(np.arange(exp_norm.shape[0]), exp_norm)
        # axes[1].bar(np.arange(sim_norm.shape[0]), sim_norm)
        # plt.show()
        print(np.sum(abs(exp_norm - sim_norm)))
        return np.sum(abs(exp_norm - sim_norm))

    # dead, living, coupled = [(180+99.13*x) for x in sim_data]
