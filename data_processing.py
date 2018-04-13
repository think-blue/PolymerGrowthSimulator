import numpy as np
import pandas as pd
import simulation as sim
import matplotlib.pyplot as plt

exp_df = pd.read_excel('polymer_20k.xlsx')
exp_molmass = exp_df[exp_df.columns[0]].values
exp_values = exp_df[exp_df.columns[1]].values
exp_chainlen = (exp_molmass/99.13).astype(int)
exp_cl_min = exp_chainlen.min()
exp_chainlen = np.concatenate((np.arange(1,exp_cl_min),exp_chainlen))
exp_cl_val = dict(zip(exp_chainlen,exp_values))
exp_cl_val.update(zip(np.arange(1,exp_cl_min+1),np.zeros(exp_cl_min)))
exp_val = np.array(list(exp_cl_val.values()))
exp_val = np.concatenate((np.zeros(exp_cl_min-1),exp_val))


sim_data = np.load('output.npy')
sim_data = sim.polymer(1000, 100000, 0.3, 0.02, 0.3, 0, 0, 1, -10000000, 0.23, 0.23, 0.5)
dead, living, coupled = sim_data
sim_data = np.concatenate((dead, living, coupled))
sim_cl_max = sim_data.max()
sim_val,sim_bins = np.histogram(sim_data,bins=np.arange(sim_cl_max+1))

diff = int(sim_cl_max - exp_val.shape[0])

if diff > 0:	 
	exp_val = np.concatenate((exp_val, np.zeros(abs(diff))))
elif sim_cl_max < exp_val.shape:
	sim_val = np.concatenate((sim_val, np.zeros(abs(diff))))

exp_val_max = exp_val.max()
exp_norm = (exp_val)/exp_val_max

sim_val_max = sim_val.max()
sim_norm = (sim_val)/sim_val_max

fig, axes = plt.subplots(ncols=2)

axes[0].bar(np.arange(exp_norm.shape[0]),exp_norm)
axes[1].bar(np.arange(sim_norm.shape[0]),sim_norm)
plt.show()
print(np.sum(abs(exp_norm - sim_norm)))

# dead, living, coupled = [(180+99.13*x) for x in sim_data]
