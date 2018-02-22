import numpy as np
import matplotlib as plt

def polymer(number_of_molecules, time_sim, p_growth, p_death, p_dead_react, kill_spawns_new, video, coloured, monomer_pool, l_exponent, d_exponent, l_naked):
#this function simulates the growth of polymers it takes;
# number_of_molecules - the number of starting chains (length 1)
# time_sim - the number of timesteps the simulation runs for
# p_growth - the probability of growth for each monomer, each time step
# p_death - the probability of a monomer dying and joining the dead pool
# p_dead_react - the probilty a dead polymer reacts with a living one
# killer_pool_size - the size of the pool of killer molecules in the system
# p_killer - the probability a killer molecule reacts with a polymer
# killer_pool_reusable - 0 means the killer pool is depleted every time a
# killer molecule is used, 1 means the killer pool is never used up
# kill_spawns_new - a binary flag as to whether a kill event means a new
# polymer of length 1 is spawned or not 1=killing event spawns a new
# polymer chain, 0 means it doesn't.
# video - binary flag on whether to output to video
# coloured - binary flag on whether to output the living, dead and coupled
# polymers as separate colours in a stacked histogram. 
# monomer_pool - the size of the monomer pool, if it is a negative number
# then no monomer pool is used and monomer supply is assumed infinite,
# p_growth is now p_growth as provided * monomer_pool size
# l_exponent - the exponent of the living length in the probability
# calculation for successful coupling
# d_exponent - the exponent of the dead length in the probability
# calculation for successful coupling
#
#
# It returns an array with the lengths of all polymers in the system
# (living and dead)


#Now let's write a subfunction to make the histogram
	def make_histogram(living, dead, coupled, coloured, initial_monomer, current_monomer, time):
		d = [living, dead, coupled]
#calculate M_n, M_w and PDI
		DPn = np.mean(d)
		DPw = np.sum(np.square(d))/(DPn*d.shape[0])
		PDI = Dpw/DPn
		conversion = 1-current_monomer/initial_monomer
		#**dlmwrite('polymerOutput.txt',[time, conversion, DPn, DPw, PDI], '-append');

		if coloured == 0:
			plt.hist(d, np.max(d)-min(d), bins='auto')
		else:
			step = np.ceil((np.max(d)-np.min(d))/1000)
			binEdges = np.arange(np.min(d)-0.5, step, np.max(d)+0.5)
			living_counts = np.histogram(living, binEdges)[0]
			dead_counts = np.histogram(dead, binEdges)[0]
			coupled_counts = np.histogram(coupled, binEdges)[0]
			midbins = binEdges[0:-1-1] + (binEdges[2:-1]-binEdges[1:-1-1])/2
			if coupled = []:
				pass
				#bar(midbins,[dead_counts;living_counts]','stacked')
            	#legend({'Dead','Living'});
            else:
            	pass
            	#bar(midbins,[coupled_counts;dead_counts;living_counts]','stacked')
            	#legend({'Terminated','Dead','Living'});
        #xlabel('Length in units')
	    #ylabel('Frequency')
    	#title(strcat(' conversion=',num2str(conversion),' time=',num2str(time), ' PDI=',num2str(PDI),' DPn=',num2str(DPn),' DPw=',num2str(DPw)))   	

#setup file for output

	#file=fopen('polymerOutput.txt','w');
	#fprintf(file,'%s\n','Time, Conversion, DPn, DPw, PDI');
	#fclose(file);
#the pool of living polymers is an array of 1's of the right size	
	living = np.ones(1, number_of_molecules)
#the pool of dead polymers 
	dead = []
	coupled = []
#if we're outputing video, set it all up
	if video == 1:
		#v=VideoWriter('polymer_distribution_video.avi');
    	#open(v);	

    initial_monomer_pool = monomer_pool
    fig = plt.figure()
    for t in range(time_sim):
    	if living.shape[0]>0:
#first make a random vector with uniform random numbers which will
#decide the fate of each living polymer	
			r = np.random.rand(living.shape[0])
#Now if the random number for a polymer is below p_growth, it will grow.
##### LIVING GROWTH #####
			if monomer_pool<0:
				living[(r<p_growth)-1] = living[(r<p_growth)-1]+1
				monomer_ratio = 1
			else:
				monomer_ratio = monomer_pool/initial_monomer_pool
				living[(r<(p_growth*monomer_ratio))-1] = living[(r<(p_growth*monomer_ratio))-1]+1
				used_monomer = np.sum(r<(p_growth*monomer_ratio))
				if used_monomer > monomer_pool:
					monomer_pool = 0
				else:
					monomer_pool = monomer_pool - used_monomer
#Next if the random number for a polymer is above p_growth but below
#p_growth+(p_death*monomer_ratio) it will die store this in new dead for now, until
#we've worked out if any of the old dead react, We multiply by
#monomer ratio, as we believe monomer is also involved here.
#However, if we have infinite monomer, monomer ratio is set to 1
#above, so then this will be just p_death)
###### KILLING ######
			if kill_spawns_new == 1:
#since the kill starts a new chain it uses up 1 monomer, so we
#decrease the monomer pool
				if monomer_pool>0:
					new_dead = living[(r<(p_growth*monomer_ratio+(p_death*monomer_ratio)) & r>=p_growth*monomer_ratio)-1]
					living[(r<(p_growth*monomer_ratio+(p_death*monomer_ratio)) & r>=p_growth*monomer_ratio)] = []
	
#So in the new system each dead chain chooses another chain from
#the system (either living or dead) to attack
        	which_chain_attacked_per_dead = np.ceil(np.random.rand(dead[0])*(dead.shape[0]+living.shape[0]));
#if the chosen chain is a dead one, we do nothing, so now only
#consider when the chosen number is above the nunber of dead chains
#let's implement this in a loop to consider each dead chain
#individually
        	
        	still_dead = np.ones(dead.shape[0])

        	for dead_counter in range(dead.shape[0]):
#if it chooses to attack a living chain...
				if which_chain_attacked_per_dead[dead_counter]>dead.shape[0]:
					which_living_attacked = which_chain_attacked_per_dead[dead_counter]-dead.shape[0]
					r_success = np.random.rand()
#calculate the probability of sucess given the formula from Bryn
					p_success = p_dead_react/(living[which_living_attacked]**np.min(living[which_living_attacked]*(l_exponent/l_naked), l_exponent)*dead[dead_counter]**np.min(dead[dead_counter]*(d_exponent/l_naked), d_exponent)) #the dead pool
					
					if r_success<p_success:
						living[which_living_attacked] = living[which_living_attacked]+dead[dead_counter]
						still_dead[dead_counter] = 0

			dead = dead[still_dead==1]
			dead = [dead new_dead]

		if video==1:
			make_histogram(living, dead, coupled, coloured, initial_monomer_pool, monomer_pool,t)	
			#frame=getframe(gcf);
        	#writeVideo(v,frame);
	distribution = [living dead coupled]
	make_histogram(living, dead, coupled, coloured, initial_monomer_pool, monomer_pool,t) 
	#if video==1:
		#close(v)
	return distribution