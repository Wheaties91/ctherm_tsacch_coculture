#!/usr/bin/python
#! python 3.9
#try to specify that we will use python version 3.9
__author__ = "Wheaton Schroeder"
#latest version: 02/17/2022
#written to run to coculture model. Also written to debug to coculture class as it is in development

#imports, again, these are used previously, not sure of their necessity
import cobra
from fva import FVA
from fba import FBA
from steadycom import SteadyCom
from datetime import datetime
from cobra import Model, Reaction, Metabolite
import re
import os
import warnings

#get the current time and date for tracking solution time
start_time = datetime.now()
print("Starting time: ",start_time)

#get the current directory to use for importing things
curr_dir = os.getcwd()

"""
#import the ctherm model
#files for cross-talk community
model1 = cobra.io.read_sbml_model(curr_dir + "/iCTH669_comm.sbml")
model2 = cobra.io.read_sbml_model(curr_dir + "/iTSA525_comm.sbml")

"""
#files for non-cross-talk community
model1 = cobra.io.read_sbml_model(curr_dir + "/iCTH669_w_GLGC_non_comm.sbml")
model2 = cobra.io.read_sbml_model(curr_dir + "/iTSA525_non_comm.sbml")


#start the steadycom object
#note: "R_" characters get removed when the model is read apparently so omit here
print("building community object...")
comm_obj = SteadyCom(model1,model2,"EXCH_",log_file='steadycom_log_test.txt')

cth_abund = 0.58125

#define relative abundances
abundance = {
    "iCTH669":cth_abund,
    "iTSA525":1-cth_abund,
}

print("building community abundances...")
#apply abundances to the community
comm_obj.define_abundance(abundance)

#note the biomass equations for each species
#note: "R_" characters get removed when the model is read apparently so omit here
biomass_eqns = {

    "iCTH669":"BIOMASS",
    "iTSA525":"biomass_target" 

}

#define bigM for the model
bigM = 1000

#define the media components
#the community will only be allowed to uptake species defined in the media, and at the rate of at most the number defined here
#positive means the maximum allowed community-level uptake
media = {

    "h_e":bigM,
    "nh4_e":bigM,
    "h2o_e":bigM,
    "ca2_e":bigM,
    "mg2_e":bigM,
    "k_e":bigM,
    "so4_e":bigM,
    "pi_e":bigM,
    "fe3_e":bigM,
    "na1_e":bigM,
    "cu2_e":bigM,
    "cellb_e":(5/2),
    "xylb_e":3

}

#apply the medium to the model
print("defining medium...")
comm_obj.define_medium(media)

#this builds the community model, makes constraints based on the relative abundances, this may take a while since all constraints need to be scaled
print("building community...")
comm_obj.build_comm_x(biomass_eqns)      #returns nothing

#solve for maximum growth rate for the given community composition
print("running steadycom...")
mu_soln = comm_obj.max_mu()

#report results
output=open('steadycom_results_test.txt','w',buffering=1)
output.write("SEADYCOM MAX MU SOLUTION REPORT\n")
output.write("----------------------------------------------------\n\n")

#check if an exception occured
if mu_soln['exception']:

    #if here, an exception occured, we have no solution
    output.write("time to solve: "+str(mu_soln['soln_time'])+"\n")
    output.write("exception occured, no solution exception: \n"+mu_soln['exception_str']+"\n")

else:

    #create a formatted string to report a large number of decimal places
    formatted_string1 = "{:.8f}".format(mu_soln['mu_objective'])
    formatted_string2 = "{:.8f}".format(mu_soln['flux_objective'])

    #if here no exception, we have a solution to write about
    #write our results to the output file
    output.write("time to solve: "+str(mu_soln['soln_time'])+"\n")
    output.write("status: "+str(mu_soln['status'])+"\n")
    output.write("objective value (mu): "+formatted_string1+"\n")
    output.write("objective value (flux sum): "+formatted_string2+"\n\n")
    output.write("Community-based reaction rates (V^k_j):")
    output.write("reaction\tlb\tlevel\tub\n")
    output.write("----------------------------------------------------------------------------------\n")

    print("objective value (mu): ",formatted_string1)
    print("objective value (flux sum): ",formatted_string2)
    print("model status: "+str(mu_soln['status']))

    #for each reaction in the first model
    for rxn in model1.reactions:

        #make the reaction id that is used in the community model
        comm_rxn_id = rxn.id+"_"+model1.id

        #write the FBA results
        output.write(model1.id+"\t"+rxn.id+"\t"+str(mu_soln[comm_rxn_id]['lb'])+"\t"+str(mu_soln[comm_rxn_id]['flux'])+"\t"+str(mu_soln[comm_rxn_id]['ub'])+"\n")

    #also write for each reaction in the secondd model
    for rxn in model2.reactions:

        #make the reaction id that is used in the community model
        comm_rxn_id = rxn.id+"_"+model2.id

        #write the FBA results
        output.write(model2.id+"\t"+rxn.id+"\t"+str(mu_soln[comm_rxn_id]['lb'])+"\t"+str(mu_soln[comm_rxn_id]['flux'])+"\t"+str(mu_soln[comm_rxn_id]['ub'])+"\n")

    #community exchanges
    output.write("\n\nCOMMUNITY EXCHANGES\n")
    output.write("metabolite\tx_c\t"+model1.id+"\t"+model2.id+"\n")
    output.write("------------------------------------------------------\n")
    
    #write the community exchange reactions
    for met in mu_soln['x_c']:

        #create a string for the contributions of the individual models
        exch_1 = "EXCH_"+met+"_"+model1.id
        exch_2 = "EXCH_"+met+"_"+model2.id

        #essentially since we don't know what fluxes are associated with reactions are associated with each exchange, use try/except structure to identify them in the simplest manner
        try:

            #write community and individual model exchanges

            #try writing both
            output.write(met+"\t"+str(mu_soln['x_c'][met])+"\t"+str(mu_soln[exch_1]['flux'])+"\t"+str(mu_soln[exch_2]['flux'])+"\n")

        except:

            #if there isn't both, figure out which one to write

            try:

                #write community and individual model exchanges
                output.write(met+"\t"+str(mu_soln['x_c'][met])+"\t0\t"+str(mu_soln[exch_2]['flux'])+"\n")

            except:

                try:

                    #write community and individual model exchanges
                    output.write(met+"\t"+str(mu_soln['x_c'][met])+"\t"+str(mu_soln[exch_1]['flux'])+"\t0\n")
                
                except:

                    #if here reaction doesn't exist, need to have a look
                    output.write(met+"\t"+str(mu_soln['x_c'][met])+"\t0\t0\tDNE\n")
        

#get the current time and date for tracking solution time
end_time = datetime.now()
print("Ending time: ",end_time)

print("success!\n")