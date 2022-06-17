#!/usr/bin/python

#try to specify that we will use python version 3.9
__author__ = "Wheaton Schroeder"
#latest version: 01/21/2022

#written to implement Flux Balance Analysis (FBA), Flux Variability Analysis (FVA), parsimonious FBA (pFBA) and other stoichiometric
#modeling tools on a community model of two species. Note that this is derived from the SteadyCom framework. 

#These are the imports I have used in past for cobrapy, so will use them here as well, not sure the the necessity of any of these
#from __future__ import absolute_import

#notes on some formatting requiresments of the models
#1) each exchange reaction should act only on a single metabolite
#2) no two models should have the same ID


from optlang.interface import OPTIMAL, FEASIBLE, INFEASIBLE, ITERATION_LIMIT, NUMERIC, SUBOPTIMAL, TIME_LIMIT
from optlang.symbolics import Zero, add

import os
import sys
import warnings
import re
import cobra
from datetime import datetime
from fba import FBA
from fva import FVA

import copy

from cobra import Model, Reaction, Metabolite, Solution

#now that we have defined the import library, let us create a class for the mintransfers algorithm
class SteadyCom(object):

    #initialization of class:
    #self - needs to be passed itself
    #model1 - first model of the coculture
    #model2 - second model of the coculture
    #exch_tag - unique string tag that identifies an exchange reaction
    #note initiation is always for two models, will add a function to add additional models later
    #note: this only really works if each exchange reaction has only one exchanged metatolite!
    def __init__(self,model1,model2,exch_tag,log_file='steadycom_log.txt',bigM=10000):

        #define an output log file which may be useful for debugging purposes
        #log file for building the community
        self.log=open(log_file,'w',buffering=1)
        
        #define a combined model, start it out as as copy of model 1
        self.combined_model = model1.copy()

        #need to sprinkle these around whenever changing the model so changes stick correctly
        self.combined_model.solver.update()
        self.combined_model.repair()

        #initialize a dictionary for exchange reaction sets
        self.exch_sets = {}

        #create a list of the models/members of the community
        #note this list will contain the models, not just the model names
        self.members = []

        #let the user define a bigM, otherwise just use a large default value
        self.bigM = bigM

        #add the members
        self.members.append(model1)
        self.members.append(model2)

        #initialize a dictionary of X^K values, where value is the relative abundance of member species, key is the model ID
        #this list will be populated later
        self.X_k = {}

        #add an attribute to the combined model noting which model the reactions originally belonged to
        for rxn in self.combined_model.reactions:

            #add origin attribute, use the original model ID for the origin attribute
            setattr(rxn,'origin',model1.id)

            #add origin tag to the reaction id, metaid, and name sp don't get "ignoring reaction since it already exists" issue
            rxn.id = rxn.id + "_" + model1.id
            rxn.name = rxn.name + " " + model1.name

            #while we are at it find exchange reactions, set attributes of if exchange reaction and exchanged metabolite
            #search using the regulat expression tag
            if bool(re.search(exch_tag,rxn.id)):

                #state that it is an exchange
                setattr(rxn,"isexch",True)

                #set metabolite which is exchanged
                rxn_mets = rxn.metabolites

                #set the metaboltes it is exchanging
                #note that this will be a dictionary so need to get bare metabolite ID

                #get the metabolite key
                met_key = str(rxn_mets.keys())

                #split the key to get just the metabolite ID
                split_key = re.findall(r'\w+', met_key)

                setattr(rxn,"exchof",split_key[2])
                setattr(rxn,"exchstoich",rxn.get_coefficient(split_key[2]))

                #add the exchange reaction to the appropriate list in the dictionary of exchange reactions
                if split_key[2] in self.exch_sets:

                    #if here, add the new reaction to the list of reactions which are exchanges of that key in the dictionary
                    self.exch_sets.append(rxn)

                else:

                    #if here, then need to create a new dictionary item for a list of reactions
                    #list will have just one reaction at this time
                    self.exch_sets[split_key[2]] = [rxn]

            else:

                #state that the reactions is not an exchange
                setattr(rxn,"isexch",False)

                #give a blank for the exchanged metabolites
                setattr(rxn,"exchof","")

                #give a blank dictionary for the exchanged metabolites
                setattr(rxn,"exchstoich",0)

        #need to sprinkle these around whenever changing the model so changes stick correctly
        self.combined_model.solver.update()
        self.combined_model.repair()

        #also give origins and updated ids to metabolites
        for met in self.combined_model.metabolites:

            #add origin attribute, use the original model ID for the origin attribute
            setattr(met,'origin',model1.id)

            #add origin tag to the reaction id, metaid, and name sp don't get "ignoring reaction since it already exists" issue
            met.id = met.id + "_" + model1.id
            met.name = met.name + " " + model1.name

        #need to sprinkle these around whenever changing the model so changes stick correctly
        self.combined_model.solver.update()
        self.combined_model.repair()

        #add the second model
        self.combined_model.merge(copy.copy(model2))

        #need to sprinkle these around whenever changing the model so changes stick correctly
        self.combined_model.solver.update()
        self.combined_model.repair()

        #for each reaction now in the combined model
        for rxn in self.combined_model.reactions:

            #if it does not yet have an origin, then the origin is the second model
            if not hasattr(rxn,'origin'):

                #add origin attribute, use the original model ID for the origin attribute
                setattr(rxn,'origin',model2.id)

                #add origin tag to the reaction id, metaid, and name sp don't get "ignoring reaction since it already exists" issue
                rxn.id = rxn.id + "_" + model2.id
                rxn.name = rxn.name + " " + model2.name

            #while we are at it find exchange reactions, set attributes of if exchange reaction and exchanged metabolite
            #search using the regulat expression tag
            #not need to only make notes on exchange if it was part of the second model originally
            if bool((re.search(exch_tag,rxn.id))) and (bool(re.match(model2.id,rxn.origin))):

                #state that it is an exchange
                setattr(rxn,"isexch",True)

                #set metabolite which is exchanged
                rxn_mets = rxn.metabolites

                #set the metaboltes it is exchanging
                #note that this will be a dictionary
                
                #get the metabolite key
                met_key = str(rxn_mets.keys())

                #split the key to get just the metabolite ID
                split_key = re.findall(r'\w+', met_key)

                setattr(rxn,"exchof",split_key[2])
                setattr(rxn,"exchstoich",rxn.get_coefficient(split_key[2]))

                #add the exchange reaction to the appropriate list in the dictionary of exchange reactions
                if split_key[2] in self.exch_sets:

                    #if here, add the new reaction to the list of reactions which are exchanges of that key in the dictionary
                    self.exch_sets[split_key[2]].append(rxn)

                else:

                    #if here, then need to create a new dictionary item for a list of reactions
                    #list will have just one reaction at this time
                    self.exch_sets[split_key[2]] = [rxn]

            else:

                #state that the reactions is not an exchange
                setattr(rxn,"isexch",False)

                #give a blank for the exchanged metabolites
                setattr(rxn,"exchof","")

                #give a blank dictionary for the exchanged metabolites
                setattr(rxn,"exchstoich",0)

        #need to sprinkle these around whenever changing the model so changes stick correctly
        self.combined_model.solver.update()
        self.combined_model.repair()
        
        #also give origins and updated ids to metabolites
        for met in self.combined_model.metabolites:

            #if it does not yet have an origin, then the origin is the second model
            if not hasattr(met,'origin'):
                
                #add origin attribute, use the original model ID for the origin attribute
                setattr(met,'origin',model2.id)

                #add origin tag to the reaction id, metaid, and name sp don't get "ignoring reaction since it already exists" issue
                met.id = met.id + "_" + model2.id
                met.name = met.name + " " + model2.name

        #need to sprinkle these around whenever changing the model so changes stick correctly
        self.combined_model.solver.update()
        self.combined_model.repair()

        #set the name and ID of the combined model
        new_name = ""
        new_id = ""

        #make a count for the number of names, used for formatting
        num_names = 0
        
        #for each model
        for model in self.members:

            if num_names == 0: 

                #add that model to the community model name
                new_name = new_name+model.name
                new_id = new_id+model.id

            else: 

                #add that model to the community model name, with appropriate formatting
                new_name = new_name+"&"+model.name
                new_id = new_id+"&"+model.id
        
            num_names = num_names + 1
            
        self.combined_model.id = new_id+"_community"
        self.combined_model.name = new_name+"_community"

        #need to sprinkle these around whenever changing the model so changes stick correctly
        self.combined_model.solver.update()
        self.combined_model.repair()

        self.log.write("current model id: "+self.combined_model.id+"\n")
        self.log.write("current model name: "+self.combined_model.name+"\n")
        self.log.write("current model members: "+str(self.members)+"\n\n")

        self.log.write("current model exchange reaction sets: "+self.combined_model.id+"\n")
        #write the exchange dictionary to make sure I have done this right
        for key in self.exch_sets:

            self.log.write("key: "+key+"; value: "+str(self.exch_sets[key])+"\n")

        #need to sprinkle these around whenever changing the model so changes stick correctly
        self.combined_model.solver.update()
        self.combined_model.repair()

    #sets up the equations related to the medium based on the 
    def define_medium(self,media):

        """
        THIS SECTION DEALS WITH DEFINING u^c_i AND e^c_i VARIABLES AND THE CONSTRAINT RELATED TO CLASS-LEVEL 
        BECAUSE u^c_i AND e^c_i CAN TAKE ANY VALUE WITH RESPECT TO EACH OTHER WE WILL INSTEAD DEFINE A NEW VARIABLE
        x^c_i WHICH REPRESENTS THE COMMUNITY EXCHANGE OF METABOLITE i
        """

        #included for formatting of the log file
        self.log.write("\n\n")

        #we determine which metabolites need community constraints by self.exch_sets
        #recall that the keys of self.exch_sets are metabolites
        for met in self.exch_sets:

            #for each metabolite, we define one new variable and one new constraint to relate its exchanges to a community exchange
            #initiate x_c for the metabolite, default assumption it is not in the medium
            x_met = self.combined_model.problem.Variable(name='x_c_{}'.format(met),lb=-self.bigM,ub=0)

            #define x_c based on its presence/absence in the medium
            if met in media.keys():

                #first make the variable, this depends on if it is in the medium
                x_met.ub = media[met]

            #add the new variable to the model
            self.combined_model.add_cons_vars([x_met], sloppy=False)

            #apparently this is needed otherwise the new constraint won't register
            self.combined_model.solver.update()
            self.combined_model.repair()

            #now defined the new constraint
            exch_const = self.combined_model.problem.Constraint(x_met,lb=0,ub=0,name='exch_const_{}'.format(met),sloppy=False)

            #need to add constraint to the model before I can chang their coefficients
            self.combined_model.add_cons_vars([exch_const], sloppy=False)

            #apparently this is needed otherwise the new constraint won't register
            self.combined_model.solver.update()
            self.combined_model.repair()

            #need add each exchange reaction in the set to the constraint
            for exch_rxn in self.exch_sets[met]:

                #add flux expression of that reaction to the constraint
                exch_const.set_linear_coefficients({exch_rxn.forward_variable: 1 * self.X_k[exch_rxn.origin]})
                exch_const.set_linear_coefficients({exch_rxn.reverse_variable: -1 * self.X_k[exch_rxn.origin]})

            #by this point the exchange constraint should be written
            self.log.write("\nCommunity exchange constraint for "+met+":\n"+str(exch_const)+"\n")
            self.log.write("x_c bounds, lb: "+str(x_met.lb)+"\tub: "+str(x_met.ub)+"\n\n")

        #need to sprinkle these around whenever changing the model so changes stick correctly
        self.combined_model.solver.update()
        self.combined_model.repair()

    #this adds another member to the community
    #modeln - the model of the nth member to add to the community
    #exch_tag - as in __init__, the exchange tag is what 
    def add_member(self,modeln,exch_tag):

        self.members.append(modeln)

        #add the second model
        self.combined_model.merge(modeln.copy())

        #need to sprinkle these around whenever changing the model so changes stick correctly
        self.combined_model.solver.update()
        self.combined_model.repair()

        #for each reaction now in the combined model
        for rxn in self.combined_model.reactions:

            #if it does not yet have an origin, then the origin is the second model
            if not hasattr(rxn,'origin'):

                #add origin attribute, use the original model ID for the origin attribute
                setattr(rxn,'origin',modeln.id)

                #add origin tag to the reaction id, metaid, and name sp don't get "ignoring reaction since it already exists" issue
                rxn.id = rxn.id + "_" + modeln.id
                rxn.name = rxn.name + " " + modeln.name

            #while we are at it find exchange reactions, set attributes of if exchange reaction and exchanged metabolite
            #search using the regulat expression tag
            #not need to only make notes on exchange if it was part of the second model originally
            if (bool(re.search(exch_tag,rxn.id))) and (bool(re.match(modeln.id,rxn.origin))):

                #state that it is an exchange
                setattr(rxn,"isexch",True)

                #set metabolite which is exchanged
                rxn_mets = rxn.metabolites

                #set the metaboltes it is exchanging
                #note that this will be a dictionary
                
                #get the metabolite key
                met_key = str(rxn_mets.keys())

                #split the key to get just the metabolite ID
                split_key = re.findall(r'\w+', met_key)

                setattr(rxn,"exchof",split_key[2])
                setattr(rxn,"exchstoich",rxn.get_coefficient(split_key[2]))

                #add the exchange reaction to the appropriate list in the dictionary of exchange reactions
                if split_key[2] in self.exch_sets:

                    #if here, add the new reaction to the list of reactions which are exchanges of that key in the dictionary
                    self.exch_sets[split_key[2]].append(rxn)

                else:

                    #if here, then need to create a new dictionary item for a list of reactions
                    #list will have just one reaction at this time
                    self.exch_sets[split_key[2]] = [rxn]

            else:

                #state that the reactions is not an exchange
                setattr(rxn,"isexch",False)

                #give a blank for the exchanged metabolites
                setattr(rxn,"exchof","")

                #give a blank dictionary for the exchanged metabolites
                setattr(rxn,"exchstoich",0)

        #need to sprinkle these around whenever changing the model so changes stick correctly
        self.combined_model.solver.update()
        self.combined_model.repair()

        #also give origins and updated ids to metabolites
        for met in self.combined_model.metabolites:

            #if it does not yet have an origin, then the origin is the second model
            if not hasattr(met,'origin'):
                
                #add origin attribute, use the original model ID for the origin attribute
                setattr(met,'origin',modeln.id)

                #add origin tag to the reaction id, metaid, and name sp don't get "ignoring reaction since it already exists" issue
                met.id = met.id + "_" + modeln.id
                met.name = met.name + " " + modeln.name

        #need to sprinkle these around whenever changing the model so changes stick correctly
        self.combined_model.solver.update()
        self.combined_model.repair()

        #set the name and ID of the combined model
        new_name = ""
        new_id = ""

        #make a count for the number of names, used for formatting
        num_names = 0
        
        #for each model
        for model in self.members:

            if num_names == 0: 

                #add that model to the community model name
                new_name = new_name+model.name
                new_id = new_id+model.id

            else: 

                #add that model to the community model name, with appropriate formatting
                new_name = new_name+"&"+model.name
                new_id = new_id+"&"+model.id
        
            num_names = num_names + 1
            
        self.combined_model.id = new_id+"_community"
        self.combined_model.name = new_name+"_community"

        #need to sprinkle these around whenever changing the model so changes stick correctly
        self.combined_model.solver.update()
        self.combined_model.repair()

        self.log.write("current model id: "+self.combined_model.id)
        self.log.write("current model name: "+self.combined_model.name)
        self.log.write("current model members: "+str(self.members)+"\n")

        self.log.write("current model exchange reaction sets: "+self.combined_model.id)
        #write the exchange dictionary to make sure I have done this right
        for key in self.exch_sets:

            self.log.write("key: "+key+"; value: "+str(self.exch_sets[key]))

    #method to define the abundances of the species members
    #note - length of X_k needs to be the same as the number of members of the community for the assignement to work
    #keys should be model ID's, values abundances
    #also - abundances should add up to one
    #returns false if at least one of the above is violated, true if assignment was made
    def define_abundance(self,X_k):
        
        #get the sum of abundances
        sum_abund = 0

        for key in X_k:

            sum_abund = sum_abund + X_k[key]
    
        #check if abundanecs add up to one and have the same number in the list as members
        if (len(X_k) == len(self.members)) and (sum_abund == 1):

            #if so, make the abundances assignement
            self.X_k = X_k
            return True

        else:

            print("Wrong number of abundances given")
            self.log.write("Wrong number of abundances given")
            return False
        
    #self - needs to be passed itself
    #x - mass fraction of total community mass which is of the species of model 1. This is necessary because
    #biomass_dict - dictionary of biomass reaction identifiers
    #all of the flux rates are normalize by gDW of the modeled species, so to have community exchanges will need
    #to convert to the basis of the weight of the total community
    #note that the created reactions for the community lack metabolites, but are tied directy to the exchanges of the community members
    #allows tracking community-based exchange rates without mass imbalances that might get flagged because "gDW" is different for each
    #community member and the community as a whole
    #note this is pretty much SteadyCom 
    def build_comm_x(self,biomass_dict): 

        #check if abundances add up to one and have the same number in the list as members
        if not len(biomass_dict) == len(self.members):

            print("\n\nWrong number of biomass equations given!")
            self.log.write("\n\nWrong number of biomass equations given!")
            return False 

        #keep the biomass_dict for later use in other functions
        self.biomass_dict = biomass_dict

        #if here, then the right number of biomass equations given

        self.log.write("\n\nBegin log for building the following community model of with constant community composition: \n")

        #since variable number of models, need to next writing these models and abundanes in a loop
        num_models = 0

        for model in self.members:

            #write what the component models are abundances are to the log file
            self.log.write("Model "+str(num_models)+": "+model.id+", (abundance: "+str(self.X_k[model.id])+")\n")

            num_models = num_models + 1

        self.log.write("Community Model: "+self.combined_model.id+"\n")

        """
        THIS SECTION DEALS WITH CHANGING v^k_j TO V^k_j (E.G. SCALING EVERYTHING BY X^K)
        THIS CAN BE DONE MOST SIMPLY BY MULTIPYING ALL COEFFICIENTS OF v^k_j BY X^K AND SIMPLY NOTING THAT v^k_j IS NOW V^k_j
        FOR EACH REAECTION IN THE MODEL
        """

        #number of metabolite constraints changed
        num_mets_done = 0

        #stores the number of metabolites in the model
        num_mets = len(self.combined_model.metabolites)

        #determines if the first time above a certain percentage
        done_10 = False
        done_20 = False
        done_30 = False
        done_40 = False
        done_50 = False
        done_60 = False
        done_70 = False
        done_80 = False
        done_90 = False
        done_100 = False

        #note that the [1:2] bit means that the list has only one element, big time saver when debugging
        #for met in self.combined_model.metabolites[1:2]:
        for met in self.combined_model.metabolites:

            #get the constraints it is involved with
            constraint = met.constraint
            
            #get the constraint expression
            const_expr = constraint.expression

            #use the expression, split into terms
            #note that pulling coefficients using "get_linear_coefficients" always showed coefficients of zero, so have to go through expression to get 
            #them using find all since there should be multiple matches
            #this captures sign, stoichiometry, and reaction in order, ignores the rest of the framework
            
            self.log.write("\nconstraint to update\n")
            self.log.write("constraint: "+str(const_expr)+"\n")

            expr_terms = re.findall(r"(?P<sign>\-|\+)*\s*(?P<stoich>\d+\.\d+)\*(?P<rxn>.+?)(\s|$)", str(const_expr))

            self.log.write("split constraint: "+str(expr_terms)+"\n")

            #loop through each term getting stoichiometry and reaction
            for term in expr_terms:

                #multiply abundance by old stoichiometry to get the new stoichiometry
                new_stoich = float(term[1])*self.X_k[met.origin]

                #get the variable
                #specify an empty string as a placeholder
                term_var = ""

                #for each variable
                for var in constraint.variables:

                    #check if name matches the term[2] value
                    if bool(re.match("^"+term[2]+"$",var.name)):

                        #then var is the variable we are looking for
                        term_var = var

                #aupdate constraint coefficient for the new term
                #needed to concatenate the sign with the coefficient when passing it to update the stoichiometry
                constraint.set_linear_coefficients({term_var: float(str(term[0])+str(new_stoich))})  

                #update the bounds of the constraint based on species abundance
                constraint.lb = constraint.lb * self.X_k[met.origin]
                constraint.ub = constraint.ub * self.X_k[met.origin]

            #after this, the constraint in question should be updated and reaction rates will now be V^k_j
            #here, we write those results to the log for inspection in case that is necesssary 
            self.log.write("new mass balance constraint on: "+met.id+"\n"+str(constraint)+"\n\n")

            #decide if need to report on progress
            num_mets_done += 1

            if num_mets_done >= 0.1 * num_mets and not done_10:

                done_10 = True

                print("10% complete")

            elif num_mets_done >= 0.2 * num_mets and not done_20:

                done_20 = True

                print("20% complete")

            elif num_mets_done >= 0.3 * num_mets and not done_30:

                done_30 = True

                print("30% complete")

            elif num_mets_done >= 0.4 * num_mets and not done_40:

                done_40 = True

                print("40% complete")    

            elif num_mets_done >= 0.5 * num_mets and not done_50:

                done_50 = True

                print("50% complete")    

            elif num_mets_done >= 0.6 * num_mets and not done_60:

                done_60 = True

                print("60% complete")    

            elif num_mets_done >= 0.7 * num_mets and not done_70:

                done_70 = True

                print("70% complete")    

            elif num_mets_done >= 0.8 * num_mets and not done_80:

                done_80 = True

                print("80% complete")    

            elif num_mets_done >= 0.9 * num_mets and not done_90:

                done_90 = True

                print("90% complete")    

            elif num_mets_done >= 1 * num_mets and not done_100:

                done_100 = True

                print("100% complete")      

        #need to sprinkle these around whenever changing the model so changes stick correctly
        self.combined_model.solver.update()
        self.combined_model.repair()
        
        """
        THIS SECTION DEALS WITH DEFINING A BIOMASS CONSTRAINT
        Note that, we have already set V^K = X^k * mu for each species, but mu will be a variable. So here we add the variable mu
        and link each biomass to mu
        """
        
        #first create the variable mu
        mu_var = self.combined_model.problem.Variable('mu',lb=0,ub=self.bigM)

        #we will set mu 
        mu_var.set_bounds(0,self.bigM)

        #add the new variable to the model
        self.combined_model.add_cons_vars([mu_var], sloppy=False)

        #apparently this is needed otherwise the new constraint won't register
        self.combined_model.solver.update()

        #for each model reaction
        for rxn in self.combined_model.reactions:

            #check if it is a biomass reaction, recall added 
            if bool(re.fullmatch(rxn.id,biomass_dict[rxn.origin]+"_"+rxn.origin)):

                #define a new constraint based on biomass
                bio_const = self.combined_model.problem.Constraint(rxn.flux_expression - self.X_k[rxn.origin] * mu_var,lb=0,ub=0,sloppy=False,name='bio_const_{}'.format(rxn.origin))

                #add the new constraint to the model
                self.combined_model.add_cons_vars([bio_const], sloppy=False)

                #apparently this is needed otherwise the new constraint won't register
                self.combined_model.solver.update()

                #by this point the exchange constraint should be written
                self.log.write("Biomass constraint for "+rxn.id+" (model: "+rxn.origin+"):\n"+str(bio_const)+"\n\n")

        #need to sprinkle these around whenever changing the model so changes stick correctly
        self.combined_model.solver.update()
        self.combined_model.repair()
        
        """
        THIS SECTION SETS MAXIMIZING BIOMASS AS THE OBJECTIVE OF THE MODEL
        """

        #start by creating a new objective equations
        self.combined_model.objective = self.combined_model.problem.Objective(mu_var, direction='max')
        
        #I think this is all that is needed, lets check
        self.log.write("Ojective equation for "+self.combined_model.id+":\n"+str(self.combined_model.objective)+"\n")

        #need to sprinkle these around whenever changing the model so changes stick correctly
        self.combined_model.solver.update()
        self.combined_model.repair()

    #pass a string to set the solver to that string
    def set_solver(self,solver):

        #set the solver to the passed string
        try:

            self.model.solver = solver

        #if an exception occurs, store as "e"
        except Exception as e:

            print("solver assignement unsuccessful, exception: "+str(e))

            #if here, it didn't work, return false
            return False
        
        #need to sprinkle these around whenever changing the model so changes stick correctly
        self.combined_model.solver.update()
        self.combined_model.repair()

        #if here then it worked, return true
        return True

    #the job of this method is to find the maximum growth rate (mu) which the model can achieve
    #everything is set up already, so just need to solve
    #media - a dictionary of metabolites which comprises allowed community uptake metabolites and the max uptake rate
    def max_mu(self,fixed_rates=dict()):

        #create a duplicate model for adding constraints without affecting the base model
        max_mu_model = self.combined_model.copy()

        #give the problem a name
        max_mu_model.problem.name = "Find parsimonious maximum growth sum"

        #try to make sure the model is good to go for solving
        max_mu_model.solver.update()
        max_mu_model.repair()

        #initialize an empty dictionary for returning with results
        mu_results = { }

        #add to the results the list of exchange sets
        mu_results['ex_sets'] = self.exch_sets
    
        #at this point, everything should be set up to maximize for mu
        self.log.write("\n\nAttempting to solve "+max_mu_model.id+" for maximum growth rate\n\n")

        #fix the rates that need to be fixed, if any
        for rxn in max_mu_model.reactions:

            if rxn.id in fixed_rates.keys():

                rxn.lower_bound = fixed_rates[rxn.id]
                rxn.upper_bound = fixed_rates[rxn.id]

        #solve, but put in a try/except framework in case there is an error
        try:

            start_time_mu = datetime.now()

            #need to sprinkle these around whenever changing the model so changes stick correctly
            max_mu_model.solver.update()
            max_mu_model.repair()

            #try to solve, here is where the error may get thrown
            #this will be to maximize mu
            mu_soln = max_mu_model.optimize()

            #get the sum of biomass reaction rates so we can fix them
            max_mu = mu_soln.objective_value

            mu_results['mu_objective'] = max_mu

            print("solver status: \n"+str(mu_soln.status)+"\n")
            print("Objective value (mu): \n"+str(max_mu)+"\n\n")
            self.log.write("solver status: \n"+str(mu_soln.status)+"\n")
            self.log.write("Objective value (mu): \n"+str(max_mu)+"\n\n")

            #fix the value of mu based on this solution so that biomass rates must be maintained while minimizing reaction rates

            #do this by fixing the bounds
            max_mu_model.variables.mu.lb = max_mu
            max_mu_model.variables.mu.ub = max_mu

            #create a new objective equation minimizing the sum of flux rates
            #set a dummy objective to add coefficients for each reaction to
            max_mu_model.objective = max_mu_model.problem.Objective(Zero, direction='min')

            #go through each reaction, see which matches the identifier
            for rxn in max_mu_model.reactions:

                #make a variable to store the absolute value of each reaction rate
                v_plus_rxn = max_mu_model.problem.Variable(name='v_+_{}'.format(rxn.id),lb=0,ub=2*self.bigM)

                #create two constraints to get back the absolute value
                v_plus_const_1 = max_mu_model.problem.Constraint(v_plus_rxn - rxn.flux_expression,lb=0,ub=2*self.bigM,name='v_+_1_{}'.format(rxn.id),sloppy=False)
                v_plus_const_2 = max_mu_model.problem.Constraint(v_plus_rxn + rxn.flux_expression,lb=0,ub=2*self.bigM,name='v_+_2_{}'.format(rxn.id),sloppy=False)

                #add these constraints to the model
                max_mu_model.add_cons_vars([v_plus_const_1, v_plus_const_2], sloppy=False)

                #set the objective so that we are minimizing the sum of absolute values
                max_mu_model.objective.set_linear_coefficients({v_plus_rxn: 1})

            #need to sprinkle these around whenever changing the model so changes stick correctly
            max_mu_model.solver.update()
            max_mu_model.repair()

            #solve with a fixed growth rate, minimizing sum of reaction fluxes
            mu_soln = max_mu_model.optimize()

            print("solver status: \n"+str(mu_soln.status)+"\n")
            print("Objective value (mu): \n"+str(mu_results['mu_objective'])+"\n")
            print("Objective value (flux sum): \n"+str(mu_soln.objective_value)+"\n\n")
            self.log.write("solver status: \n"+str(mu_soln.status)+"\n")
            self.log.write("Objective value (mu): \n"+str(mu_results['mu_objective'])+"\n")
            self.log.write("Objective value (flux sum): \n"+str(mu_soln.objective_value)+"\n\n")

            end_time_mu = datetime.now()

            #get the total solve time
            total_time_mu = end_time_mu - start_time_mu

            mu_results['solve_time'] = total_time_mu

            #state that no exception occured
            mu_results['exception']=False

            mu_results['status'] = mu_soln.status

            #return the solution time in the dictionary
            mu_results['soln_time'] = str(total_time_mu)

            if mu_soln.status == INFEASIBLE:

                #no objective to return
                mu_results['mu_objective'] = 0
                mu_results['flux_objective'] = 0

                #write and store the lower bound, flux, and upper bound for each reaction
                for rxn in max_mu_model.reactions:

                    #initialize element to nest
                    mu_results[rxn.id] = { }
                
                    mu_results[rxn.id]['lb'] = rxn.lower_bound
                    mu_results[rxn.id]['flux'] = 0
                    mu_results[rxn.id]['ub'] = rxn.upper_bound

                #initialize x_c dictionary
                mu_results['x_c'] = { }
                
                #save the community exchange reactions
                for met in self.exch_sets.keys():
                    
                    mu_results['x_c'][met] = 0

            else:

                #return the objective
                mu_results['mu_objective'] = max_mu
                mu_results['flux_objective'] = mu_soln.objective_value

                #store the lower bound, flux, and upper bound for each reaction
                for rxn in max_mu_model.reactions:

                    #store the results to the output dictionary
                    #I believe this would be most useful set up as a nested dictionary

                    #initialize element to nest
                    mu_results[rxn.id] = { }

                    #nest bounds and flux rates
                    mu_results[rxn.id]['lb'] = rxn.lower_bound
                    mu_results[rxn.id]['flux'] = mu_soln.fluxes[rxn.id]
                    mu_results[rxn.id]['ub'] = rxn.upper_bound

                #initialize x_c dictionary
                mu_results['x_c'] = { }

                #save the community exchange reactions
                for met in self.exch_sets.keys():
                    
                    mu_results['x_c'][met] = max_mu_model.solver.primal_values.get('x_c_{}'.format(met))

        #if an exception occurs, store as "e"
        except Exception as e:

            #get the timein information
            end_time_mu = datetime.now()

            #print the total solve time
            total_time_mu = end_time_mu - start_time_mu
            
            #state that an exception occured
            mu_results['exception']=True

            mu_results['status'] = "exception occurred"

            #save the exception string to return
            mu_results['exception_str']=str(e)

            #return the solution time in the dictionary
            mu_results['soln_time'] = str(total_time_mu)
            mu_results['solve_time'] = str(total_time_mu)

            #return objective value of NaN since the problem was not solved
            mu_results['mu_objective'] = 0
            mu_results['flux_objective'] = 0

            #write and store the lower bound, flux, and upper bound for each reaction
            for rxn in max_mu_model.reactions:

                #initialize element to nest
                mu_results[rxn.id] = { }

                mu_results[rxn.id]['lb'] = rxn.lower_bound
                mu_results[rxn.id]['flux'] = 0
                mu_results[rxn.id]['ub'] = rxn.upper_bound

            #initialize x_c dictionary
            mu_results['x_c'] = { }
            
            #save the community exchange reactions
            for met in self.exch_sets.keys():
                    
                mu_results['x_c'][met] = 0
        
        #return the solution that it got
        return mu_results

    #this function will seek to maximize the sum of species biomasses
    #will do this on a copy of the combined mode to avoid messing up the combined
    #model. Need to pass in a dictionary of biomass equations
    #media - an array of metabolite ids which are allowed to be uptaken by the community
    def max_sum(self,biomass_dict):

        """
        This section deals with initial checks and setting the biomass sum as the objective equation        
        """

        #check if abundances add up to one and have the same number in the list as members
        if not len(biomass_dict) == len(self.members):

            print("\n\nWrong number of biomass equations given!")
            self.log.write("\n\nWrong number of biomass equations given!")
            return False 

        max_sum_model = self.combined_model.copy()

        #initialize an empty dictionary for returning with results
        max_results = { }

        #add to the results the list of exchange sets
        max_results['ex_sets'] = self.exch_sets

        #need to sprinkle these around whenever changing the model so changes stick correctly
        max_sum_model.solver.update()
        max_sum_model.repair()

        #set a dummy objective to add biomass equations to
        max_sum_model.objective = max_sum_model.problem.Objective(Zero, direction='max')

        #go through each reaction, see which matches the identifier
        for rxn in max_sum_model.reactions:

            for model in biomass_dict:

                if bool(re.fullmatch(str(rxn.id),biomass_dict[model]+"_"+model)):

                    #set the linear coefficient
                    max_sum_model.objective.set_linear_coefficients({rxn.forward_variable: 1})
                    max_sum_model.objective.set_linear_coefficients({rxn.reverse_variable: -1})
        
        #try solving
        #solve, but put in a try/except framework in case there is an error
        try:

            start_time_max = datetime.now()

            #need to sprinkle these around whenever changing the model so changes stick correctly
            max_sum_model.solver.update()
            max_sum_model.repair()
            
            #need to sprinkle these around whenever changing the model so changes stick correctly
            max_sum_model.solver.update()
            max_sum_model.repair()
            
            #solve in a pfba-like manner, start by finding the maximum value of mu
            max_soln = max_sum_model.optimize()

            #get the sum of biomass reaction rates so we can fix them
            max_bio_sum = max_soln.objective_value

            print("solver status: \n"+str(max_soln.status)+"\n")
            print("Objective value (bio sum): \n"+str(max_bio_sum)+"\n")
            self.log.write("solver status: \n"+str(max_soln.status)+"\n")
            self.log.write("Objective value (bio sum): \n"+str(max_bio_sum)+"\n")
            self.log.write("Individual biomass flux rates:\n")
            
            #create a constraint to ensure the next solution has the same sum of biomass rates
            bio_sum_const = max_sum_model.problem.Constraint(Zero,lb=max_bio_sum,ub=max_bio_sum,name='bio_sum_const',sloppy=False)

            #add the new constraint to the model so we can change its coefficients
            max_sum_model.add_cons_vars(bio_sum_const)

            #need to sprinkle these around whenever changing the model so changes stick correctly
            max_sum_model.solver.update()
            max_sum_model.repair()

            #find the biomass equations, give them a "-1" coefficient
            for rxn in max_sum_model.reactions:

                for model in biomass_dict:

                    if bool(re.fullmatch(str(rxn.id),biomass_dict[model]+"_"+model)):

                        #set the linear coefficient
                        bio_sum_const.set_linear_coefficients({rxn.forward_variable: 1})
                        bio_sum_const.set_linear_coefficients({rxn.reverse_variable: -1})

            self.log.write("biomass constraint for second solve: "+str(bio_sum_const)+"\n\n")

            #need to sprinkle these around whenever changing the model so changes stick correctly
            max_sum_model.solver.update()
            max_sum_model.repair()

            #create a new objective equation minimizing the sum of flux rates
            #set a dummy objective to add coefficients for each reaction to
            max_sum_model.objective = max_sum_model.problem.Objective(Zero, direction='min')
            
            #go through each reaction, see which matches the identifier
            for rxn in max_sum_model.reactions:

                #make a variable to store the absolute value of each reaction rate
                #note that the max needs to be 2*bigM in case any reactions are at bigM for rate
                v_plus_rxn = max_sum_model.problem.Variable(name='v_+_{}'.format(rxn.id),lb=0,ub=2*self.bigM)

                #create two constraints to get back the absolute value
                v_plus_const_1 = max_sum_model.problem.Constraint(v_plus_rxn - rxn.flux_expression,lb=0,ub=2*self.bigM,name='v_+_1_{}'.format(rxn.id),sloppy=False)
                v_plus_const_2 = max_sum_model.problem.Constraint(v_plus_rxn + rxn.flux_expression,lb=0,ub=2*self.bigM,name='v_+_2_{}'.format(rxn.id),sloppy=False)

                #add these constraints to the model
                max_sum_model.add_cons_vars([v_plus_const_1, v_plus_const_2], sloppy=False)

                #set the objective so that we are minimizing the sum of absolute values
                max_sum_model.objective.set_linear_coefficients({v_plus_rxn: 1})

                #need to sprinkle these around whenever changing the model so changes stick correctly
                max_sum_model.solver.update()
                max_sum_model.repair()

            #need to sprinkle these around whenever changing the model so changes stick correctly
            max_sum_model.solver.update()
            max_sum_model.repair()

            #solve with a fixed growth rate, minimizing sum of reaction fluxes
            max_soln = max_sum_model.optimize()

            #go through each reaction, get a sum for the biomass reaction rates
            bio_sum_2 = 0

            for rxn in max_sum_model.reactions:

                for model in biomass_dict:

                    if bool(re.fullmatch(str(rxn.id),biomass_dict[model]+"_"+model)):

                        #write the biomass flux rate
                        bio_sum_2 = bio_sum_2 + max_soln.fluxes[rxn.id]

            print("solver status: \n"+str(max_soln.status)+"\n")
            print("Objective value (bio sum): \n"+str(bio_sum_2)+"\n\n")
            print("Objective value (flux sum): \n"+str(max_soln.objective_value)+"\n")
            self.log.write("solver status: \n"+str(max_soln.status)+"\n")
            self.log.write("Objective value (bio sum): \n"+str(bio_sum_2)+"\n")
            self.log.write("Objective value (flux sum): \n"+str(max_soln.objective_value)+"\n\n")

            end_time_max = datetime.now()

            #get the total solve time
            total_time_max = end_time_max - start_time_max

            max_results['solve_time'] = total_time_max

            #state that no exception occured
            max_results['exception']=False

            max_results['status'] = max_soln.status

            #return the solution time in the dictionary
            max_results['soln_time'] = str(total_time_max)

            if max_soln.status == INFEASIBLE:

                #no objective to return
                max_results['flux_objective'] = 0
                max_results['bio_objective'] = 0

                #write and store the lower bound, flux, and upper bound for each reaction
                for rxn in max_sum_model.reactions:

                    #initialize element to nest
                    max_results[rxn.id] = { }
                
                    max_results[rxn.id]['lb'] = rxn.lower_bound
                    max_results[rxn.id]['flux'] = 0
                    max_results[rxn.id]['ub'] = rxn.upper_bound

                #initialize x_c dictionary
                max_results['x_c'] = { }
                
                #save the community exchange reactions dummy values
                for met in self.exch_sets.keys():
                    
                    max_results['x_c'][met] = 0

            else:

                #return the objective
                max_results['flux_objective'] = max_soln.objective_value
                max_results['bio_objective'] = max_bio_sum

                #store the lower bound, flux, and upper bound for each reaction
                for rxn in max_sum_model.reactions:

                    #store the results to the output dictionary
                    #I believe this would be most useful set up as a nested dictionary

                    #initialize element to nest
                    max_results[rxn.id] = { }

                    #nest bounds and flux rates
                    max_results[rxn.id]['lb'] = rxn.lower_bound
                    max_results[rxn.id]['flux'] = max_soln.fluxes[rxn.id]
                    max_results[rxn.id]['ub'] = rxn.upper_bound

                #initialize x_c dictionary
                max_results['x_c'] = { }
                
                #save the community exchange reactions
                for met in self.exch_sets.keys():
                    
                    max_results['x_c'][met] = max_sum_model.solver.primal_values.get('x_c_{}'.format(met))

        #if an exception occurs, store as "e"
        except Exception as e:

            #get the timein information
            end_time_max = datetime.now()

            #print the total solve time
            total_time_max = end_time_max - start_time_max
            
            #state that an exception occured
            max_results['exception']=True

            max_results['status'] = "exception occurred"

            #save the exception string to return
            max_results['exception_str']=str(e)

            #return the solution time in the dictionary
            max_results['soln_time'] = str(total_time_max)

            max_results['solve_time'] = total_time_max

            #return objective value of NaN since the problem was not solved
            max_results['bio_objective'] = 0
            max_results['flux_objective'] = 0

            #write and store the lower bound, flux, and upper bound for each reaction
            for rxn in max_sum_model.reactions:

                #initialize element to nest
                max_results[rxn.id] = { }

                max_results[rxn.id]['lb'] = rxn.lower_bound
                max_results[rxn.id]['flux'] = "NaN"
                max_results[rxn.id]['ub'] = rxn.upper_bound

            #initialize x_c dictionary
            max_results['x_c'] = { }
            
            #save the community exchange reactions dummy values
            for met in self.exch_sets:
            
                max_results['x_c'][met] = "NaN"
        
        #return the solution that it got
        return max_results
