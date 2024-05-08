# -*- coding: utf-8 -*-

from cobra.util.solver import linear_reaction_coefficients
import json 
import multiprocessing
from multiprocessing import Pipe, Process
import time 

from mqc.config import Config
from mqc.utils import *

class Quantitative():
    """ """
    def __init__(self,cfg:Config) -> None:
        """ """
        self.cfg = cfg

    def find_net_initial_flux(self, model_info, model, model_control_info):
        """
        Determining whether there is a net formation of substances
        """
        model_control_info["check_metabolite_production"]["score"] = 1
        r_met, net_generation, only_pro_mets, met_num = [], [], [], 0
        initial_net_flux, temp = [], {}
        for rxn in model_info["reactions"]:  # Find the metabolites of the reactions with the lowest bound less than 0 in the exchange reaction
            if rxn["id"] in model_info["exchange_rxns"] and rxn["bounds"][0] < 0:
                r_met += rxn["all_mets"]
            if len(rxn['reactants_mets']) == 0 and len(rxn['products_mets']) != 0:
                only_pro_mets += rxn['products_mets']
        for met in model.metabolites:
            met_num += 1
            if met.name not in (FREE_METS + only_pro_mets):
                if re.search("C[A-Z]",str(met.formula)) or re.search("C[\d]",str(met.formula)): # Carbon followed by capital letters or numbers         
                    with model:
                        if met.name in ACCOA and met.compartment in C_COMPARTMENT:
                            objectiveId = add_accoa(model, str(met.id))
                        else:
                            objectiveId = add_demand(model_info, model, met.id)
                        model.objective = objectiveId
                        if model.slim_optimize() > 1e-5 and met.name not in r_met:
                            net_generation.append(met.id)
                            initial_net_flux.extend([f"{met.id} : {round(model.slim_optimize(),3)}"])
                            temp[f"{met.id}"] = round(model.slim_optimize(),3)
        model_info['net_generation'] = net_generation
        if len(net_generation) != 0:
            model_control_info["check_metabolite_production"]["score"] = 0
        model_info["total_number_of_metabolites:"] = met_num
        model_info["total_net_matter:"] = len(net_generation)
        model_info["net_matter:"] = ', '.join(net_generation)
        # if len(initial_net_flux) != 0:
        #     model_control_info['initial_flux']["initial_net"] = ', '.join(initial_net_flux)
        # else:
        #     model_control_info['initial_flux']["initial_net"] = "no net matter"
        # model_control_info['initial_flux']["initial_net"] = temp
        model_info['initial_flux']["initial_net"] = temp


    def get_energy_initial_flux(self, model_info, model, model_control_info):
        """"""
        xtp, xtp_flux = ['A','C','G','U','I'], []
        for k in xtp:
            if f"ADD_{k}TPM" in model.reactions:
                atpm_id = "ADD_ATPM"
            else:
                atpm_id = get_atpm_rxn(model_info, model, f'{k}TP', f'{k}DP')
            model.reactions.get_by_id(atpm_id).bounds = (0,1000)
            set_model_objective(model_info, model, atpm_id)
            model_info['initial_flux'][f"initial_{k}tp"] = model_info["now_obj_rxn"]
            # model_control_info['initial_flux'][f"initial_{k}tp"] = round(model_info['now_obj_rxn']['rxn_flux'],3)
            xtp_flux.extend([f"{k}TP : {round(model_info['now_obj_rxn']['rxn_flux'],3)}"])
        # model_control_info['initial_flux']["initial_atp"] = ', '.join(xtp_flux)


    def check_yield_formula(self, model, metId):
        """"""
        try:
            met = model.metabolites.get_by_id(metId)
            met.elements.keys()
        except:
            print(met.id,'...................................')
            return 0
        calculate_met = ['C','H','O','N','P','S']
        met = model.metabolites.get_by_id(metId)
        if not set(met.elements.keys()).issubset(set(calculate_met)): # formula含有[CHONPS]以外的不算
            return 0
        if pd.isna(met.formula) or 'R' in met.formula or 'X' in met.formula or met.formula == 'nan' or not met.formula or met.formula == 'null':
            return 0
        return 1

    def find_yield_generation_subprocess(self, model_info, model, metabolites,conn):
        """Get all carbonaceous species that exceed the maximum theoretical yield"""
        yield_generation, c_num = [], 0
        initial_yield_flux, temp = [], {}
        
        pid = multiprocessing.current_process().pid
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),f'pid {pid} model.metabolites num:',len(metabolites))
        
        for met in metabolites:
            if met.id in model_info['need_set_carbon_mets']:
                continue
            if re.search("C[A-Z]",str(met.formula)) or re.search("C[\d]",str(met.formula)): # Carbon followed by capital letters or numbers
                c_num += 1
                if self.check_yield_formula(model, met.id) == 0:
                    print(met,': formula不明确,不计算')
                    continue
                try:
                    with model:
                        if met.name in ACCOA and met.compartment in C_COMPARTMENT:
                            objectiveId = add_accoa(model, str(met.id))
                        else:
                            objectiveId = add_demand(model_info, model, met.id)
                        model.objective = objectiveId
                        # write_flux_file(model_info, model, objectiveId,self.cfg)
                        max_rate = max_reduced_degree_rate(model, met.id)
                        if max_rate == 0: # 含碳物质，最大得率如果是0不用算；formula含有[CHONPS]以外的不算
                            continue
                        
                        model_opt = model.slim_optimize()
                        print(f'pid {pid} met: ',met.id,' object: ',model_opt,'  max_rate: ',max_rate)
                        if round(model_opt,2) > round(max_rate,2):
                            yield_generation.append(met.id)
                            initial_yield_flux.extend([f"{met.id} : {round(model_opt,3)}"])
                            temp[f"{met.id}"] = round(model_opt,3)
                except:
                    pass 
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),f' pid {pid} yield_generation:',yield_generation)
        
        conn.send({
            "yield_generation":yield_generation,
            "c_num":c_num,
            "initial_yield_flux":initial_yield_flux,
            "met_flux":temp,
        })
        conn.close()
        

    def find_yield_initial_flux(self,model_info, model, model_control_info):
        """
        
        """
        print ("Running in parallel")

        # get all model metabolites
        metabolites = model.metabolites
        print("metabolites num: ",len(metabolites))
        metabolites_num = len(metabolites)
        PROCESS_NUM = 10
        
        average = metabolites_num // 10
        remains = metabolites_num % 10
        # create a list to keep all processes
        processes = []

        # create a list to keep connections
        parent_connections = []
        
        # create a process per instance
        for i in range(PROCESS_NUM):
            if i == PROCESS_NUM -1:
                tmp_metabolites = metabolites[i*average:]
            else:
                tmp_metabolites =  metabolites[i*average:(i+1)*average]     
            # create a pipe for communication
            parent_conn, child_conn = Pipe()
            parent_connections.append(parent_conn)

            # create the process, pass instance and connection
            process = Process(target=self.find_yield_generation_subprocess, args=(model_info,model,tmp_metabolites, child_conn,))
            processes.append(process)

        # start all processes
        for process in processes:
            process.start()

        # make sure that all processes have finished
        for process in processes:
            process.join()

        yield_generation = []
        c_num = 0
        initial_yield_flux = []
        temp_met_flux = {}
        for parent_connection in parent_connections:
            response = parent_connection.recv()
            yield_generation += response["yield_generation"]
            c_num += response["c_num"]
            initial_yield_flux += response["initial_yield_flux"]
            temp_met_flux.update(response["met_flux"])
            
            
        model_info['yield_generation'] = yield_generation
        if len(yield_generation) != 0:
            model_control_info["check_metabolite_yield"]["score"] = 0
        model_info["total_carbon_content"] = c_num
        model_info["yield_exceeds_the_theoretical_maximum_number"] = len(yield_generation)
        model_info["yield_net_matter"] = ', '.join(yield_generation)
        # if len(initial_yield_flux) != 0:
        #     model_control_info['initial_flux']["initial_yield"] = ', '.join(initial_yield_flux)
        # else:
        #     model_control_info['initial_flux']["initial_yield"] = 'no yield net substance'
        # model_control_info['initial_flux']["initial_yield"] = temp_met_flux
        


    def get_model_initial_flux(self, model_info, model, model_control_info):
        """"""
        with model:
            model_info['initial_flux'] = {}
            model_info['initial_flux']["initial_biomass"] = model_info["initial_rxn"]
            # model_control_info['initial_flux']["initial_biomass"] = round(model_info["initial_rxn"]["biomass_rxn_flux"],3)
            if "ADD_NADH" in model.reactions:
                nadh_id = "ADD_NADH"
            else:
                nadh_id = add_nadh_rxn(model)
            model.reactions.get_by_id(nadh_id).bounds = (0,1000)
            set_model_objective(model_info, model, nadh_id)    
            model_info['initial_flux']["initial_nadh"] = model_info["now_obj_rxn"]
            # model_control_info['initial_flux']["initial_nadh"] = round(model_info["now_obj_rxn"]["rxn_flux"],3)
            # model.reactions.remove('ADD_NADH')
            self.get_energy_initial_flux(model_info, model, model_control_info)
            # self.find_net_initial_flux(model_info, model, model_control_info)
            model_control_info["check_metabolite_yield"]["score"] = 1
            # model_control_info['initial_flux']["initial_yield"] = ""
            if is_autotrophic(model_info):
                model_control_info["check_metabolite_yield"]["autotrophic"] = "This is an autotrophic strain, no need to calculate the yield"
                # model_control_info['initial_flux']["initial_yield"] = "This is an autotrophic strain, no need to calculate the yield"
                return 0
            else:
                get_c_source_supply(model_info, model)
                self.find_yield_initial_flux(model_info, model, model_control_info)

    def get_initial(self, model_info, model, model_control_info):
        """""" 
        # model_control_info['initial_flux'] = {}
        self.get_model_initial_flux(model_info, model, model_control_info)
        