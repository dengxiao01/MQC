# -*- coding: utf-8 -*-

from cobra.util.solver import linear_reaction_coefficients
import json 
import multiprocessing
from multiprocessing import Pipe, Process
import time 

from mqc.config import Config
from mqc.control.rules_control import Rules
from mqc.utils import *


class InitialPreprocess():
    """ """
    def __init__(self,cfg:Config) -> None:
        """ """
        self.cfg = cfg
    

    def get_exchange(self, model, model_info, check_model):
        """"""
        temp, exchange_info = [], {}
        modify_rxn = list(set(model_info['modify_rxn']))
        for ids in modify_rxn:
            rxns = model.reactions.get_by_id(ids)
            check_rxns = check_model.reactions.get_by_id(ids)
            all_mets = [m.id for m in rxns.metabolites]
            if len(all_mets) != 1:
                # formulas = model.metabolites.get_by_id(all_mets[0]).formula if len(all_mets) == 1 else ""
                exchange_info = {"Reaction_id":ids,
                                "Equation(ID)":rxns.build_reaction_string(),
                                "Formula":"",
                                "Original":str(check_rxns.bounds),
                                "L_b":rxns.lower_bound,
                                "U_b":rxns.upper_bound}
                temp.append(exchange_info)
        # need_fluxes=model.optimize().fluxes[abs(model.optimize().fluxes)>1e-6]
        for rxn in model.reactions:
            all_mets = [m.id for m in rxn.metabolites]
            if len(all_mets) == 1:
                check_rxns = check_model.reactions.get_by_id(rxn.id)
                mets = model.metabolites.get_by_id(all_mets[0])
                exchange_info = {"Reaction_id":rxn.id,
                                "Equation(ID)":rxn.build_reaction_string(use_metabolite_names=True),
                                 "Formula":mets.formula,
                                 "Original":str(check_rxns.bounds),
                                 "L_b":rxn.lower_bound,
                                 "U_b":rxn.upper_bound}
                temp.append(exchange_info)
        return temp

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
                        if not objectiveId:
                            continue
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

    def get_reducing_power_initial_flux(self, model_info, model, model_control_info):
        """"""
        nadh_name, nadh_flux = ['NADH','NADPH','FADH2','FMNH2','Q8H2','MQL8','DMMQL8'], []
        reas = [NADH_NAME,NADPH,FADH2,FMNH2,Q8H2_NAME,MQL8,DMMQL8]
        pro1s = [NAD_NAME,NADP,FAD,FMN,Q8_NAME,MQN8,DMMQ8]
        pro2 = H_NAME
        for i, name in enumerate(nadh_name):
            nadh_id = get_other_nadh_rxn(model, reas[i], pro1s[i], pro2, name)
            if not nadh_id:
                continue
            model.reactions.get_by_id(nadh_id).bounds = (0,1000)
            set_model_objective(model_info, model, nadh_id)   
            initial_nadh = f"initial_{name}".lower()
            model_info['initial_flux'][initial_nadh] = model_info["now_obj_rxn"]
            # model_control_info['initial_flux'][initial_nadh] = round(model_info["now_obj_rxn"]["rxn_flux"],3)
            nadh_flux.extend([f"{name} : {round(model_info['now_obj_rxn']['rxn_flux'],3)}"])
        # model_control_info['initial_flux']["initial_reducing_power"] = ', '.join(nadh_flux)


    def get_energy_initial_flux(self, model_info, model, model_control_info):
        """"""
        xtp, xtp_flux = ['A','C','G','U','I'], []
        for k in xtp:
            atpm_id = get_atpm_rxn(model_info, model, f'{k}TP', f'{k}DP')
            if not atpm_id:
                continue
            model.reactions.get_by_id(atpm_id).bounds = (0,1000)
            set_model_objective(model_info, model, atpm_id)
            model_info['initial_flux'][f"initial_{k}tp"] = model_info["now_obj_rxn"]
            # model_control_info['initial_flux'][f"initial_{k}tp"] = round(model_info['now_obj_rxn']['rxn_flux'],3)
            xtp_flux.extend([f"{k}TP : {round(model_info['now_obj_rxn']['rxn_flux'],3)}"])
        # model_control_info['initial_flux']["initial_atp"] = ', '.join(xtp_flux)


    def check_yield_formula(self, model, metId):
        """"""
        calculate_met = ['C','H','O','N','P','S']
        met = model.metabolites.get_by_id(metId)
        if met.elements and (not set(met.elements.keys()).issubset(set(calculate_met))):
            return 0
        if pd.isna(met.formula) or 'R' in met.formula or 'X' in met.formula or met.formula == 'nan' or not met.formula or met.formula == 'null':
            return 0
        return 1

    def find_yield_generation_subprocess(self, model_info, model, metabolites,conn):
        """Get all carbonaceous species that exceed the maximum theoretical yield"""
        yield_generation, c_num = [], 0
        initial_yield_flux, temp, substance_flux = [], {}, []
        temps = []
        
        pid = multiprocessing.current_process().pid
        # print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),f'pid {pid} model.metabolites num:',len(metabolites))
        
        for met in metabolites:
            if met.id in model_info['need_set_carbon_mets']:
                continue
            if re.search("C[A-Z]",str(met.formula)) or re.search("C[\d]",str(met.formula)): # Carbon followed by capital letters or numbers
                c_num += 1
                if self.check_yield_formula(model, met.id) == 0:
                    # print(met,': formula不明确,不计算')
                    continue
                try:
                    with model:
                        if met.name in ACCOA and met.compartment in C_COMPARTMENT:
                            objectiveId = add_accoa(model, str(met.id))
                        else:
                            objectiveId = add_demand(model_info, model, met.id)
                        model.objective = objectiveId
                        # write_flux_file(model_info, model, objectiveId,self.cfg)
                        max_rate, temp = max_reduced_degree_rate(model, met.id, model_info)
                        model_info['yield_note'].extend(temp)
                        if max_rate == 0 or 'R' in met.formula or 'X' in met.formula: # 含碳物质，最大得率如果是0不用算；formula含有[CHONPS]以外的不算
                            continue
                        
                        model_opt = model.slim_optimize()
                        # print(f'pid {pid} met: ',met.id,' object: ',model_opt,'  max_rate: ',max_rate)
                        if round(model_opt,2) > round(max_rate,2):
                            # print(met.id,round(model_opt,2),round(max_rate,2),'..................')
                            yield_generation.append(met.id)
                            initial_yield_flux.extend([f"{met.id} : {round(model_opt,3)}"])
                            temp[f"{met.id}"] = round(model_opt,3)
                            approach_info = get_net_approach(model, met.id, max_rate)
                            temps.append(approach_info)
                            # pfba_solution = pfba(model)  # 限定通量值v大于1e-6的反应
                            # need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]
                            # substance_flux = get_substance_approach(model, need_fluxes,  met.id)
                            # get_visualiza_html(model, objectiveId, need_fluxes, self.cfg)
                except:
                    pass 
        # print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),f' pid {pid} yield_generation:',yield_generation)
        
        conn.send({
            "yield_generation":yield_generation,
            "c_num":c_num,
            "initial_yield_flux":initial_yield_flux,
            "met_flux":temp,
            "initial_substance_approach" : temps
            # "substance_flux":substance_flux
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
        initial_substance_approach = []
        substance_flux = []
        for parent_connection in parent_connections:
            response = parent_connection.recv()
            yield_generation += response["yield_generation"]
            c_num += response["c_num"]
            initial_yield_flux += response["initial_yield_flux"]
            temp_met_flux.update(response["met_flux"])
            initial_substance_approach += response["initial_substance_approach"]
            # substance_flux += response["substance_flux"]
            
            
        model_info['yield_generation'] = yield_generation
        # model_info['substance_flux'] = substance_flux
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
        model_control_info["check_metabolite_yield"]["description0"] = "Initial information on all substances exceeding the maximum theoretical yield"
        model_control_info["check_metabolite_yield"]["initial_substance_approach"] = initial_substance_approach


    def get_model_initial_flux(self, model_info, model, model_control_info, model_check_info):
        """"""
        with model:
            model_info['initial_flux'] = {}
            model_info['initial_flux']["initial_biomass"] = model_info["initial_rxn"]
            # model_control_info['initial_flux']["initial_biomass"] = round(model_info["initial_rxn"]["biomass_rxn_flux"],3)
            # if "ADD_NADH" in model.reactions:
            #     nadh_id = "ADD_NADH"
            # else:
            #     nadh_id = add_nadh_rxn(model)
            # model.reactions.get_by_id(nadh_id).bounds = (0,1000)
            # set_model_objective(model_info, model, nadh_id)    
            # model_info['initial_flux']["initial_nadh"] = model_info["now_obj_rxn"]
            # model_control_info['initial_flux']["initial_nadh"] = round(model_info["now_obj_rxn"]["rxn_flux"],3)
            # # model.reactions.remove('ADD_NADH')
            self.get_reducing_power_initial_flux(model_info, model, model_control_info)
            self.get_energy_initial_flux(model_info, model, model_control_info)
            # self.find_net_initial_flux(model_info, model, model_control_info)
            model_control_info["check_metabolite_yield"]["score"] = 1
            # model_control_info['initial_flux']["initial_yield"] = ""
            if is_autotrophic(model_info):
                model_control_info["check_metabolite_yield"]["autotrophic"] = "This is an autotrophic strain, no need to calculate the yield"
                model_check_info["check_metabolite_yield"]["autotrophic"] = "This is an autotrophic strain, no need to calculate the yield"
                # model_control_info['initial_flux']["initial_yield"] = "This is an autotrophic strain, no need to calculate the yield"
                return 0
            else:
                get_c_source_supply(model_info, model)
                self.find_yield_initial_flux(model_info, model, model_control_info)

    def initial_control(self, model_info, model, check_model, model_control_info, model_check_info):
        """""" 
        # model_check_info['boundary_information']["score"] = 1
        model_info["control_analysis"].append(model.id)
        # model_control_info['initial_flux'] = {}
        del_model_boundary(model_info, model, check_model, model_check_info)
        model_check_info['boundary_information']['biomass_boundary'] = model_info['biomass_boundary']
        model_control_info['boundary_information']['carbon_source_supply'] = ', '.join(model_info['carbon_source_supply'])
        print('初始biomass值:',model.slim_optimize(),'........................................111.......................')

        set_c_source_supply(model_info, model, 'initial', check_model)
        model_info["biomass_by_c_source"] = round(model.slim_optimize(), 3)
        print('设置碳源后biomass值:',model.slim_optimize(),'.......................................222........................')
    
        # self.get_model_initial_flux(model_info, model, model_control_info, model_check_info)
        model_info['initial_flux'] = {}
        set_special_boundary_rxn(model_info, model, check_model)
        set_bio_boundary(model_info, model)
        set_no2_no3_bounds(model_info, model)
     
        model_control_info['boundary_information']['non_carbon_source_limiting_substrate'] = get_not_carbon_source_limiting_metabolites(model_info, model, check_model)
        # model_check_info['boundary_information']['pre_modify_bounds'] = ','.join(model_info['pre_modify_bounds'])
        model_check_info['boundary_information']['special_boundary_reaction'] = ""
        model_check_info['boundary_information']['min_carbon_source_set'] = []
        model_info["biomass_by_special_boundary"] = round(model.slim_optimize(), 3)
        print('设置特殊边界后biomass值:',model.slim_optimize(),'.......................................333........................')
        # all_met = ['cpd11493_c0','cpd00039_c0','cpd00051_c0','cpd00028_c0','cpd00003_c0','cpd00060_c0','cpd00357_c0','cpd00087_c0','cpd00053_c0','cpd00557_c0','cpd00041_c0','cpd00161_c0','cpd15748_c0','cpd15667_c0','cpd00104_c0','cpd00345_c0','cpd00002_c0','cpd00264_c0','cpd00035_c0','cpd00241_c0','cpd00017_c0','cpd15352_c0','cpd00023_c0','cpd00056_c0','cpd15793_c0','cpd00220_c0','cpd00052_c0','cpd15500_c0','cpd15668_c0','cpd00118_c0','cpd15766_c0','cpd00201_c0','cpd00033_c0','cpd00119_c0','cpd15540_c0','cpd00132_c0','cpd00084_c0','cpd00042_c0','cpd00107_c0','cpd00050_c0','cpd15757_c0','cpd00356_c0','cpd15560_c0','cpd11459_c0','cpd00015_c0','cpd15669_c0','cpd15533_c0','cpd00010_c0','cpd00115_c0','cpd00065_c0','cpd00054_c0','cpd00062_c0','cpd00156_c0','cpd00037_c0','cpd15775_c0','cpd00166_c0','cpd00038_c0','cpd15665_c0','cpd00129_c0','cpd00006_c0','cpd00322_c0']
        # for a_met in all_met:
        #     with model:
        #         mets=model.metabolites.get_by_id(a_met)
        #         demand_rxn = model.add_boundary(mets, type = 'demand') # 添加demand反应
        #         model.objective = demand_rxn.id
        #         print(model.slim_optimize(),demand_rxn.id)
        # exit()
        set_auto_source_supply(model_info, model, check_model, 'initial')
        # pfba_solution = pfba(model)
        # need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]
        # with open('/home/dengxiao/mqc/tmp/bigg/iJN746/flux1.txt', 'w') as flux:  
        #     for r,v in need_fluxes.items():
        #         rxns = model.reactions.get_by_id(r) 
        #         try:
        #             check = rxns.check_mass_balance()  
        #             flux.write(f"{r}\t{round(v,5)}\t{rxns.reaction}\t{rxns.build_reaction_string(use_metabolite_names=True)}\t{rxns.bounds}\t{check}\n") 
        #         except:
        #             flux.write(f"{r}\t{round(v,5)}\t{rxns.reaction}\t{rxns.build_reaction_string(use_metabolite_names=True)}\t{rxns.bounds}\tno_check_mass_balance\n") 
        model_control_info['boundary_information']['carbon_source_boundary'] = ','.join(model_info["carbon_source_boundary"])
        exchange_rxns = self.get_exchange(model, model_info, check_model)
        model_check_info['boundary_information']['exchange_reaction_boundary'] =  exchange_rxns
        model_control_info['boundary_information']['reaction_with_special_boundary'] = ','.join(model_info['special_boundary_reaction'])
        model_control_info['boundary_information']['exchange_reaction_boundary'] =  exchange_rxns
        
        # set_special_boundary_rxn2(model_info, model)
        # min_carbon_source_set = model_info['min_carbon_source_set']
        # if min_carbon_source_set:
        #     if type(min_carbon_source_set[0]) != str:
        #         model_check_info['boundary_information']['min_carbon_source_set'] = min_carbon_source_set
        #     else:
        #         model_check_info['boundary_information']['min_carbon_source_set'] = ', '.join(min_carbon_source_set)
        # else:
        #     model_check_info['boundary_information']['min_carbon_source_set'] = []
       
        print(model.slim_optimize(),'...............................................................')
        rules = Rules(self.cfg) 
        run_rules(rules, model_info, model)
        # rules.get_all_rules(model_info, model)
        # check_rxn_balance(model_info, model)
        # check_C_balance(model_info, model)
        # net_penalty_points(model_info)
        
      
        print(model.slim_optimize(),'...............................................................')
        # exit()
        # write_flux_file2(model_info, model, 'BIOMASS__3',self.cfg)
        # model.reactions.get_by_id('EX_glc__D_e').bounds=(0,0) # 0
        # print(model.slim_optimize(),',,,')   # 0.25
        # model.reactions.get_by_id('DNA').bounds=(0,0)
        # print(model.slim_optimize(),'...')
        # model.reactions.get_by_id('DNA').bounds=(0,1000)
        # print(model.slim_optimize(),'///')
        # exit()
        
        
        
        