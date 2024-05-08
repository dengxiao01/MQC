# -*- coding: utf-8 -*-
import re ,time

from mqc.config import Config
from mqc.utils import *
import multiprocessing
from multiprocessing import Pipe, Process

class Yields2():
    """"""

    def __init__(self,cfg:Config) -> None:
        self.cfg = cfg



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


    def find_yield_generation(self, model_info, model, model_control_info):
        """Get all carbonaceous species that exceed the maximum theoretical yield"""
        yield_generation, c_num = [], 0
        initial_yield_flux = []
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),' model.metabolites num:',len(model.metabolites))
        
        for met in model.metabolites:
            if met.id in model_info['need_set_carbon_mets']:
                continue
            if re.search("C[A-Z]",str(met.formula)) or re.search("C[\d]",str(met.formula)): # Carbon followed by capital letters or numbers
                c_num += 1
                if self.check_yield_formula(model, met.id) == 0:
                    print(met,': formula不明确,不计算')
                    continue
                
                with model:
                    objectiveId = add_demand(model_info, model, met.id)
                    model.objective = objectiveId
                    # write_flux_file(model_info, model, objectiveId,self.cfg)
                    
                    max_rate = max_reduced_degree_rate(model, met.id)
                    
                    if max_rate == 0 or 'R' in met.formula or 'X' in met.formula: # 含碳物质，最大得率如果是0不用算；formula含有[CHONPS]以外的不算
                        continue
                    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),' model_opt start.',met.id)
                    model_opt = model.slim_optimize()
                    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),' model_opt end.',met.id)
                    print('met: ',met.id,' object: ',model_opt,'  max_rate: ',max_rate)
                    if round(model_opt,2) > round(max_rate,2):
                        yield_generation.append(met.id)
                        initial_yield_flux.extend([f"{met.id} : {round(model_opt,2)}"])
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),' yield_generation:',yield_generation)
        model_info['yield_generation'] = yield_generation
        if len(yield_generation) != 0:
            model_control_info["check_metabolite_yield"]["score"] = 0
        model_control_info["check_metabolite_yield"]["total_carbon_content:"] = c_num
        model_control_info["check_metabolite_yield"]["yield_exceeds_the_theoretical_maximum_number:"] = len(yield_generation)
        model_control_info["check_metabolite_yield"]["yield_net_matter:"] = ', '.join(yield_generation)
        model_control_info['initial_flux']["initial_yield"] = ', '.join(initial_yield_flux)

    def find_yield_generation_subprocess(self, model_info, model, metabolites,conn):
        """Get all carbonaceous species that exceed the maximum theoretical yield"""
        yield_generation, c_num = [], 0
        initial_yield_flux, temp, substance_flux = [], {}, []
        
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
                        if not objectiveId:
                            continue
                        model.objective = objectiveId
                        # write_flux_file(model_info, model, objectiveId,self.cfg)
                        max_rate = max_reduced_degree_rate(model, met.id)
                        if max_rate == 0: # 含碳物质，最大得率如果是0不用算；
                            continue
                        
                        model_opt = model.slim_optimize()
                        # print(f'pid {pid} met: ',met.id,' object: ',model_opt,'  max_rate: ',max_rate)
                        if round(model_opt,2) > round(max_rate,2):
                            yield_generation.append(met.id)
                            initial_yield_flux.extend([f"{met.id} : {round(model_opt,3)}"])
                            temp[f"{met.id}"] = round(model_opt,3)
                            pfba_solution = pfba(model)  # 限定通量值v大于1e-6的反应
                            need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]
                            substance_flux = get_substance_approach(model, need_fluxes,  met.id)
                            get_visualiza_html(model, objectiveId, need_fluxes, self.cfg)
                except:
                    pass 
        # print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),f' pid {pid} yield_generation:',yield_generation)
        
        conn.send({
            "yield_generation":yield_generation,
            "c_num":c_num,
            "initial_yield_flux":initial_yield_flux,
            "met_flux":temp,
            "substance_flux":substance_flux
        })
        conn.close()
        

    def find_yield_generation_multiprocessing(self,model_info, model, model_control_info):
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
        substance_flux = []
        for parent_connection in parent_connections:
            response = parent_connection.recv()
            yield_generation += response["yield_generation"]
            c_num += response["c_num"]
            initial_yield_flux += response["initial_yield_flux"]
            temp_met_flux.update(response["met_flux"])
            substance_flux += response["substance_flux"]
            
            
        model_info['yield_generation'] = yield_generation
        model_info['substance_flux'] = substance_flux
        if len(yield_generation) != 0:
            model_control_info["check_metabolite_yield"]["score"] = 0
        model_info["total_carbon_content"] = c_num
        model_info["yield_exceeds_the_theoretical_maximum_number"] = len(yield_generation)
        model_info["yield_net_matter"] = ', '.join(yield_generation)
        # if len(initial_yield_flux) != 0:
        #     model_control_info['initial_flux']["initial_yield"] = ', '.join(initial_yield_flux)
        # else:
        #     model_control_info['initial_flux']["initial_yield"] = 'no yield net substance'
        model_control_info['initial_flux']["initial_yield"] = temp_met_flux


    def close_max_value_yield(self, model_info, model_control_info, temp, model, object_rxn_id, yieldId, max_rate, check_model, model_check_info):
        """"""
        grate_delete_rxn, count, fluxs, objects, count2, rely_rxn = [], 0, 0, '', 0, []
        if round(model.slim_optimize(),2) <= round(max_rate,2):
            # print(model,yieldId,'________得率_yes')
            pass
        while round(model.slim_optimize(),2) > round(max_rate,2): 
            count += 1
            print(count,' -- ',object_rxn_id)
            if objects == object_rxn_id and fluxs == round(model.slim_optimize(),3):
                count2 += 1
            else:
                fluxs = round(model.slim_optimize(),3)
                objects = object_rxn_id
                count2 = 0
            if count2 > 3:
                print('error: 得率无限循环')
                model_control_info["check_metabolite_yield"]["score"] = 0
                temp_related_reaction = get_other_net_infinite_info(model_info, model, need_fluxes)
                temp['related_reaction'].extend(temp_related_reaction)
                model_control_info["check_metabolite_yield"]["model_revision"] = temp["model_revision"]
                model_control_info["check_metabolite_yield"]["related_reaction"] = temp["related_reaction"]
                model_check_info["check_metabolite_yield"] = temp
                del model_check_info["check_metabolite_yield"]["model_revision"]
                del model_check_info["check_metabolite_yield"]["related_reaction"]
                model_check_info["check_metabolite_yield"]["error"] = f"error: Yield infinite loop, biomass growth must depend on {'、'.join(rely_rxn)}, but this must be closed"
                raise RuntimeError(f"error: Yield infinite loop")
            # write_flux_file(model_info, model, object_rxn_id,self.cfg)
            # write_flux_file3(model_info, model, object_rxn_id,self.cfg)
            fba_rxn_dic = {}
            pfba_solution = pfba(model)  # 限定通量值v大于1e-6的反应
            need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]
            # need_fluxes=model.optimize().fluxes[abs(model.optimize().fluxes)>1e-6]
            for ids,v in need_fluxes.items():
                for rxn in model_info['reactions']:
                    if ids == rxn['id']:
                        if ids == object_rxn_id: 
                            continue
                        fba_rxn_dic[ids] = rxn['net_penalty_points']
            with model:
                for ids,v in need_fluxes.items():
                    model.objective = model_info["initial_rxn"]["biomass_rxn_id"]
                    if ids == object_rxn_id:
                        continue
                    model.reactions.get_by_id(ids).bounds = (0,0)
                    if model.slim_optimize() < 1e-6:
                        rely_rxn.append(ids)
                    model.reactions.get_by_id(ids).bounds = check_model.reactions.get_by_id(ids).bounds

            for ids,penalty_points in fba_rxn_dic.items():
                if penalty_points == max(fba_rxn_dic.values()): #关掉每次结果中的所有最大值
                    if fba_rxn_dic[ids] != 0 and ids not in rely_rxn:  
                        model.reactions.get_by_id(ids).bounds = (0,0)
                        grate_delete_rxn.append(ids)
                        print(ids,' ',model.reactions.get_by_id(ids).build_reaction_string(use_metabolite_names=True),fba_rxn_dic[ids])
            max_rate = max_reduced_degree_rate(model, yieldId)
            print('close :  objective_value:',model.slim_optimize(),' max_rate: ',max_rate)
        return grate_delete_rxn

    def old_get_associated_substances(self, modelInfo, model, temp, infinite_rxn, all_yield, check_model, model_control_info):
        """"""
        temp_yield = []
        
        for yieldId in all_yield:
            with model:
                met = model.metabolites.get_by_id(yieldId) 
                if met.name in ACCOA and met.compartment in C_COMPARTMENT:
                    object_rxn_id = add_accoa(model, str(met.id))
                else:
                    object_rxn_id = add_demand(modelInfo, model, yieldId)
                if not object_rxn_id:
                    continue
                model.objective = object_rxn_id
                temp_max_rate = max_reduced_degree_rate(model, yieldId)
                if round(model.slim_optimize(),2) <= round(temp_max_rate,2):
                    temp_yield.append(yieldId)
        # nets_dict[','.join(infinite_rxn)] = temp_net  # 每个净物质处理后都会影响部分其他净物质，nets_dict字典记录了相关联的物质{'for_c': ['for_c', 'hco3_c', 'mmcoa__S_c', 'succoa_c', 'cbp_c', 'urea_c', 'allphn_c', 'agps23_c'], 'malcoa_c': ['malcoa_c', 'prpncoa_c', '3hpcoa_c', 'ppcoa_c', '3opcoa_c'], 'ag160_c': ['ag160_c', '2agpe160_c']}
        temp["dependence_on_carbonaceous_substances"][f"{','.join(infinite_rxn)}"] = ','.join(temp_yield)
        return temp_yield
    
    def get_associated_substances(self, modelInfo, model, temp, all_yield, check_model):
        """"""
        temp_yield = []
        for rxnId in modelInfo['modified_rxn_id'].keys():
            model.reactions.get_by_id(rxnId).bounds = check_model.reactions.get_by_id(rxnId).bounds
        for rxnId,bounds in modelInfo['modified_rxn_id'].items():
            model.reactions.get_by_id(rxnId).bounds = bounds
            for yieldId in all_yield:
                with model:
                    met = model.metabolites.get_by_id(yieldId) 
                    if met.name in ACCOA and met.compartment in C_COMPARTMENT:
                        object_rxn_id = add_accoa(model, str(met.id))
                    else:
                        object_rxn_id = add_demand(modelInfo, model, yieldId)
                    if not object_rxn_id:
                        continue
                    model.objective = object_rxn_id
                    temp_max_rate = max_reduced_degree_rate(model, yieldId)
                    if round(model.slim_optimize(),2) <= round(temp_max_rate,2):
                        temp_yield.append(yieldId)
            # nets_dict[','.join(infinite_rxn)] = temp_net  # 每个净物质处理后都会影响部分其他净物质，nets_dict字典记录了相关联的物质{'for_c': ['for_c', 'hco3_c', 'mmcoa__S_c', 'succoa_c', 'cbp_c', 'urea_c', 'allphn_c', 'agps23_c'], 'malcoa_c': ['malcoa_c', 'prpncoa_c', '3hpcoa_c', 'ppcoa_c', '3opcoa_c'], 'ag160_c': ['ag160_c', '2agpe160_c']}
            if temp_yield:
                temp["dependence_on_carbonaceous_substances"][rxnId] = ','.join(temp_yield)
                # for k in temp_yield:
                #     all_yield.remove(k)
                temp_yield = []
 


    



    def yield_control(self, model_info, model, check_model, model_control_info, model_check_info):
        """"""
        all_yield, temp, all_need_fluxes = [], {}, {}
        set_c_source_supply(model_info, model, 'yields', check_model)
        set_auto_source_supply(model_info, model, check_model, 'yields')
        model_control_info["check_metabolite_yield"]["score"] = 1
        # model_control_info['initial_flux']["initial_yield"] = ""
        if is_autotrophic(model_info):
            model_control_info["check_metabolite_yield"]["autotrophic"] = "This is an autotrophic strain, no need to calculate the yield"
            model_control_info["model_check_info"]["autotrophic"] = "This is an autotrophic strain, no need to calculate the yield"
            model_control_info['initial_flux']["initial_yield"] = "This is an autotrophic strain, no need to calculate the yield"
            return 0
        else:     
            # self.find_yield_generation(model_info, model, model_control_info)
            self.find_yield_generation_multiprocessing(model_info, model, model_control_info)
            temp["total_carbon_content"] = model_info["total_carbon_content"]
            temp["yield_exceeds_the_theoretical_maximum_number"] = model_info["yield_exceeds_the_theoretical_maximum_number"]
            temp["yield_net_matter"] = model_info["yield_net_matter"]
            all_yield.extend(model_info['yield_generation'])
            temp["dependence_on_carbonaceous_substances"] = {}
            # model_control_info['check_metabolite_yield']['model_revision'] = []
            temp['model_revision'], temp['related_reaction'] = [], []
            for yieldId in model_info['yield_generation']:
                print('yieldId:',yieldId)          
                with model:
                    met = model.metabolites.get_by_id(yieldId)
                    if met.name in ACCOA and met.compartment in C_COMPARTMENT:
                        object_rxn_id = add_accoa(model, str(met.id))
                    else:
                        object_rxn_id = add_demand(model_info, model, yieldId)
                    if not object_rxn_id:
                        continue
                    model.objective = object_rxn_id 
                    
                    max_rate = max_reduced_degree_rate(model, yieldId)
                    approach_info = get_net_approach(model, yieldId, max_rate)
                    temp['substance_approach'].append(approach_info)
                    print('objective_value: ',model.slim_optimize(),' max_rate: ',max_rate)
                    grate_delete_rxn = self.close_max_value_yield(model_info, model_control_info, temp, model, object_rxn_id, yieldId, max_rate, check_model, model_check_info)
                    if len(grate_delete_rxn) == 0 :
                        continue
                    need_flux, infinite_rxn = add_back_reactions_net(model_info, model, check_model, grate_delete_rxn, 'yields', max_rate, self.cfg)
                    all_need_fluxes.update(need_flux)
                    temp_model_revision = get_net_infinite_info(model_info, model, 'yields')
                    temp['model_revision'].extend(temp_model_revision)
                    temp_related_reaction = get_other_net_infinite_info(model_info, model, all_need_fluxes)
                    temp['related_reaction'].update({','.join(infinite_rxn) : temp_related_reaction})
                    # temp_net = self.get_associated_substances(model_info, model, temp, infinite_rxn, all_yield, check_model, model_control_info)
                    # for k in temp_net:
                    #     all_yield.remove(k)
                boundary_restoration(model_info, model, 'yields')
            # model_info['temp_final_flux']["final_yield"] = get_final_yield_fluxes(model_info, model)
            # temp['related_reaction'] = get_other_net_infinite_info(model_info, model, all_need_fluxes)
            self.get_associated_substances(model_info, model, temp, all_yield, check_model)
            if model_control_info["check_metabolite_yield"]["score"] == 0:
                # temp["score"] = 0
                model_control_info["check_metabolite_yield"]["substance_approach"] = temp["substance_approach"]
                model_control_info["check_metabolite_yield"]["substance_related_reaction"] = model_info['substance_flux']
                model_control_info["check_metabolite_yield"]["model_revision"] = temp["model_revision"]
                model_control_info["check_metabolite_yield"]["related_reaction"] = temp["related_reaction"]
                model_check_info["check_metabolite_yield"] = temp
                del model_check_info["check_metabolite_yield"]["substance_approach"] 
                del model_check_info["check_metabolite_yield"]["model_revision"]
                del model_check_info["check_metabolite_yield"]["related_reaction"]
                # model_check_info["check_metabolite_yield"] = temp
            else:
                model_control_info["check_metabolite_yield"]["summarize"] = f"There are a total of {temp['total_carbon_content']} carbonaceous substances, none of which exceed the maximum theoretical yield, so pass"
                model_check_info["check_metabolite_yield"]["summarize"] = f"There are a total of {temp['total_carbon_content']} carbonaceous substances, none of which exceed the maximum theoretical yield, so pass"
