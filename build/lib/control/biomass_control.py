# -*- coding: utf-8 -*-

import cobra
import time 
import numpy as np
import cplex 
from collections import OrderedDict
from mqc.config import Config
from mqc.utils import *
from mqc.control.preprocessing_control import Preprocess
from mqc.control.model_control import ModelPreprocess


class Biomasses():


    def __init__(self,cfg:Config):
        """
        Get model information and check_model object.
        """
        self.cfg = cfg

    def recover_bio_objective(self, model_info, model):
        """"""
        model.objective = model_info["initial_rxn"]["biomass_rxn_id"]
  


    def test(self, model, check_model):
        """"""
        # model.reactions.get_by_id('EX_h_e').bounds = (0,1000)
        # model.reactions.get_by_id('FDR').bounds =(0,0)
        # model.reactions.get_by_id('G3POACTF').bounds = (0,0)
        # model.reactions.get_by_id('MADH').bounds = (0,0)
   

    def get_biomass_info(self, model_info, model_control_info):
        """"""
        model_info['biomass_original'] = model_info['initial_rxn']['biomass_rxn_flux']
        temp = [f"{model_info['initial_rxn']['biomass_rxn_id']} : {model_info['initial_rxn']['biomass_rxn_exp']}\n"]
        temp.extend([f"Initial growth rate : {model_info['biomass_original']}\n"])
        temp.extend([model_info["initial_rxn"]["check_biomass_composition"]])
        write_biomass(model_control_info, temp, "Biomass_information")
        

    def add_rely_biomass_rxn(self, model_info, model, check_model, model_control_info, initial_rxn_id):
        """"""
        temp, rely_rxn = [f"biomass-dependent responses:"], []
        with model:
            for rxn in model_info["reactions"]:
                if 'nadhs_modify' in rxn.keys() or 'atps_modify' in rxn.keys() or 'nets_modify' in rxn.keys() or 'yields_modify' in rxn.keys():
                    model.reactions.get_by_id(rxn['id']).bounds = check_model.reactions.get_by_id(rxn['id']).bounds
            for rxn in model_info["reactions"]:
                if 'nadhs_modify' in rxn.keys() or 'atps_modify' in rxn.keys() or 'nets_modify' in rxn.keys() or 'yields_modify' in rxn.keys():
                    model.objective = initial_rxn_id
                    model.reactions.get_by_id(rxn['id']).bounds = (0,0)
                    print(rxn['id'],model.reactions.get_by_id(rxn['id']).bounds,model.slim_optimize(),'//////////////////////////')
                    if model.slim_optimize() <= 1e-6:
                        rely_rxn.append(rxn['id'])
                        rxn["rely_biomass"] = "true"
                    model.reactions.get_by_id(rxn['id']).bounds = check_model.reactions.get_by_id(rxn['id']).bounds
                    # for ids in model_info['modified_rxn_id'].keys():
                    #     check_model.reactions.get_by_id(ids).bounds = (0,0)
                    #     if check_model.slim_optimize() <= 1e-6:
                    #         rely_rxn.append(rxn['id'])
                    #         rely_rxn.append(ids)
                    #         rxn["rely_biomass"] = "true"
            # if 'nadhs_modify' not in rxn.keys() and 'atps_modify' not in rxn.keys() and 'nets_modify' not in rxn.keys() and 'yields_modify' not in rxn.keys() and 'carbon_source' not in rxn.keys():
            #     if rxn['id'] != initial_rxn_id:
            #         model.reactions.get_by_id(rxn['id']).bounds = check_model.reactions.get_by_id(rxn['id']).bounds
        print('rely_biomass_rxn:............................')
        rely_rxn = list(set(rely_rxn))
        for rxnId in rely_rxn:
            model.reactions.get_by_id(rxnId).bounds = check_model.reactions.get_by_id(rxnId).bounds
            print(rxnId,model.reactions.get_by_id(rxnId).bounds)
        # model.reactions.get_by_id('GLCD').bounds = check_model.reactions.get_by_id('GLCD').bounds
        if rely_rxn:
            temp.extend(rely_rxn)
            model_control_info["check_biomass_production"]["rely_biomass_rxn"] = ' '.join(temp)
            # write_biomass(model_control_info, temp, "preprocessing")


    def get_general_library(self):
        """"""
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':读取总库')
        # ks=[]
        # controler = Preprocess(self.cfg.general_library_path,self.cfg)
        # with controler.model:
        #     model_pre_process = ModelPreprocess(self.cfg)
        #     model_pre_process.preprocess_met(controler)
        #     for rxn in controler.model.reactions:
        #         if 'C' in rxn.check_mass_balance():
        #             ks.append(rxn.id)
        # print(ks,'...................')
        model = cobra.io.read_sbml_model(self.cfg.general_library_path)
        # for k in ks:
        #     model.reactions.get_by_id(k).remove_from_model()
        # cobra.io.write_sbml_model(model, "/home/dengxiao/mqc/mqc/summary/test.xml")
        # return cobra.io.load_json_model(self.cfg.general_library_path)
        return model


    def add_sink_rxn(self, big_model, model_info, check_model):
        """"""
        sink_rxn = []
        initial_rxn_id = model_info["initial_rxn"]["biomass_rxn_id"]
        initial_rxn = big_model.reactions.get_by_id(initial_rxn_id)
        for met in initial_rxn.products:
            if met.name not in (ADP_NAME+H_NAME+PI_NAME+H2O_NAME) and 'biomass' not in met.name.lower() and f'SK_{met.id}' not in big_model.reactions:
                big_model.add_boundary(metabolite=met, type="sink")
                # check_model.add_boundary(metabolite=met, type="sink")
                sink_rxn.append(f'SK_{met.id}')
                # big_model.reactions.get_by_id(f'SK_{met.id}').bounds = (0,1000)
                # check_model.reactions.get_by_id(f'SK_{met.id}').bounds = (0,1000)
                print('添加sink反应的代谢物:',met,'............................')
        model_info["sink_rxn"] = sink_rxn
        # return sink_rxn


    def get_big_general_library(self, model_info, model, general_library):
        """"""
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':获取大模型')
        initial_rxn_id = model_info["initial_rxn"]["biomass_rxn_id"]
        model_file = model_info["model_file"]
        if self.cfg.big_general_library:
            return 0
        with model:
            for rxns in model.reactions:  # 给系数扩大1000倍,添加的demand反应不改
                temp_dict = {}
                for met in rxns.metabolites:
                    temp_dict[met] = rxns.get_coefficient(met) * (1000)
                rxns.add_metabolites(temp_dict) # 相同时会覆盖掉原代谢物。乘以1000是为了pfba计算
            # for rxns in model.reactions:
            #     if rxns.id not in general_library.reactions and len(rxns.products) != 0:
            #         general_library.add_reaction(rxns) 
            add_rxn = []
            for rxns in general_library.reactions:
                if rxns.id not in model.reactions and rxns.id != initial_rxn_id and len(rxns.metabolites) != 1 and 'C' not in rxns.check_mass_balance():
                    model.add_reaction(rxns)  
                    add_rxn.append(rxns.id)
            # self.convert_nan_to_null(model)
            delect_rxn = []
            for rxn in model.reactions: # 添加的反应有些代谢物的formula会变，导致不平衡，所有在添加后还需要删一次
                if rxn.id in add_rxn and 'C' in rxn.check_mass_balance():
                    delect_rxn.append(rxn.id)
            for rxnId in delect_rxn:
                model.reactions.get_by_id(rxnId).remove_from_model()
            convert_nan_to_null(model)
            if model_file.endswith(".json"):
                cobra.io.save_json_model(model, self.cfg.big_general_library_json)
                self.cfg.big_general_library = self.cfg.big_general_library_json
            elif model_file.endswith(".yaml"):
                cobra.io.save_yaml_model(model,self.cfg.big_general_library_yaml)
                self.cfg.big_general_library = self.cfg.big_general_library_yaml
            elif model_file.endswith(".mat"):
                cobra.io.save_matlab_model(model,self.cfg.big_general_library_mat)
                self.cfg.big_general_library = self.cfg.big_general_library_mat
            else:
                cobra.io.write_sbml_model(model, self.cfg.big_general_library_xml)
                self.cfg.big_general_library = self.cfg.big_general_library_xml
            # try:
            #     cobra.io.write_sbml_model(model, self.cfg.big_general_library_xml)
            #     # file_path = self.cfg.big_general_library_xml
            #     self.cfg.big_general_library = self.cfg.big_general_library_xml
            # except:
            #     cobra.io.save_json_model(model, self.cfg.big_general_library_json)
            #     # file_path = self.cfg.big_general_library_json
            #     self.cfg.big_general_library = self.cfg.big_general_library_json
            # print('大模型路径:',self.cfg.big_general_library)
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':获取大模型完毕')
        return model

    def get_big_model(self, model_info, model, general_library):
        """"""
        # if not os.path.exists(self.cfg.big_general_library):
        big_model = self.get_big_general_library(model_info, model, general_library)
        model_file = model_info["model_file"]
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':读大模型')
        if model_file.endswith(".json"):
            big_model = cobra.io.load_json_model(self.cfg.big_general_library)
        elif model_file.endswith(".yaml"):
            big_model = cobra.io.load_yaml_model(self.cfg.big_general_library)
        elif model_file.endswith(".mat"):
            big_model = cobra.io.load_matlab_model(self.cfg.big_general_library)
        else:
            big_model = cobra.io.read_sbml_model(self.cfg.big_general_library)
        # try:     
        #     big_model = cobra.io.read_sbml_model(self.cfg.big_general_library)
        # except:
        #     big_model = cobra.io.load_json_model(self.cfg.big_general_library)
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':读大模型完毕')
        return big_model



    def bio_precursor_gapfilling(self, model_info, big_model, no_synthesized_mets, model_control_info):
        """"""
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':看大模型中有多少个物质不能合成')
        to_gap, no_gap, model_info['no_gap_metid'] = [], [], []
        # test = ["cpd15665_c0","cpd00166_c0","cpd15668_c0","cpd15669_c0","cpd15667_c0"]
        temp = [f"After quality control, the initial rate of biomass is 0, the metabolite {','.join(no_synthesized_mets)} of biomass components cannot be synthesized"]
        with big_model:
            for metId in no_synthesized_mets:
                objectiveId = add_demand(model_info, big_model, metId)
                big_model.objective = objectiveId
                ini_flux = big_model.slim_optimize()      
                if pd.isna(ini_flux) or ini_flux <= 1e-6:
                    no_gap.append(metId)
                    continue
                # if metId in test:
                #     write_flux_file(model_info, big_model, objectiveId, metId)
                to_gap.append(metId)
        # if no_gap:
        #     temp.extend([f"In the original model, the precursor substance {','.join(no_synthesized_mets)} of biomass cannot be synthesized, but in the larger model,{len(no_gap)} metabolites cannot be synthesized by gapfilling: {','.join(no_gap)}"])
        # else:
        #     temp.extend([f"In the original model, the precursor substance {','.join(no_synthesized_mets)} of biomass cannot be synthesized, but in the larger model, all substances can be synthesized by gapfilling"])
        # temp.extend([f"Gap is applied to the substances that can be synthesized: {to_gap}"])
        write_biomass(model_control_info, temp, "Gapfilling_of_biomass_components")
        model_info['no_gap_metid'] = no_gap
        return to_gap

    def get_exchange_rxn(self, model_info, big_model, no_synthesized_mets, model_control_info):
        """"""
        # to_gap_metid = self.bio_precursor_gapfilling(model_info, big_model, no_synthesized_mets, model_control_info)
        # sink_rxn = self.add_sink_rxn(big_model, model_info)
        to_gap_metid = self.bio_precursor_gapfilling(model_info, big_model, no_synthesized_mets, model_control_info)
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':获取能合成物质的交换反应')
        if to_gap_metid:
            reaction = cobra.Reaction("EX_gap")
            for met in to_gap_metid:
                reaction.add_metabolites({big_model.metabolites.get_by_id(met): -1})
            reaction.bounds = (0, 1)
            big_model.add_reaction(reaction)
            return reaction.id 
        else:
            return ''
        

    def get_new_way_rxn(self, model_info, big_model, no_synthesized_mets, model_control_info):
        """"""
        exchange_rxnid = self.get_exchange_rxn(model_info, big_model, no_synthesized_mets, model_control_info)
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':初始biomass为0时以该交换反应为目标获取大模型的途径')
        new_way_rxn = []
        if exchange_rxnid:
            big_model.objective = exchange_rxnid
            big_model.reactions.get_by_id(exchange_rxnid).bounds = (0,1)  # 限制demand反应的通量
            if big_model.slim_optimize() > 0:
                pfba_solution = pfba(big_model) 
                need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]
                # need_fluxes=big_model.optimize().fluxes[abs(big_model.optimize().fluxes)>1e-6]
                new_way_rxn = list(need_fluxes.keys())
            return new_way_rxn
        else:
            return []


    def gapfilling(self, model, model_info, demand_id, general_library):
        """"""
        big_model = self.get_big_model(model_info, model, general_library)
        print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}:在大模型中{demand_id}的途径")
        new_way_rxn = []
        with big_model:
            rxns = model.reactions.get_by_id(demand_id)
            big_model.add_reaction(rxns) 
            big_model.objective = demand_id
            big_model.reactions.get_by_id(demand_id).bounds = (0,1)  # 限制demand反应的通量
            if model.slim_optimize() > 0:
                pfba_solution = pfba(big_model)  # 将模型和总库的反应融合到一起后，再次计算当前目标下的通量结果文件；随后把生成目标的途径添加到模型中
                need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]
                # need_fluxes=model.optimize().fluxes[abs(model.optimize().fluxes)>1e-6]
                new_way_rxn = list(need_fluxes.keys())
        print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}:在大模型中{demand_id}的途径计算完毕,途径为{new_way_rxn}")
        return new_way_rxn

    def biomass_gapfilling(self, model, model_info, initial_rxn_id, general_library):
        """"""
        big_model = self.get_big_model(model_info, model, general_library)
        print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}:在大模型中{initial_rxn_id}的途径")
        new_way_rxn = []
        with big_model:
            print(initial_rxn_id,'....................')
            big_model.objective = initial_rxn_id
            big_model.reactions.get_by_id(initial_rxn_id).bounds = (0,0.1)  # 限制demand反应的通量
            # print(big_model.reactions.get_by_id('ALAB').bounds,'......................')
            print(big_model.slim_optimize(),'..............')
            if model.slim_optimize() > 0:
                pfba_solution = pfba(big_model)  # 将模型和总库的反应融合到一起后，再次计算当前目标下的通量结果文件；随后把生成目标的途径添加到模型中
                need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]
                print('........................')
                # need_fluxes=model.optimize().fluxes[abs(model.optimize().fluxes)>1e-6]
                new_way_rxn = list(need_fluxes.keys())
        print(f"{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}:在大模型中{initial_rxn_id}的途径计算完毕")
        return new_way_rxn

    def add_met_name(self, model, need_add_rxn, general_library, obj):
        """"""
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':添加gap反应的物质名称、formula')
        seed_met, meta_met, virtual_met, kegg_met, other_met = [], [], [], [], []
        pattern = r'^C\d{5}$' 
        for rxnId in need_add_rxn:
            # if rxnId not in model.reactions and rxnId in general_library.reactions:
            rxns = general_library.reactions.get_by_id(rxnId)
            for met in rxns.metabolites:
                if met.id.startswith("cpd"):  # seed
                    seed_met.append(met.id)
                elif "@" in met.id:  # meta
                    meta_met.append(met.id)
                elif "[" in met.id and "]" in met.id:  # virtual
                    virtual_met.append(met.id)
                elif re.match(pattern, met.id):  # kegg
                    kegg_met.append(met.id)
                else:
                    other_met.append(met.id)
        if len(seed_met) != 0:
            obj.add_seed_formula(seed_met)
        if len(meta_met) != 0:
            obj.add_meta_formula(meta_met)
        # if len(virtual_met) != 0:
        #     pass
        if len(kegg_met) != 0:
            obj.add_kegg_formula(kegg_met)
        if len(other_met) != 0:
            obj.add_other_formula(other_met)
      


    def add_gapfilling_rxn(self, model, model_info, check_model, general_library, new_way_rxn, obj):
        """"""
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':开始添加gapfilling反应:')
        need_add_rxn, model_info["gaped_rxn"] = [], []
        if not new_way_rxn:
            gap_info = ["Can not able to gapfill !!!"]
        else:    
            t1 = time.time()
            gap_info = ["The gapfilling process added the reaction:"]
            for rxnId in new_way_rxn:
                if rxnId not in model.reactions and rxnId in general_library.reactions:
                    rxns = general_library.reactions.get_by_id(rxnId)
                    model.add_reaction(rxns)
                    need_add_rxn.append(rxnId)
                    # rxn_exp = [f"{rxns}"]
                    # if model_info['model_identifier'] == 'modelseed' or model_info['model_identifier'] == 'metaNetx':
                    #     rxn_exp = [f"{rxns.id}: {rxns.build_reaction_string(use_metabolite_names=True)}"]
                    # gap_info.extend(rxn_exp)
            print('need_add_rxn: ',need_add_rxn)
            t2 = time.time()
            print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':正在gapfilling添加反应', 'time:',t2-t1)
            self.add_met_name(model, need_add_rxn, general_library, obj)
            for rxnId in need_add_rxn:
                # if rxnId not in model.reactions and rxnId in general_library.reactions:
                rxns = model.reactions.get_by_id(rxnId)
                rxn_exp = [f"{rxns.id}: {rxns.build_reaction_string(use_metabolite_names=True)}"]
                gap_info.extend(rxn_exp)
            if len(need_add_rxn) == 0:
                gap_info.extend(["no gap response,since it already exists!"])
        model_info["gaped_rxn"].extend(need_add_rxn)
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':添加gapfilling反应完毕')
        return need_add_rxn, gap_info


    def add_gap_rxn(self, model, need_add_rxn, general_library, check_model):
        """"""
        for rxnId in need_add_rxn:
            rxns = general_library.reactions.get_by_id(rxnId)
            if rxnId not in model.reactions:
                model.add_reaction(rxns)
                check_model.add_reaction(rxns)
      

    def check_bio_component_is_zero(self, model_info, model, model_control_info, initial_rxn_id):
        """"""
        t1 = time.time()
        no_synthesized_mets, synthesized_mets, num = [], [], 0
        ATP_NADH = [ATP_NAME+CTP_NAME+GTP_NAME+UTP_NAME+ITP_NAME+NADH_NAME+NADPH+FADH2+Q8H2_NAME+FMNH2+MQL8+DMMQL8]
        for met in model.reactions.get_by_id(initial_rxn_id).reactants:
            num += 1
            if met.name in ATP_NADH:
                pass 
            with model:
                if met.name in ACCOA and met.compartment in C_COMPARTMENT:
                    objectiveId = add_accoa(model, str(met.id))
                else:
                    objectiveId = add_demand(model_info, model, str(met.id))
                model.objective = objectiveId
                ini_flux = model.slim_optimize()
                synthesized_mets.extend([f"{str(met.id)} : {round(ini_flux,3)}"])
                # ini_flux = model.optimize().fluxes[objectiveId]
                if pd.isna(ini_flux) or ini_flux <= 1e-6:
                    no_synthesized_mets.append(str(met.id))
                # print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':检查生物量组成物质哪些为0', objectiveId,ini_flux)
        temp = [f"The biomass equation consists of {num} total components"]
        if no_synthesized_mets:
            temp.extend([f"among which the synthesis rate of component {','.join(no_synthesized_mets)} was 0 "])
        else:
            temp.extend([f"none of these precursor species have a rate of 0, namely: {','.join(synthesized_mets)}"])
        # model_control_info["check_biomass_production"]["Check_synthesizability_of_biomass_components"] = ', '.join(temp)
        write_biomass(model_control_info, temp, "Check_synthesizability_of_biomass_components")
        t2 = time.time()
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':检查生物量组成物质哪些为0', t2-t1,'s')
        return no_synthesized_mets


    def gap_initial_zero_bio(self, model_info, model, model_control_info, general_library, no_synthesized_mets, big_model, obj, check_model):
        """"""
        temp = []
        new_way_rxn = self.get_new_way_rxn(model_info, big_model, no_synthesized_mets, model_control_info)
        need_add_rxn, gap_info = self.add_gapfilling_rxn(model, model_info, check_model, general_library, new_way_rxn, obj)
        model_info["bio_modified"].extend(need_add_rxn)
        temp.extend(gap_info)    
        # self.add_sink_rxn(model, model_info)
        # if need_add_rxn: self.add_gap_rxn(model, need_add_rxn, general_library)
        write_biomass(model_control_info, temp, "Gapfilling_of_biomass_components")
        return no_synthesized_mets


    def get_no_synthesized_mets_flux(self, model, model_info, no_synthesized_mets):
        """"""
        no_synthesized_mets_dic = {}
        with model:
            for metId in no_synthesized_mets:
                objectiveId = add_demand(model_info, model, metId)
                model.objective = objectiveId
                now_flux = model.slim_optimize()   
                no_synthesized_mets_dic[metId] = round(now_flux,3)
        return no_synthesized_mets_dic

    def reducing_carbon_source_input(self, model_info, model, check_model):
        """"""
        model_info["no_gap_met_input"] = []
        for rxn in check_model.reactions:
            if len(rxn.metabolites) == 1:
                r_met = list(rxn.metabolites.keys())[0]
                if r_met.id in model_info['no_gap_metid']:
                    # if (rxn.reactants and rxn.lower_bound < 0) or (rxn.products and rxn.upper_bound > 0):
                        print('reducing_carbon_source_input',model.reactions.get_by_id(rxn.id).bounds,rxn.bounds,rxn.id,r_met.id,rxn)
                        # model.reactions.get_by_id(rxn.id).bounds = rxn.bounds
                        model.reactions.get_by_id(rxn.id).bounds = (-1000,1000)
                        model_info["no_gap_met_input"].append(r_met.id)
                        model_info["no_gap_rxn_input"].append(rxn.id)
        for rxn in model_info["reactions"]:
            if rxn['id'] in model_info["no_gap_rxn_input"]:
                rxn['bounds'] = (-1000,1000) 
        print('biomass',model.slim_optimize())
        # write_flux_file(model_info, model, 'bio1', 'test_biomasss', self.cfg)
        # # small_demand_id = add_demand(model_info, model, 'cpd00166_c0')
        # mets = model.metabolites.get_by_id('cpd00166_c0') 
        # demand_rxn = model.add_boundary(mets, type = 'demand')
        # model.objective = demand_rxn.id
        # print(demand_rxn)
        # print(demand_rxn,model.reactions.get_by_id(demand_rxn.id).bounds)
        # print('cpd00166_c0',model.slim_optimize())
        # # sink_rxn = self.add_sink_rxn(model, model_info)
        # print('cpd00166_c0',model.slim_optimize())
        # exit()

    def update_no_gap_metid(self, model_info, model, no_synthesized_mets, flag):
        """"""
        no_gap = []
        with model:
            for metId in no_synthesized_mets:
                if flag == 0:
                    objectiveId = add_demand(model_info, model, metId)
                else:
                    mets = model.metabolites.get_by_id(metId)
                    try:
                        demand_rxn = model.add_boundary(mets, type = 'demand') # 添加demand反应
                        objectiveId = demand_rxn.id
                    except:
                        objectiveId = f"DM_{metId}"
                model.objective = objectiveId
                ini_flux = model.slim_optimize()     
                # print(metId,objectiveId,ini_flux,'ccccccccccccc') 
                if pd.isna(ini_flux) or ini_flux <= 1e-6:
                    no_gap.append(metId)
        model_info['no_gap_metid'] = no_gap
        return no_gap
    
    def find_min_carbon_source_set(self, model_info, model, check_model):
        """"""
        carbon_source_rxn, flag = [], 0
        for rxn in model_info["reactions"]:
            if 'carbon_source' in rxn.keys():
                carbon_source_rxn.append(rxn['id'])
        # model.reactions.get_by_id('EX_fe2(e)').bounds=(-1000,1000)
        # print('EX_fe2(e):',model.slim_optimize(),'........................................111.......................')
        min_set_dict = get_min_carbon_source_set(model_info, model, check_model, carbon_source_rxn, flag)
        print('min_set_dict',min_set_dict,'sssssssssssssssss')
        if min_set_dict:
            for k,v in min_set_dict.items():
                min_set_key = k
                min_set_all = v
            glucose_rxnId = model_info['glucose_rxnId']
            if min_set_all:
                if 'min_set' == min_set_key: # 这是第一步逐一打开获取到的列表，在列表中任取一个作为最小碳源都可以满足生长，先找葡萄糖、蔗糖等，没有就任选一个作为碳源
                    if glucose_rxnId in min_set_all:
                        min_set = glucose_rxnId
                    else:
                        carbon_rxnId = find_other_carbon_source(model_info, model)
                        min_set = carbon_rxnId if carbon_rxnId in min_set_all else min_set_all[0]  
                elif 'min_set1' == min_set_key: # 第二步的列表，其中包含多个子列表，每个子列表都是一对碳源，任选一对碳源打开都可以满足生长
                    min_set = min_set_all[0]
                else: # 第三步得到的列表，其中所有物质都需要打开最为最小碳源集
                    min_set = min_set_all
            if type(min_set) == str:
                model_info["min_set2"] = [min_set]
            else:
                model_info["min_set2"] = min_set
        else:
            model_info["min_set2"] = ""
        for rxn in model_info["reactions"]:
            for rxnId in model_info["min_set2"]:
                if rxn['id'] == rxnId:
                    rxns = model.reactions.get_by_id(rxnId)
                    if rxns.lower_bound == 0 or rxns.lower_bound < -100:
                        rxns.bounds = (-10,1000)
                    rxn['bounds'] = rxns.bounds # model_info中的反应信息同步更新


    def add_no_gap_sink_rxn(self, model_info, model):
        """"""
        initial_rxn_id = model_info["initial_rxn"]["biomass_rxn_id"]
        initial_rxn = model.reactions.get_by_id(initial_rxn_id)
        if len(model_info['no_gap_metid']) <= (len(initial_rxn.reactants)*0.2):
            for metId in model_info['no_gap_metid']:
                mets = model.metabolites.get_by_id(metId)
                model.add_boundary(mets, type = 'sink')
            

    def preprocess_initial_zero_bio(self, model_info, model, model_control_info, general_library, no_synthesized_mets, obj, check_model):
        """"""
        print('bio',model.slim_optimize())
        if model.slim_optimize() <= 1e-6:
            big_model = self.get_big_model(model_info, model, general_library)
            model_control_info['check_biomass_production']["score"] = 0
            no_synthesized_mets = self.gap_initial_zero_bio(model_info, model, model_control_info, general_library, no_synthesized_mets, big_model, obj, check_model)
            model.objective = model_info["initial_rxn"]["biomass_rxn_id"]
            print('bio22222',model.slim_optimize())
            if model.slim_optimize() <= 1e-6:
                no_gap = self.update_no_gap_metid(model_info, model, no_synthesized_mets, 0)
                self.reducing_carbon_source_input(model_info, model, check_model) # 将含有交换反应的物质打开
                if model.slim_optimize() <= 1e-6:
                    new_no_gap = self.update_no_gap_metid(model_info, model, no_synthesized_mets, 1)
                    self.find_min_carbon_source_set(model_info, model, check_model) # 将最小碳源集打开
                    print(model.slim_optimize(),'sssssssssssssss')
                    if model.slim_optimize() <= 1e-6:
                        self.add_no_gap_sink_rxn(model_info, model) # 添加所有不能合成物质的sink反应
                        print(model.slim_optimize(),'aaaaaaaaaaaaaaa',model_info['no_gap_metid'],new_no_gap,'dddddddd',no_gap)
                        if model.slim_optimize() <= 1e-6:
                            temp = [f"After gap, the model still cannot simulate growth, the biomass components {','.join(model_info['no_gap_metid'])} cannot be synthesized, please check the model boundary conditions and monomer synthesis pathway"]
                            write_biomass(model_control_info, temp, "Gapfilling_of_biomass_components")
                            return -1
                        else:
                            temp = [f"After gap, the model still cannot simulate growth, the biomass components {','.join(model_info['no_gap_metid'])} cannot be synthesized"]
                            temp.extend([f"But in the initial model, ({','.join(model_info['no_gap_met_input'])}) can be directly used as a carbon source input, so supply ({','.join(model_info['no_gap_met_input'])}) as carbon source,\n the model still cannot simulate growth, then open the minimum carbon source set ({','.join(model_info['min_set2'])}), the model still cannot simulate growth;\n so add the sink reaction of the non synthesized components ({','.join(model_info['no_gap_metid'])}) yields a biomass value of {round(model.slim_optimize(),3)}"])
                            write_biomass(model_control_info, temp, "Gapfilling_of_biomass_components")
                            model_info["initial_rxn"]["biomass_rxn_flux"] = model.slim_optimize()
                    else:
                        temp = [f"After gap, the model still cannot simulate growth, the biomass components {','.join(model_info['no_gap_metid'])} cannot be synthesized"]
                        temp.extend([f"But in the initial model, ({','.join(model_info['no_gap_met_input'])}) can be directly used as a carbon source input, so supply ({','.join(model_info['no_gap_met_input'])}) as carbon source, the model still cannot simulate growth,\n then open the minimum carbon source set ({','.join(model_info['min_set2'])}) yields a biomass value of {round(model.slim_optimize(),3)}"])
                        write_biomass(model_control_info, temp, "Gapfilling_of_biomass_components")
                        model_info["initial_rxn"]["biomass_rxn_flux"] = model.slim_optimize()
                else:
                    temp = [f"After gap, the model still cannot simulate growth, the biomass components {','.join(no_gap)} cannot be synthesized"]
                    temp.extend([f"But in the initial model, ({','.join(model_info['no_gap_met_input'])}) can be directly used as a carbon source input, so supply ({','.join(model_info['no_gap_met_input'])}) as carbon source yields a biomass value of {round(model.slim_optimize(),3)}"])
                    write_biomass(model_control_info, temp, "Gapfilling_of_biomass_components")
                    model_info["initial_rxn"]["biomass_rxn_flux"] = model.slim_optimize()
            else:
                no_synthesized_mets_dic = self.get_no_synthesized_mets_flux(model, model_info, no_synthesized_mets)
                temp = [f"After gapfilling, {','.join(no_synthesized_mets)} can be synthesized as follows: {no_synthesized_mets_dic}"]
                temp.extend([f"The growth rate of this model is {round(model.slim_optimize(),3)} after adding the above reaction"])
                write_biomass(model_control_info, temp, "Gapfilling_of_biomass_components")
                model_info["initial_rxn"]["biomass_rxn_flux"] = model.slim_optimize()
        else:
            temp = ["All of the biomass components can be synthesized."]
            write_biomass(model_control_info, temp, "Gapfilling_of_biomass_components")
        return 0
     


    def get_adp_h_pi(self, model):
        """"""
        adp_identifier, h_identifier, pi_identifier = '', '', ''
        for met in model.metabolites:
            if met.name in ADP_NAME and met.compartment in C_COMPARTMENT:
                adp_identifier = met.id 
            if met.name in H_NAME and met.compartment in C_COMPARTMENT:
                h_identifier = met.id 
            if met.name in PI_NAME and met.compartment in C_COMPARTMENT:
                pi_identifier = met.id 
        return adp_identifier, h_identifier, pi_identifier


    def is_missing_sustaining_energy(self, model_info, model, model_control_info, initial_rxn_id):
        """"""
        initial_rxns = model.reactions.get_by_id(initial_rxn_id)
        bio_text = initial_rxns.build_reaction_string(use_metabolite_names=False)
        missing_identifier, sustain, is_missing = [], [], 0
        products_mets = [met.id for met in model.reactions.get_by_id(initial_rxn_id).products]
        adp_identifier, h_identifier, pi_identifier = self.get_adp_h_pi(model)
        if len(set(ADP+H+PI) & set(products_mets)) == 0:
            sustain = ["Growth-related maintenance (h, adp, pi) in biomass is missing, please make up"]
            is_missing += 1
        elif len(set(ADP) & set(products_mets)) == 0 or len(set(H) & set(products_mets)) == 0 or len(set(PI) & set(products_mets)) == 0:
            if len(set(ADP) & set(products_mets)) == 0:  
                if len(set(H) & set(products_mets)) != 0:
                    coefficient = initial_rxns.get_coefficient(h_identifier)
                if len(set(PI) & set(products_mets)) != 0:
                    coefficient = initial_rxns.get_coefficient(pi_identifier)
                initial_rxns.reaction = bio_text + f' + {coefficient} {adp_identifier}'
                missing_identifier.append(adp_identifier)
            elif len(set(H) & set(products_mets)) == 0:
                if len(set(ADP) & set(products_mets)) != 0:
                    coefficient = initial_rxns.get_coefficient(adp_identifier)
                if len(set(PI) & set(products_mets)) != 0:
                    coefficient = initial_rxns.get_coefficient(pi_identifier)
                missing_identifier.append(h_identifier)
                initial_rxns.reaction = bio_text + f' + {coefficient} {h_identifier}'
            elif len(set(PI) & set(products_mets)) == 0:
                if len(set(H) & set(products_mets)) != 0:
                    coefficient = initial_rxns.get_coefficient(h_identifier)
                if len(set(ADP) & set(products_mets)) != 0:
                    coefficient = initial_rxns.get_coefficient(adp_identifier)
                missing_identifier.append(pi_identifier)
                initial_rxns.reaction = bio_text + f' + {coefficient} {pi_identifier}'
            sustain = [f"The growth-related maintenance energy in biomass is incorrect and missing, and {','.join(missing_identifier)} is added to the product side of biomass equation."]
            is_missing += 1
        if len(sustain) != 0:
            write_biomass(model_control_info, sustain, "Check_synthesizability_of_biomass_components")
        initial_bio_value = model.slim_optimize()
        initial_biomass_rxn = initial_rxns.build_reaction_string(use_metabolite_names=False)
        if is_missing != 0:
            model_control_info['check_biomass_production']["score"] = 0
            ini_bio_value = model_info["initial_rxn"]["biomass_rxn_flux"]              
            sustain_result = [f"Initial biomass growth rate (before improving Growth-Related Maintenance Energy): {round(ini_bio_value,5)}"]
            sustain_result.extend([f"Initial biomass equation (before improving Growth-Related Maintenance Energy): {bio_text}"])
            sustain_result.extend([f"Initial biomass growth rate (after improving Growth-Related Maintenance Energy): {round(initial_bio_value,5)}"])
            sustain_result.extend([f"Initial biomass equation (after improving Growth-Related Maintenance Energy): {initial_biomass_rxn}"])
        else:
            sustain_result = [f"Initial biomass growth rate (growth related maintenance can be intact): {round(initial_bio_value,5)}"]
            sustain_result.extend([f"Initial biomass equation (growth related maintenance energy intact): {initial_biomass_rxn}"])
        return sustain_result


    def get_macromolecule(self, model, initial_rxn_id, model_control_info):
        """"""
        macromolecule, special_macro = [], ['dna', 'rna', 'protein', 'lipid', 'peptido', 'cellwall']
        for metId in model.reactions.get_by_id(initial_rxn_id).reactants:
            met = model.metabolites.get_by_id(str(metId))
            if pd.isna(met.formula) or 'R' in met.formula or 'X' in met.formula or met.formula == 'nan' or not met.formula or met.formula == 'null':
                metName = met.name
                if any(metName.lower() in macro for macro in MACROMOLECULE) :  # 只要当前代谢物name在大分子列表里，就找到了大分子
                    macromolecule.append(str(metId))
                if any(macro in metName.lower() or macro in str(metId).lower() for macro in special_macro):
                    macromolecule.append(str(metId))
        macromolecule = list(set(macromolecule))
        temp = [f"Macromolecules include the following components: {macromolecule}"]
        write_biomass(model_control_info, temp, "Check_macromolecules")
        return macromolecule
    

    def get_smallmolecule(self, model_info, model, initial_rxn_id, macromolecule, model_control_info):
        """"""
        smallmolecule = []
        initial_rxn = model.reactions.get_by_id(initial_rxn_id)
        # products_mets = [met.id for met in initial_rxn.products]
        for met in initial_rxn.metabolites:
            if str(met.id) not in macromolecule and met not in initial_rxn.products and met.name not in SMALLMOLECULE:
                smallmolecule.append(str(met.id))
        temp = [f"Small molecule metabolites include : {smallmolecule}"]
        write_biomass(model_control_info, temp, "Check_small_molecules")
        return smallmolecule


    def set_atpm_bounds(self, model_info, model):
        """"""
        for atpm_id in ATPM:
            if atpm_id in model_info['all_rxn_id']:
                model.reactions.get_by_id(atpm_id).bounds = (0,1000)
                for rxn in model_info['reactions']:
                    rxn['atpms_modify'] = "false"
                    if rxn['id'] == atpm_id:
                        rxn['atpms_modify'] = "true"

    def get_formula(self, model, ids, modelInfo, check_model, obj):
        """"""
        readict, prodict, resultdict ,flag, delect_rxn, need_add_rxn = {}, {}, {}, 0, '', []
        mets = model.metabolites.get_by_id(ids)
        if pd.isna(mets.formula) or 'R' in mets.formula or 'X' in mets.formula or mets.formula == 'nan' or not mets.formula or mets.formula == 'null':
            with model:
                # demand_rxn = model.add_boundary(mets, type = 'demand')
                if f'SK_{ids}' in model.reactions:
                    model.reactions.get_by_id(f'SK_{ids}').bounds = (0,0)
                demand_rxn = add_demand(modelInfo, model, ids)
                model.objective = demand_rxn
                print(mets.id,mets.formula,model.slim_optimize(),'.................')
                if model.slim_optimize() > 0:
                    pfba_solution = pfba(model)  # 限定通量值v大于1e-6的反应
                    need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]
                else:
                    need_fluxes = {}
               
                # for r,v in need_fluxes.items():
                #     rxns = model.reactions.get_by_id(r) 
                #     if len(rxns.metabolites)==1:
                #         all_met = [m for m in rxns.metabolites][0]
                #         products_mets = [m.id for m in rxns.products]
                #         reactants_mets = [m.id for m in rxns.reactants]
                #         if (ids in reactants_mets and v < 0) or (ids in products_mets and v > 0):
                #             model.reactions.get_by_id(r).bounds = (0,0)
                #             delect_rxn = r
                #             print('delect_rxn:',delect_rxn,'.........................................................')
                #             print(mets.id,mets.formula,model.slim_optimize(),'关掉生成目标物质的反应后,速率变为上述值')
                #             if model.slim_optimize() <= 1e-6:
                #                 biomasses = Biomasses(self.cfg)
                #                 general_library = biomasses.get_general_library() 
                #                 new_way_rxn = biomasses.gapfilling(model, modelInfo, demand_rxn, general_library)
                #                 # model.reactions.get_by_id(r).bounds = (0,1000)
                #                 pfba(model)
                #                 need_add_rxn, gap_info = biomasses.add_gapfilling_rxn(model, check_model, general_library, new_way_rxn, obj) 
                #             if need_add_rxn:
                #                 biomasses.add_gap_rxn(model, need_add_rxn, general_library)
                #             model.reactions.get_by_id(r).bounds = (0,1000)
                #             pfba(model)
                #             if model.slim_optimize() > 0: # 上面删掉了一个反应，如果没有gap上，此时计算pfba就会报错：'NoneType' object has no attribute 'name'
                #                 pfba_solution = pfba(model)  # 限定通量值v大于1e-6的反应
                #                 need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]
                #             else:
                #                 need_fluxes = {}
                #             break 
                # write_flux_file(modelInfo, model, demand_rxn, ids, self.cfg)
                for r,v in need_fluxes.items():
                    if r==demand_rxn:
                        continue
                    rxns = model.reactions.get_by_id(r) 
                    if len(rxns.metabolites)==1:
                        all_met = [m for m in rxns.metabolites][0]
                        products_mets = [m.id for m in rxns.products]
                        reactants_mets = [m.id for m in rxns.reactants]
                        if pd.isna(all_met.formula) or 'R' in all_met.formula or 'X' in all_met.formula or all_met.formula == 'nan' or not all_met.formula or all_met.formula == 'null':
                            flag += 1
                        if flag > 1:
                            return 0
                        # print(rxns,v,all_met.formula,all_met.elements.items())
                        if (reactants_mets and v < 0) or (products_mets and v > 0):
                            for key, value in all_met.elements.items():
                                if key in readict: # 将两个字典合并并对相同键的值进行求和
                                    readict[key]+= value
                                else:
                                    readict[key]=value
                        if (reactants_mets and v > 0) or (products_mets and v < 0):
                            print(reactants_mets,products_mets,v,all_met.formula,'xxx')
                            for key, value in all_met.elements.items():
                                if key in prodict:
                                    prodict[key]+= value
                                else:
                                    prodict[key]=value
            print(readict,prodict)
            # if need_add_rxn:
            #     print('need_add_rxn: 关掉生成目标物质的反应后,速率为0,进行gap,gap上了这些反应:',need_add_rxn)
            #     biomasses.add_gap_rxn(model, need_add_rxn, general_library)
            #     model.reactions.get_by_id(delect_rxn).bounds = (0,0)
            for key in OrderedDict.fromkeys(list(readict.keys()) + list(prodict.keys())): # 将两个字典的所有键合并，使用OrderedDict来保持顺序
                resultdict[key] = readict.get(key,0) - prodict.get(key,0)
            print(resultdict)
            resultdict = {key:value for key,value in resultdict.items() if value !=0} # 移除值为0的键值对
            print(resultdict)
            result_string = ''.join(f'{key}' if value == 1 else f'{key}{value}' for key, value in resultdict.items())
            print(result_string,'xxx')
            print(mets.id,mets.formula,model.slim_optimize(),'.................')
            if not result_string:
                return 0
            print(result_string)
            mets.formula = result_string
            print(mets.id,mets.formula,model.slim_optimize(),'.................')
            print(mets.id,mets.formula,result_string)
       


    def calculateResult(self, model, initial_rxn_id, model_info, sink_rxn, check_model, obj):
        """"""
        result, final = 0, 0
        # no_use = ['rxn13783_c0','rxn13784_c0','rxn13782_c0'] # 这几个物质不进行计算
        if model.slim_optimize() > 0:
            pfba_solution = pfba(model)
            need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]
        else:
            need_fluxes = {}
        # need_fluxes=model.optimize().fluxes[abs(model.optimize().fluxes)>1e-6]
        for ids,v in need_fluxes.items():
            v = round(v,5)
            rxns = model.reactions.get_by_id(ids)
            all_mets = [met.id for met in rxns.metabolites]
            all_mets_name = [met.name for met in rxns.metabolites]
            if len(all_mets) == 1 and ids != initial_rxn_id and 'biomass' not in all_mets_name[0].lower() and ids not in sink_rxn and ids not in model_info["no_gap_rxn_input"]: 
                mets = model.metabolites.get_by_id(all_mets[0])
                if pd.isna(mets.formula) or 'R' in mets.formula or 'X' in mets.formula or mets.formula == 'nan' or not mets.formula or mets.formula == 'null':
                    met_formula = self.get_formula(model, all_mets[0], model_info, check_model, obj)
                    print(model_info["no_gap_rxn_input"],ids,met_formula,'ddddddddddddddddddddd')
                    if met_formula == 0:
                        return all_mets[0]
                result += v * relative_molecular_mass(model, all_mets[0])
        try :
            final = abs(result/model.slim_optimize())
        except ZeroDivisionError:
            final = 0
        return final


    def get_initial_biomass(self, model_info, model, initial_rxn_id, model_control_info, check_model, obj):
        """"""
        self.set_atpm_bounds(model_info, model)
        initial_biomass = self.calculateResult(model, initial_rxn_id, model_info, [], check_model, obj)
        model_info['biomass']['initial_biomass'] = initial_biomass
        if type(initial_biomass) != str and (990 < initial_biomass < 1010):
            temp = [f"The biomass relative molecular mass is {round(initial_biomass,3)}mg, corrected"]
            write_biomass(model_control_info, temp, "Gapfilling_of_biomass_components")
            # return 0



    def smallmolecule_iscoupling(self, model, model_info, small_demand_id, initial_rxn_id, metId, general_library, model_control_info, obj, check_model):
        """"""
        iscoupling, close_coupling_rxn, need_add_rxn, temp_list = 0, [], [], []
        sink_rxn = model_info["sink_rxn"]
        initial_flux = model.slim_optimize()
        print(initial_flux,type(initial_flux))
        if pd.isna(initial_flux) or initial_flux <= 1e-6 :
            temp_list = [f"{metId} rate is 0, so gap:"] 
            new_way_rxn = self.gapfilling(model, model_info, small_demand_id, general_library)
            need_add_rxn, gap_info = self.add_gapfilling_rxn(model, model_info, check_model, general_library, new_way_rxn, obj) 
            model_info["bio_modified"].extend(need_add_rxn)
            temp_list.extend(gap_info)
            write_biomass(model_control_info, temp_list, "Check_small_molecules")
            if len(need_add_rxn) == 0:
                temp_list = [f"In biomass generation, {metId} must be coupled with macromolecules, the gapfilling cannot be performed, and there is an error! ! ! Please correct and try again"]
                write_biomass(model_control_info, temp_list, "Check_small_molecules")
                return -1
            elif pd.isna(model.slim_optimize()):
                temp_list = [f"The compositing rate is nan after {metId} gapfilling, there is an error! ! ! Please correct and try again"]
                write_biomass(model_control_info, temp_list, "Check_small_molecules")
                return -1
            else:
                final_flux = model.slim_optimize()
                temp_list = [f"The composite rate after {metId} gapfilling is: {round(final_flux,5)}"]
                write_biomass(model_control_info, temp_list, "Check_small_molecules")
                return need_add_rxn, close_coupling_rxn
        while True :
            status = 0  # status状态符：只要本次循环没有发现biomass和大分子耦合反应，就退出循环
            if model.slim_optimize() > 0:
                pfba_solution = pfba(model)
                need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]
                # need_fluxes=model.optimize().fluxes[abs(model.optimize().fluxes)>1e-6]
            else:
                need_fluxes = {}
            for r in list(need_fluxes.keys()): # 如果biomass反应也在结果中，就要解除大分子与biomass的耦联关系
                rxns = model.reactions.get_by_id(r)
                all_mets = [m.id for m in rxns.metabolites]
                met_formula = list(rxns.metabolites.keys())[0].formula
                if r == initial_rxn_id or (len(all_mets) == 1 and (pd.isna(met_formula) or 'R' in met_formula or 'X' in met_formula or met_formula == 'nan' or not met_formula or met_formula == 'null')):
                    if r != small_demand_id and r not in sink_rxn:  # 判断当前结果文件含有biomass耦合反应 或者 交换反应中大分子耦合反应
                        print(r,sink_rxn,'xxxx')
                        met_formula = self.get_formula(model, all_mets[0], model_info, check_model, obj)
                        if met_formula != 0:
                            continue
                        status +=  1
                        iscoupling += 1
                        temp_list = [f"{metId} synthesis rate is: {round(model.slim_optimize(),5)}"]
                        model.reactions.get_by_id(r).bounds = (0,0)
                        close_coupling_rxn.append(r)
                        temp_list.extend([f"{small_demand_id} synthesis is coupled with {r}, close {r}"])
                        write_biomass(model_control_info, temp_list, "Check_small_molecules")
            if status == 0:
                break
        if iscoupling == 0 and initial_flux > 1e-6:
            return 0  
        else:
            final_flux2 = model.slim_optimize()
            if final_flux2 <= 1e-6 or pd.isna(final_flux2):
                temp_list = [f"The macromolecules coupled with the small molecule {metId} are all closed, the generation rate is 0, and the gapfilling is performed:"]
                new_way_rxn = self.gapfilling(model, model_info, small_demand_id, general_library)
                need_add_rxn, gap_info = self.add_gapfilling_rxn(model, model_info, check_model, general_library, new_way_rxn, obj) 
                model_info["bio_modified"].extend(need_add_rxn)
                temp_list.extend(gap_info)
                if len(need_add_rxn) == 0:
                    temp_list.extend([f"In biomass generation, {metId} must be coupled with macromolecules, the gapfilling cannot be performed, and there is an error! ! ! Please correct and try again"])
                    write_biomass(model_control_info, temp_list, "Check_small_molecules")
                    return -1
                else :
                    final_flux = model.slim_optimize()
                    temp_list.extend([f"The synthesis rate after gapfilling is: {round(final_flux,5)}"])
                    write_biomass(model_control_info, temp_list, "Check_small_molecules")
                    return need_add_rxn, close_coupling_rxn
            else:
                temp_list = [f"The macromolecules coupled with the small molecule {metId} are all closed, and the synthesis rate is {round(final_flux2,5)}, if it is not 0, no gap will be performed"]
                write_biomass(model_control_info, temp_list, "Check_small_molecules")
                return need_add_rxn, close_coupling_rxn


    def small_fix(self, model_info, model, initial_rxn_id, smallmolecule, general_library, model_control_info, check_model, obj):
        """"""
        right_small_met, small_need_add_rxn, close_coupling_rxn = [], [], []
        for metId in smallmolecule:
            print(metId,'xxxxxxxxxx')
            with model:
                small_demand_id = add_demand(model_info, model, str(metId))
                model.objective = small_demand_id
                # write_flux_file(model_info, model, small_demand_id, metId, self.cfg)
                small_info = self.smallmolecule_iscoupling(model, model_info, small_demand_id, initial_rxn_id, metId, general_library, model_control_info, obj, check_model)
                if small_info == -1 :
                    continue
                elif small_info == 0: 
                    right_small_met.append(metId)
                else:
                    small_need_add_rxn.extend(small_info[0])
                    close_coupling_rxn.extend(small_info[1])
            if close_coupling_rxn :
                model_info['biomass']['small_macro_modify'] = "yes"
                model_control_info['check_biomass_production']["score"] = 0
                for r in close_coupling_rxn:
                    model.reactions.get_by_id(r).bounds = (0,0)
            if small_need_add_rxn:
                model_control_info['check_biomass_production']["score"] = 0
                model_info['biomass']['small_macro_modify'] = "yes" 
                self.add_gap_rxn(model, small_need_add_rxn, general_library, check_model)
        temp = [f"The synthesis of small molecule metabolites {right_small_met} were not coupled to macromolecules, and the synthetic pathway were correct"]
        write_biomass(model_control_info, temp, "Check_small_molecules")
        model.reactions.get_by_id(initial_rxn_id).bounds = (0,1000)


    def Normalized_macromolecule(self, model, model_info, metId, initialMacromoleculeResult, macro_demand_id, initial_rxn_id, normalized_macro_time, model_control_info, general_library, obj, check_model):
        """"""
        close_coupling_rxn, initial_macromolecule_rxn, macromoleculeRxnId, need_add_rxn, close_coupling_rxn = [], '', '', [], []
        sink_rxn = model_info["sink_rxn"]
        mets = model.metabolites.get_by_id(metId)
        print(model.slim_optimize())
        print(macro_demand_id)
        print(model.reactions.get_by_id(macro_demand_id))
        if model.slim_optimize() <= 1e-6:
            temp = [f"error : Macromolecules {metId} cannot be synthesized! Unable to normalize"]
            write_biomass(model_control_info, temp, "Check_macromolecules")
            return 0
        while initialMacromoleculeResult :
            status = 0  # status状态符：只要本次循环没有发现biomass和大分子耦合反应，就退出循环    
            if model.slim_optimize() <= 1e-6 or model.slim_optimize() == 'nan':
                temp = [f"error : After the gap addition reaction, Macromolecules {metId} cannot be synthesized! Unable to normalize"]
                write_biomass(model_control_info, temp, "Check_macromolecules")
                return 0
            print('test:......................',model.slim_optimize())
            # write_final_model(model,self.cfg)
            if model.slim_optimize() > 0:
                pfba_solution = pfba(model)
                need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6]
                # need_fluxes=model.optimize().fluxes[abs(model.optimize().fluxes)>1e-6]
            else:
                need_fluxes = {}
            for r in list(need_fluxes.keys()): # 如果biomass反应也在结果中，就要解除大分子与biomass的耦联关系
                rxns = model.reactions.get_by_id(r)
                all_mets = [m.id for m in rxns.metabolites]
                met_formula = list(rxns.metabolites.keys())[0].formula
                if r == initial_rxn_id or (len(all_mets) == 1 and (pd.isna(met_formula) or 'R' in met_formula or 'X' in met_formula or met_formula == 'nan' or not met_formula or met_formula == 'null')):
                    if r != macro_demand_id and r not in sink_rxn:  # 判断当前结果文件含有biomass耦合反应 或者 交换反应中大分子耦合反应
                        print(r,sink_rxn,'xxxxxx')
                        met_formula = self.get_formula(model, all_mets[0], model_info, check_model, obj)
                        if met_formula != 0:
                            continue
                        status +=  1
                        model.reactions.get_by_id(r).bounds = (0,0)
                        close_coupling_rxn.append(r)
                        temp_Macromolecule = self.calculateResult(model, macro_demand_id, model_info, sink_rxn, check_model, obj) # 有耦合反应就把计算结果值赋给temp_Macromolecule
                        temp = [f"{macro_demand_id} composition is coupled with {r}:{rxns}"]
                        print('temp_Macromolecule: ',temp_Macromolecule)
                        if type(temp_Macromolecule) == str:
                            pass
                        else:
                            temp.extend([f"After closing {r}, the relative molecular mass of {macro_demand_id} is {round(temp_Macromolecule,5)}"])
                            write_biomass(model_control_info, temp, "Check_macromolecules")
                        if temp_Macromolecule == 0 :  # 判断当前结果是否为0，为0说明没有合成途径了，就gapfilling添加上途径
                            temp = [f"All macromolecules are closed, the relative molecular mass of {macro_demand_id} is 0, and the gap is performed:"]
                            new_way_rxn = self.gapfilling(model, model_info, macro_demand_id, general_library)
                            need_add_rxn, gap_info = self.add_gapfilling_rxn(model, model_info, check_model, general_library, new_way_rxn, obj) 
                            model_info["bio_modified"].extend(need_add_rxn)
                            temp.extend(gap_info)
                            write_biomass(model_control_info, temp, "Check_macromolecules")
                            if len(need_add_rxn) == 0:
                                temp = [f"In biomass generation, {metId} must be coupled with macromolecules, the gapfilling cannot be performed, there is an error, please correct"]
                                write_biomass(model_control_info, temp, "Check_macromolecules")
                                return 0
                            break
                        else:  # 不为0时，可能是有途径，进入下一次循环，status为0退出循环；可能是还有其他的大分子反应，进入下一轮关
                            initialMacromoleculeResult = temp_Macromolecule
            if status == 0:
                break
        if model.slim_optimize() > 0:
            pfba_solution = pfba(model)
            need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6] 
        else:
            need_fluxes = {}
        for r,v in need_fluxes.items():  # 获取合成大分子的反应以及代谢物系数
            rxns = model.reactions.get_by_id(r)
            products_mets = [m.id for m in rxns.products]
            reactants_mets = [m.id for m in rxns.reactants]
            coefficient_dict = {}
            if (mets.id in products_mets and v > 0) or (mets.id in reactants_mets and v < 0):
                macromoleculeRxnId = rxns.id
                model_info["bio_modified"].append(macromoleculeRxnId)
                initial_macromolecule_rxn = model.reactions.get_by_id(macromoleculeRxnId)
                for j in rxns.metabolites: # rxn是合成大分子的反应
                    # met_coefficient[str(j)] = rxn.get_coefficient(j)
                    if str(j) == mets.id: # 目标大分子系数不变
                        continue
                    coefficient_dict[j] = rxns.get_coefficient(j) * (1000/initialMacromoleculeResult) - rxns.get_coefficient(j) # 减去自己那一份，加回去就是所需的倍数
                rxns.add_metabolites(coefficient_dict)
                break 
        if type(initial_macromolecule_rxn) == str:
            temp = [f"error : After removing the coupled macromolecules, Macromolecules {metId} cannot be synthesized! Unable to normalize"]
            write_biomass(model_control_info, temp, "Check_macromolecules")
            return 0
        initial_macromolecule_rxn, initial_macromolecule_rxn_id = keep_decimal_places(initial_macromolecule_rxn)
        if normalized_macro_time > 1:
            temp = [f"Round {normalized_macro_time}: {metId} Equation: {initial_macromolecule_rxn}"]
            if status != 0 :
                temp.extend([f"Round {normalized_macro_time}: {metId} coupling macromolecule {','.join(close_coupling_rxn)}, and {','.join(close_coupling_rxn)} formula is not clear!!!, unable to calculate relative to initial {metId} Molecular mass"])
            else :
                temp.extend([f"Round {normalized_macro_time}: {mets.id} relative molecular mass results: {round(initialMacromoleculeResult,5)}"])
        else :
            temp = [f"Initial equation: {initial_macromolecule_rxn}"]
            if status != 0 :
                temp.extend([f"{metId} coupling macromolecule {','.join(close_coupling_rxn)}, and {','.join(close_coupling_rxn)} formula is not clear!!!, the relative molecular mass of the initial {metId} cannot be calculated"])
            else :
                temp.extend([f"Initial relative molecular mass: {round(initialMacromoleculeResult,5)}"])
        write_biomass(model_control_info, temp, "Check_macromolecules")
        return macromoleculeRxnId, need_add_rxn, close_coupling_rxn



    def recover_macromolecule_rxn(self, model, model_info, macromoleculeRxn, need_add_rxn, general_library, need_close_coupling_rxn, initial_rxn_id, model_control_info, check_model):
        """"""
        for k,v in macromoleculeRxn.items():
            model_control_info['check_biomass_production']["score"] = 0
            model_info['biomass']['small_macro_modify'] = "yes" 
            model.reactions.get_by_id(k).reaction = v
        if need_add_rxn:
            model_control_info['check_biomass_production']["score"] = 0
            model_info['biomass']['small_macro_modify'] = "yes" 
            self.add_gap_rxn(model, need_add_rxn, general_library, check_model)
        for r in need_close_coupling_rxn:
            model_control_info['check_biomass_production']["score"] = 0
            model_info['biomass']['small_macro_modify'] = "yes" 
            model.reactions.get_by_id(r).bounds = (0,0)
        model.reactions.get_by_id(initial_rxn_id).bounds = (0,1000)



    def macro_fix(self, model_info, model, initial_rxn_id, macromolecule, general_library, model_control_info, check_model, obj):
        """"""
        gapfilling_rxn, need_close_coupling_rxn, macromoleculeRxn = [], [], {}
        sink_rxn = model_info["sink_rxn"]
        for metId in macromolecule:
            temp = [f"{metId}:"]
            write_biomass(model_control_info, temp, "Check_macromolecules")
            normalized_macro_time, while_time = 0, 0
            with model:
                macro_demand_id = add_demand(model_info, model, str(metId))
                model.objective = macro_demand_id
                print(model.slim_optimize(),macro_demand_id,'xxx')
                try :
                    model.reactions.get_by_id(macro_demand_id).bounds = check_model.reactions.get_by_id(macro_demand_id).bounds # 大分子反应边界还原
                except :
                    print(f'{model.id}没有{macro_demand_id}反应,不需要还原')
                print(f'{metId}_objective_value: ',model.slim_optimize())
                # write_final_model(model,self.cfg)
                initialMacromoleculeResult = self.calculateResult(model, macro_demand_id, model_info, sink_rxn, check_model, obj)
                print(model.slim_optimize(),macro_demand_id,'xxx')
                if initialMacromoleculeResult == 0:
                    continue
                if type(initialMacromoleculeResult) != str and (990 < initialMacromoleculeResult < 1010) :
                    continue
                else:
                    normalized_macro_time += 1
                    print(model.slim_optimize(),macro_demand_id,'xxx')
                    macromoleculeInfo = self.Normalized_macromolecule(model, model_info, metId, initialMacromoleculeResult, macro_demand_id,initial_rxn_id,normalized_macro_time, model_control_info, general_library, obj, check_model) 
                    if macromoleculeInfo == -1 :
                        return -1
                    elif macromoleculeInfo == 0:
                        continue
                    else :
                        macromoleculeRxnId, need_add_rxn, close_coupling_rxn = macromoleculeInfo
                    finalMacromoleculeResult = self.calculateResult(model, macro_demand_id, model_info, sink_rxn, check_model, obj)
                    gapfilling_rxn.extend(need_add_rxn)
                    need_close_coupling_rxn.extend(close_coupling_rxn)
                    while type(finalMacromoleculeResult) == str or finalMacromoleculeResult < 990 or finalMacromoleculeResult > 1010 :
                        normalized_macro_time += 1
                        macromoleculeInfo = self.Normalized_macromolecule(model, model_info, metId, finalMacromoleculeResult, macro_demand_id,initial_rxn_id,normalized_macro_time, model_control_info, general_library, obj, check_model)
                        if macromoleculeInfo == -1 :
                            return -1
                        elif macromoleculeInfo == 0:  # 只是退出了第一个while
                            while_time += 1
                            break
                        else :
                            macromoleculeRxnId, need_add_rxn, close_coupling_rxn = macromoleculeInfo
                        finalMacromoleculeResult = self.calculateResult(model, macro_demand_id, model_info, sink_rxn, check_model, obj)
                        gapfilling_rxn.extend(need_add_rxn)
                        need_close_coupling_rxn.extend(close_coupling_rxn)
                    if while_time == 1 :
                        continue
                    finalRxn = model.reactions.get_by_id(macromoleculeRxnId)
                    finalRxn, finalRxn_id = keep_decimal_places(finalRxn)
                    final_rxn_id = model.reactions.get_by_id(macromoleculeRxnId).build_reaction_string()
                    macromoleculeRxn[macromoleculeRxnId] = final_rxn_id
                    temp = [f"Corrected equation: {finalRxn}"]
                    temp.extend([f"Corrected relative molecular mass: {round(finalMacromoleculeResult,5)}"])
                    write_biomass(model_control_info, temp, "Check_macromolecules")
            if need_close_coupling_rxn:
                for rxnId in need_close_coupling_rxn:
                    model.reactions.get_by_id(rxnId).bounds = (0,0)
            if gapfilling_rxn:
                self.add_gap_rxn(model, gapfilling_rxn, general_library, check_model)
        self.recover_macromolecule_rxn(model, model_info, macromoleculeRxn, gapfilling_rxn, general_library, need_close_coupling_rxn, initial_rxn_id, model_control_info, check_model)



    def Normalized_biomass(self, model, initial_rxn_id, final_biomass):
        """"""
        coefficient_dict = {}
        if final_biomass == 0:
            return -1
        rxns = model.reactions.get_by_id(initial_rxn_id)
        for j in rxns.metabolites:
            if (j in rxns.products) and ('biomass' in str(j).lower() or str(j).lower() in 'biomass' or j.name in BIOMASS): # 如果方程有biomass_c，那么biomass_c系数不变##########################改biomass列表
                continue
            coefficient_dict[j] = rxns.get_coefficient(j) * (1000/final_biomass) - rxns.get_coefficient(j)
        rxns.add_metabolites(coefficient_dict)
        print(rxns.reaction)
        return 0


    def biomass_fix(self, model_info, model, initial_rxn_id, model_control_info, check_model, obj):
        """"""
        with model:
            middle_bio_value = model.slim_optimize()
            print(model.slim_optimize(),'xxxxxx')
            sink_rxn = model_info["sink_rxn"]
            # temp = [f"When checking if the biomass is 1g, the gap added the following sink reaction: {','.join(sink_rxn)}"]
            # write_biomass(model_control_info, temp, "Normalization_of_biomass")
            print(model.slim_optimize(),'xxxxxxxx')
            # write_flux_file(model_info, model, initial_rxn_id, 'biomasss', self.cfg)
            modify_nor_bio = final_biomass = self.calculateResult(model, initial_rxn_id, model_info, sink_rxn, check_model, obj)
            modify_nor_bio_eq = model.reactions.get_by_id(initial_rxn_id).build_reaction_string(use_metabolite_names=False)
            model_info['biomass']['middle_bio_value'] = middle_bio_value
            model_info['biomass']['modify_nor_bio_eq'] = modify_nor_bio_eq
            model_info['biomass']['modify_nor_bio'] = modify_nor_bio
            if type(final_biomass) != str and (990 < final_biomass < 1010):
                temp = [f"The biomass is 1g after treatment with small molecules and macromolecules"]
            else:
                model_control_info['check_biomass_production']["score"] = 0
                temp = [f"When verifying the biomass weight of 1 g. The biomass did not reach 1g."]
            write_biomass(model_control_info, temp, "Normalization_of_biomass")
            if type(final_biomass) == str:
                model_control_info['check_biomass_production']["score"] = 0
                print(model.slim_optimize(),'sssssssssssssssssssss')
                temp = [f"The final biomass is still coupled with macromolecules {final_biomass}"]
                write_biomass(model_control_info, temp, "Normalization_of_biomass")
                return -1
            while final_biomass < 990 or final_biomass > 1010:
                get_final = self.Normalized_biomass(model, initial_rxn_id, final_biomass)
                print('get_final:.............',get_final)
                if get_final == -1:
                    model_control_info['check_biomass_production']["score"] = 0
                    temp = [f"The final flux final is 0, error"]
                    write_biomass(model_control_info, temp, "Normalization_of_biomass")
                    return -1
                final_biomass = self.calculateResult(model, initial_rxn_id, model_info, sink_rxn, check_model, obj)
                print('final_biomass:.................',final_biomass)
        final_biomass_equation = model.reactions.get_by_id(initial_rxn_id)
        final_biomass_equation, final_biomass_equation_ids = keep_decimal_places(final_biomass_equation)
        final_biomass_equation_id = model.reactions.get_by_id(initial_rxn_id).build_reaction_string()
        model_info['biomass']['final_biomass'] = final_biomass
        model_info['biomass']['final_biomass_equation'] = final_biomass_equation
        model.reactions.get_by_id(initial_rxn_id).reaction = final_biomass_equation_id
                

    def recover_atpm_bounds(self, model_info, model, check_model):
        """"""
        for rxn in model_info['reactions']:
            if 'atpms_modify' in rxn.keys() and rxn['atpms_modify'] == "true":
                model.reactions.get_by_id(rxn['id']).bounds == check_model.reactions.get_by_id(rxn['id']).bounds

    def convert_nan_to_null(self, model):
        """"""
        for met in model.metabolites:
            # if pd.isna(met.formula) or met.formula == 'nan' or not met.formula or met.formula == 'null' or met.formula is None:
            #     print(met,met.formula)
            #     met.formula = ''
            if pd.isna(met.formula) or met.formula == 'nan' or met.formula is None or not met.formula or met.formula == 'null':
                print(met,met.formula)
                met.formula = ''
            # if pd.isna(met.charge) or met.charge == 'nan':
            #     print(met,met.charge)
            #     met.charge = 0
        # for rxn in model.reactions:
        #     if pd.isna(rxn.bounds) or rxn.bounds == 'nan' or not rxn.bounds or rxn.bounds == 'null' or rxn.bounds is None or np.isinf(rxn.lower_bound) or np.isinf(rxn.upper_bound):
        #         print(rxn,rxn.bounds)
        #         rxn.bounds = ''
      
    


    def get_exchange(self, model, model_control_info):
        if model.slim_optimize() > 0:
            pfba_solution = pfba(model)
            need_fluxes = pfba_solution.fluxes[abs(pfba_solution.fluxes)>1e-6] 
        else:
            need_fluxes = {}
        # need_fluxes=model.optimize().fluxes[abs(model.optimize().fluxes)>1e-6]
        temp = [f"Model calculation conditions:"]
        for r,v in need_fluxes.items():
            rxn = model.reactions.get_by_id(r)
            all_mets = [m.id for m in rxn.metabolites]
            if len(all_mets) == 1:
                mets = model.metabolites.get_by_id(all_mets[0])
                temp.extend([f"{r}: {rxn.build_reaction_string(use_metabolite_names=True)}    {mets.formula}  {round(v,5)}  {rxn.bounds}"])
        write_biomass(model_control_info, temp, "model_calculation_conditions")


    def get_result(self, model_info, model, sustain_result, model_control_info):
        """"""
        final_bio_value = model.slim_optimize()
        initial_biomass = model_info['biomass']['initial_biomass']
        if type(initial_biomass) == str:
            sustain_result.extend([f"Biomass coupling macromolecule {initial_biomass}, and the formula of {initial_biomass} is unclear!!!, the relative molecular mass of the initial biomass cannot be calculated"])
        else:
            sustain_result.extend([f"Initial biomass relative molecular mass: {round(initial_biomass,5)}mg"])
        if 'small_macro_modify' in model_info['biomass'].keys():
            modify_nor_bio = model_info['biomass']['modify_nor_bio']
            modify_nor_bio_eq = model_info['biomass']['modify_nor_bio_eq']
            middle_bio_value = model_info['biomass']['middle_bio_value']
            if type(modify_nor_bio) == str:
                sustain_result.extend([f"Modifying the biomass equation after the reaction of small molecules and macromolecules: {modify_nor_bio_eq}"])
                sustain_result.extend([f"The biomass is coupled to the macromolecule {modify_nor_bio}, and the formula of {modify_nor_bio} is not clear!!!, it is impossible to calculate the relative molecular mass of the biomass after modifying the reaction of small molecules and macromolecules"])
            else:
                sustain_result.extend([f"Modifying the biomass equation after the reaction of small molecules and macromolecules: {modify_nor_bio_eq}"])
                sustain_result.extend([f"Modify the relative molecular mass results of biomass after the reaction of small molecules and macromolecules: {round(modify_nor_bio,5)}"])
            sustain_result.extend([f"After modifying the reaction of small molecules and macromolecules, the growth rate of biomass is: {round(middle_bio_value,5)}"])
        final_biomass_equation = model_info['biomass']['final_biomass_equation']
        final_biomass = model_info['biomass']['final_biomass'] 
        sustain_result.extend([f"After modifying the biomass, the final biomass equation: {final_biomass_equation}"])
        sustain_result.extend([f"After modifying the biomass, the final biomass relative molecular mass result: {round(final_biomass,5)}"])
        sustain_result.extend([f"The final relative molecular mass of biomass after modification is {round(final_biomass,5)}mg, and the final biomass growth rate is {round(final_bio_value,5)}"])
        write_biomass(model_control_info, sustain_result, "Normalization_of_biomass")
        model_info["bio_modified"].extend([model_info["initial_rxn"]["biomass_rxn_id"]])
        model_control_info["check_biomass_production"]["Final_biomass"] = f"Initial growth rate: {model_info['biomass_original']}; after setting the carbon source, the value becomes: {model_info['biomass_by_c_source']}; after setting special boundaries, the value becomes: {model_info['biomass_by_special_boundary']}; after quality control, final biomass growth rate: {round(final_bio_value,3)}.\nThe equations modified during biomass normalization correction are placed in the 'Model Revision' list, and you can make changes and edits by clicking the 'Edit' button."
        model_info["control_analysis"].append((model_info['biomass_original'],round(final_bio_value,3)))

    def modify_bios(self, model_info, model, model_control_info):
        """"""
        temp = []
        bio_modified = list(set(model_info["bio_modified"]))
        for rxnId in bio_modified:
            # print(rxnId,'xxxxxxxxxxxxxx')
            rxns = model.reactions.get_by_id(rxnId)
            equation, equation_id = keep_decimal_places(rxns)
            modify_rxn = {"reaction_id" : rxnId,
                          "equation(ID)" : equation_id, 
                          "equation(Name)" : equation,
                          "l_b" : rxns.lower_bound,
                          "u_b" : rxns.upper_bound,}
            temp.append(modify_rxn)
            # print(temp,'xxxxxxxxxxxxxxxxxx')
        model_control_info["check_biomass_production"]["model_revision"].extend(temp)


    def recover_special_boundary_rxn(self, model_info, model, check_model):
        """"""
        for rxn in model_info["reactions"]:
            if 'special_boundary' in rxn.keys():
                model.reactions.get_by_id(rxn['id']).bounds == check_model.reactions.get_by_id(rxn['id']).bounds
        for rxnId in model_info["sink_rxn"]:
            model.reactions.remove(rxnId)

    def get_final_fluxes(self, model_info, model, check_model, model_control_info, initial_rxn_id):
        """"""
        atpm_id = get_atpm_rxn(model_info, model, 'ATP', 'ADP')
        set_model_objective(model_info, model, atpm_id)
        final_atp = round(model.slim_optimize(),3)
        # write_flux_file(model_info, model, atpm_id,self.cfg)
        if "ADD_ATPM" in model.reactions:
            model.reactions.remove('ADD_ATPM')
        nadh_id = add_nadh_rxn(model)
        set_model_objective(model_info, model, nadh_id)
        final_nadh = round(model.slim_optimize(),3)
        # write_flux_file(model_info, model, nadh_id,self.cfg)
        model.reactions.remove('ADD_NADH')
        final_netld_flux = get_final_net_fluxes(model_info, model, check_model)
        if is_autotrophic(model_info):
            final_yield_flux = ''
        else:
            final_yield_flux = get_final_yield_fluxes(model_info, model)
        model.objective = initial_rxn_id
        final_biomass_flux = round(model.slim_optimize(),3)
        print('final_biomass_flux:',final_biomass_flux)
        print(model_control_info["final_flux"])
        model_control_info["final_flux"] = {'final_biomass':final_biomass_flux,
                                            'final_nadh':final_nadh,
                                            'final_atp':final_atp,          
                                            'final_netld':final_netld_flux,
                                            'final_yield':final_yield_flux}



    def biomass_control(self, model_info, model, check_model, model_control_info, obj, model_check_info):
        """"""
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':biomass start')
        model.solver='cplex'
        model_control_info['check_biomass_production'] = {}
        model_control_info['check_biomass_production']["score"] = 1
        model_control_info['check_biomass_production']["Biomass_information"] = []
        # model_control_info["check_biomass_production"]["rely_biomass_rxn"] = ""
        model_control_info["check_biomass_production"]["Check_synthesizability_of_biomass_components"] = []
        model_control_info["check_biomass_production"]["Gapfilling_of_biomass_components"] = []
        model_control_info["check_biomass_production"]["Check_small_molecules"] = []
        model_control_info["check_biomass_production"]["Check_macromolecules"] = []
        model_control_info["check_biomass_production"]["Normalization_of_biomass"] = []
        model_control_info["check_biomass_production"]["Final_biomass"] = ""
        model_control_info["check_biomass_production"]["model_revision"] = []
        # model_control_info["check_biomass_production"]["model_calculation_conditions:"] = []
        model_info['biomass'], model_info["bio_modified"], model_info["no_gap_rxn_input"] = {}, [], []
        set_c_source_supply(model_info, model, 'bio', check_model)
        set_auto_source_supply(model_info, model, check_model, 'bio')
        if not is_bio_exp_correct(model_info):
            model_control_info["check_biomass_production"]["Biomass_information"].append("Could not find the biomass equation, please check and try again!!!")
            model_control_info["check_biomass_production"]["score"] = 0
            return 0
        initial_rxn_id = model_info["initial_rxn"]["biomass_rxn_id"]
        model.reactions.get_by_id(initial_rxn_id).bounds = (0,1000)
        self.get_biomass_info(model_info, model_control_info)
        self.recover_bio_objective(model_info, model)
        self.add_sink_rxn(model, model_info, check_model)
        # self.test(model, check_model)
        # self.add_rely_biomass_rxn(model_info, model, check_model, model_control_info, initial_rxn_id)
        general_library = self.get_general_library()    
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':总库已获取,检查biomass前体物质')
        no_synthesized_mets = self.check_bio_component_is_zero(model_info, model, model_control_info, initial_rxn_id)
        if self.preprocess_initial_zero_bio(model_info, model, model_control_info, general_library, no_synthesized_mets, obj, check_model) == -1:
            # self.convert_list_to_string(model_control_info)
            model_control_info["check_biomass_production"]["Final_biomass"] = f"Initial growth rate: {model_info['biomass_original']}; after setting the carbon source, the value becomes: {model_info['biomass_by_c_source']}; after setting special boundaries, the value becomes: {model_info['biomass_by_special_boundary']}; after quality control, final biomass growth rate: {round(model.slim_optimize(), 3)}.\nThe equations modified during biomass normalization correction are placed in the 'Model Revision' list, and you can make changes and edits by clicking the 'Edit' button." 
            model_info["control_analysis"].append((model_info['biomass_original'],round(model.slim_optimize(), 3)))
            self.recover_atpm_bounds(model_info, model, check_model)
            self.modify_bios(model_info, model, model_control_info)
            self.recover_special_boundary_rxn(model_info, model, check_model)
            return 0
        # write_flux_file(model_info, model, initial_rxn_id, 'biomasss', self.cfg)
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':检查biomass维持能')
        sustain_result = self.is_missing_sustaining_energy(model_info, model, model_control_info, initial_rxn_id)
        macromolecule = self.get_macromolecule(model, initial_rxn_id, model_control_info)
        smallmolecule = self.get_smallmolecule(model_info, model, initial_rxn_id, macromolecule, model_control_info)

        self.get_initial_biomass(model_info, model, initial_rxn_id, model_control_info, check_model, obj)
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':small_fix start')
        self.small_fix(model_info, model, initial_rxn_id, smallmolecule, general_library, model_control_info, check_model, obj)
        # model.solver='cplex'
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':macro_fix start')
        self.macro_fix(model_info, model, initial_rxn_id, macromolecule, general_library, model_control_info, check_model, obj)
        print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),':biomass_fix start')
        if self.biomass_fix(model_info, model, initial_rxn_id, model_control_info, check_model, obj) == -1:
            model_control_info["check_biomass_production"]["Final_biomass"] = f"Initial growth rate: {model_info['biomass_original']}; after setting the carbon source, the value becomes: {model_info['biomass_by_c_source']}; after setting special boundaries, the value becomes: {model_info['biomass_by_special_boundary']}; after quality control, final biomass growth rate: {round(model.slim_optimize(), 3)}.\nThe equations modified during biomass normalization correction are placed in the 'Model Revision' list, and you can make changes and edits by clicking the 'Edit' button."
            model_info["control_analysis"].append((model_info['biomass_original'],round(model.slim_optimize(), 3)))
            self.recover_atpm_bounds(model_info, model, check_model)
            self.modify_bios(model_info, model, model_control_info)
            self.recover_special_boundary_rxn(model_info, model, check_model)
            return 0
        # write_flux_file(model_info, model, initial_rxn_id, 'biomasss', self.cfg)
        self.recover_atpm_bounds(model_info, model, check_model)
        self.get_result(model_info, model, sustain_result, model_control_info)
        self.modify_bios(model_info, model, model_control_info)
        # self.get_exchange(model, model_control_info)
        # self.get_final_fluxes(model_info, model, check_model, model_control_info, initial_rxn_id)
        # self.convert_list_to_string(model_control_info)
        self.recover_special_boundary_rxn(model_info, model, check_model)
        # self.convert_nan_to_null(model)
