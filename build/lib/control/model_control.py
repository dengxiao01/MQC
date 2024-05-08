# -*- coding: utf-8 -*-

from cobra.util.solver import linear_reaction_coefficients
import json 

from mqc.config import Config
from mqc.utils import *


class ModelPreprocess():
    """
    Integrate the model's attributes into a dictionary.

    Attributes
    ----------
    diff_results : python.Dictionary
        The dictionary structure of memote.MemoteResult objects.
    configuration : memote.MemoteConfiguration
        A memote configuration structure.

    """
    def __init__(self,cfg:Config) -> None:
        """
        Initialize the model dictionary.

        """
        self.cfg = cfg
        self.model_info = {}
        
    
    def old_preprocess_met(self, obj):
        """
        Add preprocessing operations such as formula and charge to metabolites.

        Notes
        -----
        The name of the seed library needs to be mapped to the table, which is mapped by modelseed_id_to_name

        """
        for met in obj.model.metabolites:
            if met.id.startswith("cpd"):  # seed
                obj.add_seed_formula()
                self.modelseed_id_to_name = obj.bind_id_name()
                break 
            elif len(met.id.split("_")[-1]) == 1:  # bigg
                obj.add_bigg_formula()
                break
            elif "@" in met.id:  # meta
                obj.add_meta_formula()
                break
            elif "[" in met.id and "]" in met.id:  # virtual
                break
            elif met.id.startswith("C"):  # kegg
                obj.add_kegg_formula()
                break 
            else:
                obj.add_metacyc_formula()
                break
                """其他库--kegg,metacyc...bigg不好区分也可放进来,做成一个大表,这个表只包含formula就行"""

    def remove_not_in_rxn_met(self, obj):
        """"""
        # 创建一个列表来存储待移除的代谢物
        metabolites_to_remove = []
        # 遍历模型中的所有代谢物
        for metabolite in obj.model.metabolites:
            # 检查代谢物是否参与了任何反应
            if len(metabolite.reactions) == 0:
                metabolites_to_remove.append(metabolite)
        # 从模型中移除待移除的代谢物
        obj.model.remove_metabolites(metabolites_to_remove)
        # 打印移除了多少个代谢物
        print(f'Removed {len(metabolites_to_remove)} metabolites..........')

    def remove_no_met_rxn(self, obj):
        """"""
        rxns_to_remove = []
        for rxn in obj.model.reactions:
            if not list(rxn.metabolites.keys()):
                rxns_to_remove.append(rxn)
        obj.model.remove_reactions(rxns_to_remove)
        print(f'Removed {len(rxns_to_remove)} rxns..........')

    def preprocess_met(self, obj):
        """
        Add preprocessing operations such as formula and charge to metabolites.

        Notes
        -----
        The name of the seed library needs to be mapped to the table, which is mapped by modelseed_id_to_name

        """
        seed_met, meta_met, virtual_met, kegg_met, other_met = [], [], [], [], []
        pattern = r'^C\d{5}$' 
        for met in obj.model.metabolites:
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
        if len(virtual_met) != 0:
            obj.add_virtual_formula(virtual_met)
        if len(kegg_met) != 0:
            obj.add_kegg_formula(kegg_met)
        if len(other_met) != 0:
            obj.add_other_formula(other_met)
        

    def modify_met_formula(self, obj):
        """"""
        for met in obj.model.metabolites:
            if pd.isna(met.formula) or met.formula == "nan":
                met.formula = ""
            if pd.isna(met.charge) or met.charge == "nan":
                met.charge = 0
            # if '*' in met.formula: # C26H36O24*2
            #     met.formula = met.formula.split('*')[0]
            # if '(' in met.formula:
            #     met.formula = met.formula.split('(')[0]
                
    def modify_met_name(self, obj):
        """"""
        flag = 0
        for met in obj.model.metabolites:
            if pd.isna(met.name):
                met.name = ''
        for met in obj.model.metabolites:
            if '_' in str(met.name) and len(str(met.name).split('_')[-1]) == 2:
                flag += 1
        if flag == len(obj.model.metabolites):
            for met in obj.model.metabolites:
                met.name = met.name.split('_')[0]

    def get_met_name(self, obj, metId):
        met = obj.model.metabolites.get_by_id(metId)
        return met.name  # IDs that do not belong to the four libraries


    def preprocess_rxn(self, obj):  
        """
        Additive exchange reaction without exchange reaction.

        Notes
        -----
        Find the intersection of all reaction IDs with the exchange reaction list

        """
        met_name = []
        for rxn in obj.model.reactions:
            all_mets = [m.id for m in rxn.metabolites]
            if len(all_mets)==1:
                met_name.append(self.get_met_name(obj, all_mets[0]))
        if len(set(met_name) & set(H2O_NAME)) == 0:
            obj.add_exchange_rxn("add_h2o_exchange")
        if len(set(met_name) & set(NH4_NAME)) == 0:
            obj.add_exchange_rxn("add_nh4_exchange")
        if len(set(met_name) & set(PI_NAME)) == 0:
            obj.add_exchange_rxn("add_pi_exchange")
        if len(set(met_name) & set(H_NAME)) == 0:
            obj.add_exchange_rxn("add_h_exchange")
        if len(set(met_name) & set(O2_NAME)) == 0:
            obj.add_exchange_rxn("add_o2_exchange")
        if len(set(met_name) & set(CO2_NAME)) == 0:
            obj.add_exchange_rxn("add_co2_exchange")
        if len(set(met_name) & set(SO4_NAME)) == 0:
            obj.add_exchange_rxn("add_so4_exchange")


    def get_identifier(self, obj):
        """"""
        identifier = ""
        for met in obj.model.metabolites:
            if met.id.startswith("cpd"):  # seed
                identifier = "modelseed"
                break 
            elif len(met.id.split("_")[-1]) == 1:  # bigg
                identifier = "bigg"
                break 
            elif "@" in met.id:  # metaNetx
                identifier = "metaNetx"
                break 
            elif "[" in met.id and "]" in met.id:  # virtual
                identifier = "virtual"
                break
            else:
                identifier = "other"
                break 
        return identifier


    def get_met_info(self, obj):
        """
        Collect metabolite information.

        Returns
        -------
        met_info : list
            metabolite information

        """
        met_info = []
        for met in obj.model.metabolites:
            if "[" in met.id and "]" in met.id:  # virtual
                metInfo = {"id" : met.id,
                            "name" : met.id.split('[')[0],
                            "charge" : met.charge,
                            "formula" : met.formula}
            else:  # IDs that do not belong to the four libraries
                metInfo = {"id" : met.id,  
                            "name" : met.name,
                            "charge" : met.charge,
                            "formula" : met.formula}
            met_info.append(metInfo)
        return met_info

    

    def get_rxn_info(self, obj):
        """
        Collect reaction information.

        Returns
        -------
        rxn_info : list
            reaction information

        """
        rxn_info = []
        for rxn in obj.model.reactions:
            rxnInfo = {"id" : rxn.id,
                        "name" : rxn.name,
                        "bounds" : rxn.bounds,
                        # "lower_bound" : rxn.lower_bound,
                        # "upper_bound" : rxn.upper_bound,
                        "rxn_exp_id" : rxn.build_reaction_string(),
                        "rxn_exp_name" : rxn.build_reaction_string(use_metabolite_names = True),
                        "annotation" : rxn.annotation}
            rxnInfo["all_mets"] = [self.get_met_name(obj, met.id) for met in rxn.metabolites]   
            rxnInfo["reactants_mets"] = [self.get_met_name(obj, met.id) for met in rxn.reactants] 
            rxnInfo["products_mets"] = [self.get_met_name(obj, met.id) for met in rxn.products] 
            rxnInfo["rules"] = {}
            rxnInfo["rules2"] = {}
            rxn_info.append(rxnInfo)
        return rxn_info


    def get_exchange_rxns(self):
        """
        Get exchange reaction ID.

        Returns
        -------
        exchange_rxns : list
            exchange reaction id

        """
        exchange_rxns = []
        for rxn in self.model_info["reactions"]:
            if len(rxn["all_mets"]) == 1:
                exchange_rxns.append(rxn["id"])
        return exchange_rxns

    def get_transport_rxns(self):
        """
        Get the transport reaction ID.

        Returns
        -------
        transport_rxns : list
            transport reaction id

        """
        transport_rxns = []
        for rxn in self.model_info["reactions"]:
            if sorted(rxn["reactants_mets"]) == sorted(rxn["products_mets"]):
                transport_rxns.append(rxn["id"])
            if len(rxn["all_mets"]) == 2 and len(set(H_NAME) & set(rxn['all_mets'])) == 2:
                transport_rxns.append(rxn["id"])
        return transport_rxns


    def get_initial_obj_info(self, obj):
        """"""
        biomassIds = []
        try:
            initial_rxn_id = list(linear_reaction_coefficients(obj.model).keys())[0].id
        except:
            initial_rxn_id = ''
        first_biomass_rxn_id = first_get_biomass_rxn(obj.model, initial_rxn_id)
        if first_biomass_rxn_id:
            biomass_rxn_id = first_biomass_rxn_id
            print('first_biomass_rxn_id:', first_biomass_rxn_id,'.......................')
        else:
            all_biomassIds = get_biomass_rxn(obj.model)
            print('all_biomassIds:', all_biomassIds,'.......................')
            if len(all_biomassIds) > 1: # 如果找到了多个biomass方程，需要进行判断
                with obj.model as model:
                    for rxnId in all_biomassIds:
                        model.objective = rxnId
                        if model.slim_optimize() > 1e-6:
                            biomassIds.append(rxnId)
            else:
                biomassIds = all_biomassIds
            # biomassIds=['BIOMASS_glyc']
            print('biomassIds:', biomassIds,'.......................')
            if biomassIds:
                if initial_rxn_id in biomassIds:
                    biomass_rxn_id = initial_rxn_id
                else:
                    biomass_rxn_id = biomassIds[0]
            else:
                # biomass_rxn_id = initial_rxn_id  以前是在最后仍然找不到biomass的时候把目标反应作为biomass，这样可能会使一些非biomass反应作为biomass，所以现在改为直接告诉用户找不到
                biomass_rxn_id = ''

        if biomass_rxn_id:
            not_contained = check_biomass_rxn(obj.model, biomass_rxn_id)
            if not_contained:
                check_bio = f"Substances {', '.join(not_contained)} were missing in the biomass function."
            else:
                check_bio = f"No substances were lacking in the biomass composition equation {biomass_rxn_id}."
            obj.model.objective = biomass_rxn_id
            biomass_rxn_flux = round(obj.model.slim_optimize(),2)
            try:
                # biomass_rxn_exp = obj.model.reactions.get_by_id(biomass_rxn_id).build_reaction_string(use_metabolite_names=True)
                initial_bio_rxn = obj.model.reactions.get_by_id(biomass_rxn_id)
                biomass_rxn_exp, initial_bio_rxn_id = keep_decimal_places(initial_bio_rxn)
            except KeyError:
                biomass_rxn_exp = ''
        else:
            biomass_rxn_flux = 0
            biomass_rxn_exp = ''
            check_bio = ''
        
        InitialRxnInfo = {"model_default_initial_rxn_id" : initial_rxn_id,
                          "biomass_rxn_id" : biomass_rxn_id,
                            "biomass_rxn_flux" : biomass_rxn_flux,
                            "biomass_rxn_exp" : biomass_rxn_exp,
                            "check_biomass_composition" : check_bio}
        # initial_rxn_info.append(InitialRxnInfo)
        return InitialRxnInfo

    def find_glucose(self):
        """"""
        glucose_rxnId = ''
        for rxn in self.model_info["reactions"]:
            for rxnId in GLUCOSE_RXN:
                if rxn['id'] in self.model_info['exchange_rxns'] and rxn['id'] == rxnId:
                    glucose_rxnId = rxnId
                    # print(glucose_rxnId,rxn['id'])
        if not glucose_rxnId:
            for rxn in self.model_info["reactions"]:
                if rxn['id'] in self.model_info['exchange_rxns'] and ('glucose' == rxn['all_mets'][0].lower() or 'd-glucose' == rxn['all_mets'][0].lower()):
                    glucose_rxnId = rxn['id']
        return glucose_rxnId
        


    def get_model_info(self, controler):
        """
        Integrate model properties together.

        Parameters
        ----------
        file_path : str
            model file

        """
        
        self.remove_not_in_rxn_met(controler)
        self.remove_no_met_rxn(controler)
        self.preprocess_met(controler)
        self.modify_met_formula(controler)
        self.modify_met_name(controler)
        self.preprocess_rxn(controler)
        self.model_info = {"model_id" : controler.model.id,
                        #    "model_identifier" : self.get_identifier(controler),
                            "metabolites" : self.get_met_info(controler),
                            "reactions" : self.get_rxn_info(controler)}
        self.model_info["initial_rxn"] = self.get_initial_obj_info(controler)
        self.model_info["all_rxn_id"] = [rxn.id for rxn in controler.model.reactions]
        self.model_info["exchange_rxns"] = self.get_exchange_rxns()
        self.model_info["transport_rxns"] = self.get_transport_rxns()
        self.model_info["glucose_rxnId"] = self.find_glucose()
        self.model_info['need_close_rely_rxn'] = ""
        self.model_info["yield_revision"] = []
        self.model_info["control_analysis"] = []
        modelInfo = json.dumps(self.model_info, ensure_ascii=False, allow_nan = True, indent=1)
        with open(self.cfg.model_info, "w", newline='',) as f:
            f.write(modelInfo)
        
        



