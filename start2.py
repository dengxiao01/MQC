#-*- coding:utf-8 -*-
'''
Author: wangruoyu, wangry@tib.cas.cn
Date: 2023-04-20 06:01:56
LastEditors: wangruoyu
LastEditTime: 2023-04-20 06:49:23
Description: file content
FilePath: /mqc/start.py
'''
import argparse
import os, sys
from pathlib import Path
import time 
import csv 
from multiprocessing.pool import Pool

from mqc.config import Config
from mqc.utils import *
from mqc.control.preprocessing_control import Preprocess
from mqc.control.model_control import ModelPreprocess
from mqc.control.initial_control import InitialPreprocess
from mqc.control.nadh_control import Nadhs
from mqc.control.atp_control import Atps
from mqc.control.net_control import Nets
from mqc.control.yield_control import Yields
from mqc.control.biomass_control import Biomasses
from mqc.control.quantitative_control import Quantitative
from mqc.control.yield_control2 import Yields2
from mqc.control.check_control import Check
from mqc.control.rules_control import Rules
# parser = argparse.ArgumentParser()

# parser.add_argument('--file', type=str, default='', help='model file directory')
# parser.add_argument('-o','--outputdir', type=str, default='./', help='result file directory')
# parser.add_argument('--types', type=str, default='', help='model control type')
# parser.add_argument('--rxns_job', type=list, default='', help='List of user-modified reactions')

# args = parser.parse_args()


FILE = Path(__file__).resolve()
ROOT = FILE.parents[0]
# print('FILE:',FILE,'ROOT:', ROOT)
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT)) # add ROOT to PATH


class ModelCheck():
    """
    Obtain the total output information of the model quality control.

    """
    def __init__(self, model_file:str,output_dir:str):
        """
        define output dictionary.

        """
        # input model
        self.model_file = model_file
        self.output_dir = self.create_outputdir(output_dir)
        self.cfg = Config(self.output_dir)
        self.model_control_info = {}
        self.model_check_info = {}
        self.model_check_info['boundary_information'] = {}
        self.all_data = []
        self.model_check_info["check_reducing_equivalents_production"] = {}
        self.model_check_info["check_energy_production"] = {}
        self.model_check_info["check_metabolite_production"] = {}    
        self.model_check_info["check_metabolite_yield"] = {}
        self.model_check_info["check_biomass_production"] = {}

        self.model_control_info['boundary_information'] = {}
        self.model_control_info["check_reducing_equivalents_production"] = {}
        self.model_control_info["check_energy_production"] = {}
        self.model_control_info["check_metabolite_production"] = {}    
        self.model_control_info["check_metabolite_yield"] = {}
        self.model_control_info["check_biomass_production"] = {}
        # self.model_control_info["quantitative_comparison_before_and_after_correction"] = {}

    def create_outputdir(self,output_dir):
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        return output_dir
    
    def model_check(self):
        """"""    
        t1 = time.time()
        headers = ['model', 'reducing_power', 'energy', 'metabolite', 'yield', 'biomass','all']
        controler = Preprocess(self.model_file,self.cfg)
        if not controler.model:
            model_control_info = change_model_control_info(self.model_file, self.model_control_info)
            model_control_infos = write_result2(model_control_info,self.cfg)
            t2 = time.time()
            print("Total time: ", t2-t1)
            return model_control_infos
        model_pre_process = ModelPreprocess(self.cfg)
        model_pre_process.get_model_info(controler)
        checks = Check(self.cfg)
        all_data = checks.check_control(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info, self.model_check_info)
        print(all_data)
        # model_name = f'{controler.model}'
        # if not f'{controler.model}':
        #     model_name = self.model_file.split('/')[-1].split('.')[0]
        # if '0' in all_data:
        #     all_data.append('0')
        # else:
        #     all_data.append('1')
        # # 创建一个数据框（DataFrame）
        # df = pd.DataFrame([all_data], columns=headers)
        # # 尝试读取现有文件
        # try:
        #     existing_df = pd.read_excel(f'/home/dengxiao/mqc/tmp/GCF_modelseed_check/{model_name}.xlsx', sheet_name='Sheet1')
        #     # 将新的数据追加到现有数据框
        #     existing_df = existing_df.append(df, ignore_index=True)
        # except FileNotFoundError:
        #     # 如果文件不存在，直接使用新数据框
        #     existing_df = df
        # # 将数据框写入Excel文件
        # existing_df.to_excel(f'/home/dengxiao/mqc/tmp/GCF_modelseed_check/{model_name}.xlsx', index=False, header=True, sheet_name='Sheet1')
        del self.model_check_info['boundary_information']
        model_control_infos = write_result(self.model_check_info,self.cfg)
        t2 = time.time()
        print("Total time: ", t2-t1)
        return model_control_infos, all_data

    def model_check2(self):
        """
        The overall process of model quality control.

        """
        t1 = time.time()
        model_control_infos, final_model = '', ''
        self.model_check_info['boundary_information'] = {}
        controler = Preprocess(self.model_file,self.cfg)
        if not controler.model:
            model_control_info = change_model_control_info(self.model_file, self.model_control_info)
            model_control_infos = write_result2(model_control_info,self.cfg)
            final_model = self.model_file
            t2 = time.time()
            print("Total time: ", t2-t1)
            return model_control_infos
        model_pre_process = ModelPreprocess(self.cfg)
        model_pre_process.get_model_info(controler)
        model_pre_process.model_info["model_file"] = self.model_file
        try:
            initial_pre_process = InitialPreprocess(self.cfg)
            initial_pre_process.initial_control(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info, self.model_check_info) 
        # quantitative = Quantitative(self.cfg)
        # quantitative.get_initial(model_pre_process.model_info, controler.check_model, self.model_control_info)
            # model_control_infos = write_result(self.model_control_info,self.cfg)
            # final_model = write_final_model(controler.model,self.cfg)
            # write_model_info(model_pre_process.model_info,self.cfg)
            # t2 = time.time()
            # print("Total time: ", t2-t1)
            # return model_control_infos, final_model
        # try:
            nadhs = Nadhs(self.cfg)
            nadhs.nadh_control(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info, self.model_check_info, controler)
            # write_model_info(model_pre_process.model_info,self.cfg)
            atps = Atps(self.cfg)
            atps.atp_control(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info, self.model_check_info, controler) 
            # write_result3(self.model_control_info,self.cfg)
            # write_model_info(model_pre_process.model_info,self.cfg)
            nets = Nets(self.cfg)
            nets.net_control(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info, self.model_check_info, controler)
            # write_result(self.model_control_info,self.cfg)
            # write_model_info(model_pre_process.model_info,self.cfg)
            yields = Yields(self.cfg)  
            yield_one = yields.yield_control(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info, controler, self.model_check_info) 
            if yield_one == 1:
                return_restoration(model_pre_process.model_info, controler.model, controler.check_model, Biomasses(self.cfg))
                yield_two = yields.yield_control(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info, controler, self.model_check_info) 
                # self.model_control_info["check_metabolite_yield"]["model_revision"].extend(model_pre_process.model_info["yield_revision"])
            # write_result(self.model_control_info,self.cfg)
            # write_model_info(model_pre_process.model_info,self.cfg)
            # final_model = write_final_model(controler.model,self.cfg)
            # if final_model.endswith(".json"):
            #     controler.model = cobra.io.load_json_model(final_model)
            # else:
            #     controler.model = cobra.io.read_sbml_model(final_model)
            biomasses = Biomasses(self.cfg)
            biomasses.biomass_control(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info, controler, self.model_check_info) 
            # write_model_info(model_pre_process.model_info,self.cfg)
            # write_result(self.model_control_info,self.cfg)
            # with open("result.json", "w") as f:
            #     json.dump(self.model_control_info, f, ensure_ascii=False)
            convert_nan_to_null(controler.model)
            convert_nan_to_null(controler.check_model)
            # get_final_fluxes(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info)
            # comparison_quantitative(self.model_control_info, model_pre_process.model_info)
            convert_list_to_string(self.model_control_info, self.model_check_info)
            # print(controler.model.reactions.get_by_id('FDR').check_mass_balance(),'xxxxxxxxxxxxxx')
            # final_model = f"/home/dengxiao/mqc/tmp/bigg/test.xml"
            # cobra.io.write_sbml_model(controler.model,final_model)
            final_model = write_final_model(controler.model,self.cfg, self.model_file)
            check_model = write_check_model(controler.check_model,self.cfg, self.model_file)   
        except RuntimeError as e:
            final_model = write_final_model(controler.model,self.cfg, self.model_file)
            check_model = write_check_model(controler.check_model,self.cfg, self.model_file) 
            model_control_infos = write_result2(self.model_control_info,self.cfg)
            print(repr(e),'.............')
            # raise
        except Exception as e:
            model_control_infos = write_result2(self.model_control_info,self.cfg)
            final_model = write_final_model(controler.model,self.cfg, self.model_file)
            check_model = write_check_model(controler.check_model,self.cfg, self.model_file) 
            print(repr(e),'.............')
            raise 
        print('0000000000000000000000000000000000000000000000000000000000')
        control_analysis = model_pre_process.model_info["control_analysis"]
        headers=['model','NADH','NADPH','FADH2','FMNH2','Q8H2','MQL8','DMMQL8','ATP','CTP','GTP','UTP','ITP','metabolite','yield','biomass','carbon source supply','restricted metabolites']
        if len(control_analysis) < 16:
            for i in range(16-len(control_analysis)):
                control_analysis.append("")
        control_analysis.append(model_pre_process.model_info["carbon_source_boundary"])
        control_analysis.append(model_pre_process.model_info["limiting_metabolites"])
        model_name = f'{controler.model}'
        if not f'{controler.model}':
            model_name = self.model_file.split('/')[-1].split('.')[0]
        # 创建一个数据框（DataFrame）
        df = pd.DataFrame([control_analysis],columns=headers)
        # 尝试读取现有文件
        try:
            existing_df = pd.read_excel(f'/home/dengxiao/mqc/tmp/GCF_modelseed/{model_name}.xlsx', sheet_name='Sheet1')
            # 将新的数据追加到现有数据框
            existing_df = existing_df.append(df, ignore_index=True)
        except FileNotFoundError:
            # 如果文件不存在，直接使用新数据框
            existing_df = df
        # 将数据框写入Excel文件
        existing_df.to_excel(f'/home/dengxiao/mqc/tmp/GCF_modelseed/{model_name}.xlsx', index=False, header=True, sheet_name='Sheet1')
        modelInfo = write_model_info(model_pre_process.model_info,self.cfg)
        # model_check_infos = write_result(self.model_check_info,self.cfg)
        # if not model_control_infos:
        model_control_infos = write_result3(self.model_control_info,self.cfg)
        t2 = time.time()
        print(controler.model.slim_optimize(),'...............................................................')
        print("Total time: ", t2-t1)
        return model_control_infos, final_model, check_model


class ModelControl():
    """
    Obtain the total output information of the model quality control.

    """
    def __init__(self, model_file:str, output_dir:str, types:str, rxns_job:list, check_model:str):
        """
        define output dictionary.

        """
        # input model
        self.output_dir = self.create_outputdir(output_dir)
        self.cfg = Config(self.output_dir)
        self.model_file = model_file
        self.types = types
        self.rxns_job = rxns_job
        self.check_model = check_model
        self.model_check_info = {}
        self.model_check_info['boundary_information'] = {}
        self.model_check_info["check_reducing_equivalents_production"] = {}
        self.model_check_info["check_energy_production"] = {}
        self.model_check_info["check_metabolite_production"] = {}    
        self.model_check_info["check_metabolite_yield"] = {}
        self.model_check_info["check_biomass_production"] = {}
        # self.model_check_info["quantitative_comparison_before_and_after_correction"] = {}
        self.model_control_info = {}
        self.model_control_info['boundary_information'] = {}
        self.model_control_info["check_reducing_equivalents_production"] = {}
        self.model_control_info["check_energy_production"] = {}
        self.model_control_info["check_metabolite_production"] = {}    
        self.model_control_info["check_metabolite_yield"] = {}
        self.model_control_info["check_biomass_production"] = {}
        # self.model_control_info["quantitative_comparison_before_and_after_correction"] = {}
        # self.model_control_info['initial_flux'] = {}
        self.temp_model_control_info = {}
        self.model_control_info["check_biomass_production"]["Check_synthesizability_of_biomass_components"] = []
        self.model_control_info["check_biomass_production"]["Gapfilling_of_biomass_components"] = []

    def create_outputdir(self,output_dir):
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        return output_dir
    
    def get_initial_check_model(self, controler, check_model):
        """"""
        if check_model.endswith(".json"):
            controler.check_model = cobra.io.load_json_model(check_model)
        elif check_model.endswith(".yaml"):
            controler.check_model = cobra.io.load_yaml_model(check_model)
        elif check_model.endswith(".mat"):
            controler.check_model = cobra.io.load_matlab_model(check_model)
        else:
            controler.check_model = cobra.io.read_sbml_model(check_model)
    

    def model_control(self):
        """
        The overall process of model quality control.

        """
        t1 = time.time()

        model_path = get_new_model(self.model_file, self.rxns_job, self.cfg)
        model_info = get_model_info(self.cfg)
        reduct_modelinfo_modify(model_info)
        controler = Preprocess(model_path,self.cfg)
        self.get_initial_check_model(controler,self.check_model)
        # controler.model.reactions.get_by_id("MTHFC_1").bounds = (0,0)
        print(controler.model.slim_optimize(),'...............................................................')
        try:
            initial_pre_process = InitialPreprocess(self.cfg)
            rules = Rules(self.cfg) 
            run_rules(rules, model_info, controler.model)
            # initial_pre_process.initial_control(model_info, controler.model, controler.check_model, self.model_control_info, self.model_check_info) 
                # model_control_infos = write_result(self.model_control_info,self.cfg)
                # final_model = write_final_model(controler.model,self.cfg)
                # write_model_info(model_pre_process.model_info,self.cfg)
                # t2 = time.time()
                # print("Total time: ", t2-t1)
                # return model_control_infos, final_model
            if self.types == "boundary_information":
                initial_pre_process.initial_control(model_info, controler.model, controler.check_model, self.model_control_info, self.model_check_info) 
                # self.model_control_info['boundary_information']['exchange_reaction_boundary'] = initial_pre_process.get_exchange(controler.model, model_info, controler.check_model)

            if self.types == "check_reducing_equivalents_production":
                nadhs = Nadhs(self.cfg)
                nadhs.nadh_control(model_info, controler.model, controler.check_model, self.model_control_info, self.model_check_info, controler)
        
            if self.types == "check_energy_production":
                atps = Atps(self.cfg)
                atps.atp_control(model_info, controler.model, controler.check_model, self.model_control_info, self.model_check_info, controler) 
                # write_result(self.model_control_info,self.cfg)
                write_model_info(model_info,self.cfg)
        
            if self.types == "check_metabolite_production":
                nets = Nets(self.cfg)
                nets.net_control(model_info, controler.model, controler.check_model, self.model_control_info, self.model_check_info, controler)
                # write_result(self.model_control_info,self.cfg)
                # write_model_info(model_pre_process.model_info,self.cfg)
            if self.types == "check_metabolite_yield":
                yields = Yields(self.cfg)
                yields.yield_control(model_info, controler.model, controler.check_model, self.model_control_info, controler, self.model_check_info) 
                # write_result(self.model_control_info,self.cfg)
                # write_model_info(model_pre_process.model_info,self.cfg)
            if self.types == "check_biomass_production":
                biomasses = Biomasses(self.cfg)
                biomasses.biomass_control(model_info, controler.model, controler.check_model, self.model_control_info, controler, self.model_check_info) 
                convert_list_to_string(self.model_control_info, self.model_check_info)
            # if self.types == "quantitative_comparison_before_and_after_correction":
            #     quantitative = Quantitative(self.cfg)
            #     quantitative.get_initial(model_info, controler.check_model, self.model_control_info)
            #     get_final_fluxes(model_info, controler.model, controler.check_model, self.model_control_info)
            #     comparison_quantitative(self.model_control_info, model_info)
            # with open("result.json", "w") as f:
            #     json.dump(self.model_control_info, f, ensure_ascii=False)
            convert_nan_to_null(controler.model)
            final_model = write_final_model(controler.model,self.cfg, self.model_file)
        except RuntimeError as e:
            final_model = ""
            print(repr(e),'.............')
        except Exception as e:
            print(repr(e),'.............')
            raise 
        self.temp_model_control_info[self.types] = self.model_control_info[self.types]
        model_control_infos = write_result4(self.temp_model_control_info,self.cfg)
        modelInfo =write_model_info(model_info,self.cfg)
        t2 = time.time()
        print("Total time: ", t2-t1)
        return model_control_infos, final_model, modelInfo

def main(a,b):
    """"""
    # modelCheck = ModelCheck(args.file,args.outputdir)
    # a=f'/home/dengxiao/mqc/mqc/local_test_data/CARVEME_COMMEN/{file}'
    # b=f"tmp/new_CARVEME_COMMEN/{file.split('.')[0]}"
    modelCheck = ModelCheck(a,b)
    # modelCheck.model_check()
    modelCheck.model_check2()
# def main(file):
#     """"""
#     # a=f"mqc/local_test_data/literature_model/{file}"
#     # b=f"./tmp/literature2/{file.split('.')[0]}"
#     a=f"mqc/local_test_data/bigg_data/{file}"
#     b=f"tmp/bigg/IJN746"
#     print(a,'.....................')
#     modelCheck = ModelCheck(a,b)
#     result = modelCheck.model_check()
#     if type(result) == tuple:
#         result[1].append(f"{file.split('.')[0]}")
#     return result 
    
    

def main2(args):
    """"""
    args.file = "/home/dengxiao/mqc/tmp/bigg/iCN718/iCN718.xml"
    args.outputdir = "/home/dengxiao/mqc/tmp/bigg/icn718"
    args.types = "check_reducing_equivalents_production"
    args.check_model = '/home/dengxiao/mqc/tmp/bigg/icn7182/check_model.json'
    modelControl = ModelControl(args.file,args.outputdir,args.types,args.rxns_job,args.check_model)
    modelControl.model_control()

if __name__ == '__main__':
    """"""
    parser = argparse.ArgumentParser()
    # all_datas = []
    # all_datas.append(['model','reducing_power','energy','metabolite','yield','biomass','outputdir'])
    # a,b=[],[]
    parser.add_argument('--file', type=str, default='', help='model file directory')
    parser.add_argument('-o','--outputdir', type=str, default='./', help='result file directory')
    parser.add_argument('--types', type=str, default='', help='model control type')
    parser.add_argument('--rxns_job', type=list, default='', help='List of user-modified reactions')
    parser.add_argument('--check_model', type=str, default='', help='check_model file directory')

    # # files = os.listdir('mqc/local_test_data/literature_model')
    # # files.remove('human2.0.xml')
    # # files.remove('iHepatocytes2322.xml')
    # # with Pool() as pool:
    # #     result = pool.map(main,files)
    # # for file in files:
    # #     result = main(file)
    # file = 'iJN746.xml'
    # result = main(file)
    # print('................................................................................')
    # print(result)
    # # print('................................................................................')
    # # print(type(result))
    # # if type(result) == tuple:
    # #     print(result[1],'111111111111111')
    # #     all_datas.append(result[1])
    # # print(all_datas)
    # # with open('mqc/local_test_data/literature_model/summarize.csv','w',newline='') as f:
    # #     writer = csv.writer(f)
    # #     for row in all_datas:
    # #         writer.writerows(row)

    args = parser.parse_args()
    # main(args)

    # main2(args)
    alreadyrun = ['iAT601.xml','iNX1344.xml','iPSA590.xml','iSyf7152.xml','iYL619_PCP.xml','iYLI647.xml','PstM1.xml','iCAC490.xml','iEM439.xml','iHZ771.xml','iJN1462.xml',
    'human2.0.xml','iHepatocytes2322.xml','iAK888.xml','NmrFL413.xml','iBP722.xml','iPN730.xml','Yeast5.xml'] 
    # 'iBP722.xml'求解器报错，但单碳的时候没有报错 'Yeast5.xml':Segmentation fault (core dumped)内存只想错误
    # iCTH669.xml写模型报错    iAK888.xml：净物质循环  iHepatocytes2322.xml:check_mass_balance()报错
    # files = os.listdir('mqc/local_test_data/CARVEME_COMMEN')
    # yes_model = ['Aspm1282.xml','BPK282A1.xml','iATCC19606.xml','iBSCc.xml','iCbu641.xml','iCN718.xml','iBP722.xml']
    yes_model=['Aspm1282.xml', 'iATCC19606.xml', 'BPK282A1.xml', 'iBSCc.xml','iCbu641.xml', 'iCN718.xml','iES1300.xml', 'iCJ415.xml', 'iCth446.xml', 'iGD1575.xml', 
    'iHL622_2.xml', 'iGEL604.xml', 'iHY3410.xml', 'ihGlycopastoris.xml', 'iIsor850.xml', 'iHL622.xml', 'iJB1325_CBS513.88.xml', 'iHM1533.xml', 'iJB1325_ATCC1015.xml', 
    'iKS1317.xml', 'iLM.c559.xml', 'iME375.xml', 'iLME620.xml', 'iMSC1255.xml', 'iLT1021.xml', 'iMcBath.xml', 'iMM1865.xml', 'iMT1026.xml', 'iMsOB3b.xml', 'iPC815.xml',
     'iNS934.xml', 'iPS584.xml', 'iRsp1140.xml', 'iPB890.xml', 'iSyu683.xml','MetaMerge.xml', 'MTBPROM2.0.xml', 'iYY11011.xml','iCTH669.xml']
    run_file=[]
    # table = pd.read_csv('/home/dengxiao/mqc/mqc/local_test_data/embl_gems/model_list.tsv',sep='\t')
    # files = os.listdir('mqc/local_test_data/CARVEME_COMMEN')
    # for file in table['file_path'][2416:2500]:
    #     # if file in alreadyrun or file in yes_model:
    #     #     continue
    #     # file = f"/home/dengxiao/mqc/mqc/local_test_data/embl_gems/{files}"
    #     run_file.append(file)
    #     # print(run_file)
    #     main(f'/home/dengxiao/mqc/mqc/local_test_data/embl_gems/{file}',f"tmp/web_CARVEME_COMMEN/{file.split('/')[-1].split('.')[0]}")
    modelseed_files = os.listdir('/home/dengxiao/mqc/mqc/local_test_data/bigg_data')
    # for file in modelseed_files[100:]:
    #     main(f'/home/dengxiao/mqc/mqc/local_test_data/bigg_data/{file}',f"tmp/BiggAll/{file.split('.xml')[0]}") 
    # main('/home/dengxiao/mqc/mqc/local_test_data/seedpy_comon_model/Streptococcus_pneumoniae_seedpy_model.xml',"tmp/new_seedpy_comon_model/Streptococcus_pneumoniae_seedpy_model") 
    # main('/home/dengxiao/mqc/mqc/local_test_data/metanetx2/seed_Seed107806_1.COBRA-sbml3.xml',"tmp/metanetx/seed_Seed107806_1") 
    # main('/home/dengxiao/mqc/mqc/local_test_data/GCF_modelseed/GCA_000007365.1modelseed.xml',"tmp/GCF_modelseed/GCA_000007365") 
    # main('/home/dengxiao/mqc/mqc/local_test_data/embl_gems/models/f/facklamia/Facklamia_hominis_CCUG_36813.xml.gz',"tmp/CARVEME_COMMEN/Facklamia_hominis_CCUG_36813") 
    main('/home/dengxiao/mqc/mqc/local_test_data/other/iCac802_norm2.xml',"tmp/other/iCac802_norm2")
    # main('/home/dengxiao/mqc/mqc/local_test_data/other/methanolica_excel_xml.xml','tmp/other/methanolica_excel_xml')

    