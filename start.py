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

parser = argparse.ArgumentParser()

parser.add_argument('--file', type=str, default='', help='model file directory')
parser.add_argument('-o','--outputdir', type=str, default='./', help='result file directory')
parser.add_argument('--types', type=str, default='', help='model control type')
parser.add_argument('--rxns_job', type=list, default='', help='List of user-modified reactions')

args = parser.parse_args()


FILE = Path(__file__).resolve()
ROOT = FILE.parents[0]
# print('FILE:',FILE,'ROOT:', ROOT)
if str(ROOT) not in sys.path:
    sys.path.append(str(ROOT)) # add ROOT to PATH


class ModelControl():
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
        self.model_check_info = {}
        self.model_check_info['boundary_information'] = {}
        self.model_check_info["check_reducing_equivalents_production"] = {}
        self.model_check_info["check_energy_production"] = {}
        self.model_check_info["check_metabolite_production"] = {}    
        self.model_check_info["check_metabolite_yield"] = {}
        self.model_check_info["check_biomass_production"] = {}
        self.model_control_info = {}
        self.model_control_info['boundary_information'] = {}
        self.model_control_info["check_reducing_equivalents_production"] = {}
        self.model_control_info["check_energy_production"] = {}
        self.model_control_info["check_metabolite_production"] = {}    
        self.model_control_info["check_metabolite_yield"] = {}
        self.model_control_info["check_biomass_production"] = {}
        self.model_control_info["quantitative_comparison_before_and_after_correction"] = {}

    def create_outputdir(self,output_dir):
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        return output_dir
    
    

    def model_control(self):
        """
        The overall process of model quality control.

        """
        t1 = time.time()
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

        initial_pre_process = InitialPreprocess(self.cfg)
        initial_pre_process.initial_control(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info, self.model_check_info) 

        quantitative = Quantitative(self.cfg)
        quantitative.get_initial(model_pre_process.model_info, controler.check_model, self.model_control_info)
            # model_control_infos = write_result(self.model_control_info,self.cfg)
            # final_model = write_final_model(controler.model,self.cfg)
            # write_model_info(model_pre_process.model_info,self.cfg)
            # t2 = time.time()
            # print("Total time: ", t2-t1)
            # return model_control_infos, final_model
        nadh_t1 =  time.time()
        nadhs = Nadhs(self.cfg)
        nadhs.nadh_control(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info, self.model_check_info)
        nadh_t2 =  time.time()
        atp_t1 =  time.time()
        atps = Atps(self.cfg)
        atps.atp_control(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info, self.model_check_info) 
        # write_result(self.model_control_info,self.cfg)
        write_model_info(model_pre_process.model_info,self.cfg)
        atp_t2 =  time.time()
        net_t1 =  time.time()
        nets = Nets(self.cfg)
        nets.net_control(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info, self.model_check_info)
        # write_result(self.model_control_info,self.cfg)
        # write_model_info(model_pre_process.model_info,self.cfg)
        net_t2 =  time.time()
        yield_t1 =  time.time()
        yields = Yields(self.cfg)
        yields.yield_control(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info, self.model_check_info) 
        # write_result(self.model_control_info,self.cfg)
        # write_model_info(model_pre_process.model_info,self.cfg)
        yield_t2 =  time.time()
        bio_t1 =  time.time()
        biomasses = Biomasses(self.cfg)
        biomasses.biomass_control(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info, controler, self.model_check_info) 
        bio_t2 =  time.time()

        # with open("result.json", "w") as f:
        #     json.dump(self.model_control_info, f, ensure_ascii=False)
        convert_nan_to_null(controler.model)
        get_final_fluxes(model_pre_process.model_info, controler.model, controler.check_model, self.model_control_info)
        comparison_quantitative(self.model_control_info, model_pre_process.model_info)
        convert_list_to_string(self.model_control_info, self.model_check_info)
        model_control_infos = write_result(self.model_control_info,self.cfg)
        write_model_info(model_pre_process.model_info,self.cfg)
        final_model = write_final_model(controler.model,self.cfg)
        t2 = time.time()
        print("nadh time: ", nadh_t2-nadh_t1, "atp time: ", atp_t2-atp_t1, "net time: ", net_t2-net_t1, "yield time: ", yield_t2-yield_t1, "bio time: ", bio_t2-bio_t1)
        print("Total time: ", t2-t1)
        return model_control_infos, final_model
    
def main(args):

    modelControl = ModelControl(args.file,args.outputdir)
    modelControl.model_control()



if __name__ == '__main__':
    main(args)