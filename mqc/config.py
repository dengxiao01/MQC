
from pathlib import Path

FILE = Path(__file__).resolve()
# ROOT = ''
ROOT = FILE.parents[0]
# print('FILE2:',FILE,'ROOT2:', ROOT)
class Config(object):
    def __init__(self, output_dir) -> None:
        # result file
        self.modle_file = "mqc/local_test_data/bigg_data/iML1515.xml"
        self.output_dir = output_dir
        self.flux_path = f"{output_dir}/flux.txt"
        self.flux_path2 = Path(output_dir) / "flux2.txt"
        self.flux_path3 = Path(output_dir) / "flux3.txt"
        self.flux_path4 = Path(output_dir) / "flux.tsv"
        self.result_file = Path(output_dir) / "check_result.json"
        self.result_file2 = Path(output_dir) / "fail.txt"
        self.result_file3 = Path(output_dir) / "control_result.json"
        self.result_file4 = Path(output_dir) / "new_result.json"
        # self.big_general_library_xml = Path(output_dir) / "big_model.xml"
        self.big_general_library_xml = f"{output_dir}/big_model.xml"
        # self.big_general_library_json = Path(output_dir) / "big_model.json"
        self.big_general_library_json = f"{output_dir}/big_model.json"
        self.big_general_library_yaml = f"{output_dir}/big_model.yaml"
        self.big_general_library_mat = f"{output_dir}/big_model.mat"
        self.big_general_library = ""
        self.model_info = f"{output_dir}/model_info.json"

        # self.general_library_path = Path(ROOT / "summary/general_library.json")
        self.general_library_path = f"{ROOT}/summary/general_library.xml"
        self.metacyc_rules = Path(ROOT / "summary/MetacycRule.xlsx")
        self.modelseed_met = Path(ROOT / "summary/modelseed_met.xlsx")
        self.meta_met = Path(ROOT / "summary/meta_met.xlsx")
        self.kegg_met = Path(ROOT / "summary/kegg_met.xlsx")
        self.bigg_metacyc = Path(ROOT / "summary/bigg_metacyc.xlsx")
        self.modelseed_rxn = Path(ROOT / "summary/modelseed_reactions.xlsx")
        self.virtual_met = Path(ROOT / "summary/virtual_met.xlsx")
        # print('flux_path:', self.flux_path)
        self.test = 'tmp/literature/iAK888/big_model.xml'
        self.ROOT = ROOT
        self.nadh_png = f"{output_dir}/nadh.png"
        self.atp_png = f"{output_dir}/atp.png"
        self.net_png = f"{output_dir}/net.png"
        self.yield_png = f"{output_dir}/yield.png"
        self.yield2_png = f"{output_dir}/yield2.png"

    def _from_disk(self, filepath, format):
        ...

    def from_yaml(self, filepath):
        self._from_disk(filepath, format='yaml')

    def from_json(self, filepath):
        self._from_disk(filepath, format='json')

    def _to_disk(self, outpath, format):
        if format == 'csv':
            self.__dict__
        elif format == 'tsv':
            self.df.to_csv(outpath, sep='\t')
        elif format == 'json':
            self.df.to_json()
        

    def to_yaml(self, outpath):
        self._to_disk(outpath, format='yaml')

    def to_json(self, outpath):
        self._to_disk(outpath, format='json')


# cfg = Config()
# global cfg