# MQC: Genome-scale metabolic network model quality control tool

-----------------

## What is it?

**MQC** is a Genome-scale metabolic network model quality control tool

## Installing CPLEX Commercial Package

To use this program, you need to install the CPLEX commercial package. Follow these steps:

1. Visit the [IBM website](https://www.ibm.com/analytics/cplex-optimizer) and download the CPLEX installation package suitable for your operating system.

2. Install CPLEX by following the instructions provided with the downloaded installation package.

3. After the installation is complete, ensure that the correct environment variables are set so that your program can find CPLEX. You may need to add the installation directory of CPLEX to the PATH environment variable.


If you encounter any issues, refer to the [IBM CPLEX documentation](https://www.ibm.com/support/knowledgecenter/SSSA5P_12.10.0/COS_KC_home.html) or contact their support team.

OR 

```
$ cd /hpcfs/fhome/mao_zt/project/MCModel/cplex_new/CPLEX221/python  
$ python setup.py install
```
## Install this program. 

You can do this via pip by running the following command:

   ```
   pip install mqc
   import main 
   main.main(model_path, output_path)
   ```
## JSON 转网页

通过安装mqc包进行质控，得到的json文件可以转为网页，便于可视化观看，下面提供3种途径：

1. Windows端打开https://d1e9ilirjteg2w.cloudfront.net/tools/mqc-viewer.zip， 下载mqc-viewer.zip压缩包，解压后打开其中的exe文件，上传质控后的json文件即可

2. 通过网页https://mqc.biodesign.ac.cn/view 上传

3. Mac端打开https://d1e9ilirjteg2w.cloudfront.net/tools/mqc.dmg， 下载macOS的安装包mqc.dmg进行上传

## Run the example

```sh
# get projects
git clone http://172.16.25.29/dengxiao/mqc.git

cd mqc

# create env

conda env create -f environments.yaml

conda activate mqc


# run
python start.py -h

# 结果输出在当前路径
python start.py --file ./mqc/test_data/bigg_data/iCN718.xml 

# 指定结果输出路径
python start.py --file ./mqc/test_data/bigg_data/iCN718.xml -o /tmp/test1
```

input : model file

output : 
    result.json :  self.result_file
    new model file : f"{cfg.ROOT}/{model.id}.xml" or f"{cfg.ROOT}/{model.id}.json" (mqc/utils.py - write_final_model())


