import setuptools
with open("README.md", "r") as fh:
    long_description = fh.read()
setuptools.setup(
    name="mqc",  # 模块名称
    version="2.0",  # 当前版本
    author="dengxiao",  # 作者
    author_email="dengxg@tib.cas.cn",  # 作者邮箱
    description="MQC is a Genome-scale metabolic network model quality control tool",  # 模块简介
    long_description=long_description,  # 模块详细介绍
    long_description_content_type="text/markdown",  # 模块详细介绍格式
    url="http://172.16.25.29/dengxiao/mqc/",  # 模块GitLab地址
    packages=setuptools.find_packages(where="mqc"),  # 自动找到项目中导入的模块
    package_dir={"": "mqc"},
    package_data={
        "mqc.summary": ["*.xlsx"],
        "mqc.summary": ["*.json"],
        "mqc.summary": ["*.xml"]
    },
    include_package_data=True,
    # 模块相关的元数据
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    # 依赖模块
    install_requires=[
        "httpx==0.23.0",
        "numpy==1.22.4",
        "cobra==0.25.0",
        "matplotlib==3.7.2",
        "openpyxl"
    ],
    python_requires='>=3',
    entry_points={
    'console_scripts': [
        'mqc = mqc.main:main'
    ]
},
)