# 匿名复现仓库

本仓库用于应急无人机部署、有限电池池与条件性设备价值的模型复现。代码和21个可分发合成基准输入随仓库提供；定日和Rasuwa仅提供公开数据获取、处理和模型构造脚本，不分发外部原始数据。没有受限第三方工作簿或无法确认重用许可的损毁建筑矢量，没有预存求解方案。

## 环境

已在Python 3.11.2验证。建议Python 3.11，依赖固定在requirements.txt；MILP使用SciPy附带的HiGHS，不需要商业求解器。安装：

```powershell
python -m venv .venv
.venv\Scripts\python -m pip install -r requirements.txt
```

Linux/macOS可使用虚拟环境的`bin/python`；脚本只依赖仓库相对目录，不依赖特定电脑。

## 一键复现

先核验21任务、264个配送请求及输入SHA：

```powershell
python reproduce.py --mode verify
```

完整公开基准复算与独立验收：

```powershell
python reproduce.py --mode benchmark --budget 60
```

共252配置：168有限池、42等基数48点选址、42即时就绪控制。即时就绪控制放松充电准备度，不代表真实充电技术。运行可能需数小时，终端逐配置报告进度，输出位于忽略跟踪的results/reports目录。每个求解方案由另一进程重新验收能量、充电、资源占用与连续通信。

`delivery request`是具有指定载荷和位置的不可分割包裹，随后简称request；task表示整套独立基准场景。task-aware选址不是计量单位。定日communication profile表示一个地面终端及完整的规定服务窗口。

固定预算、机器速度和求解时序可能改变可行下界；重新搜索不保证逐字节返回论文中的同一方案。有效上界是固定运输条件的必要条件界，不是受限MILP对偶。未闭合gap不能称为全问题最优或最优设备增量。

## 公开案例数据处理

合法取得源文件后，目录需含CN.zip/NP.zip（或解压txt）、AOI01/AOI02的初始分级GPKG，以及相关GLO-30 TIFF。可下载公开源：

```powershell
python case_processing/fetch_public_data.py --output public_data
python reproduce.py --mode cases --data-root public_data --budget 60
```

案例命令先处理源数据，再建立定日10个通信剖面、Rasuwa7个标准配送请求，两假设基地、三个部署域、四个库存条件，共36配置。源数量或字段与既定协议不符会停止，不补造点位。对现有合法裁切也可提供`dingri_DEM.tif`与`rasuwa_DEM.tif`。官方源变更或下载失败需要重新取得正确资料，代码不会改用随机数据。

人口需求构造：

```powershell
python case_processing/prepare_population_demands.py --exposure public_data/UNOSAT_population.xlsx --buildings source_data/Rasuwa确认受损建筑点.csv
```

用UNOSAT行政暴露量和受损建筑比例构造20/96/192/479请求，以及300/1800/7200秒窗口；假定每栋建筑暴露均匀，不是村级观测人口、救援订单或实际行动效果。

## 数据来源与许可

- 合成21任务由解析地形和固定种子生成；自有代码及随附合成输入按MIT条款开放。真实事件原始源不随本仓库分发。
- GeoNames公开国家导出：https://download.geonames.org/export/dump/；CC BY 4.0。公开地名不证明建筑受损或村庄孤立。
- Copernicus GLO-30：https://doi.org/10.5270/ESA-c5d3d65；专用DEM免费使用条款及必要署名见SOURCE_NOTICES.md，不称为MIT数据。
- Copernicus EMSR927：https://data.humdata.org/dataset/npl-flood-emsr927；HDX记录CC BY，另遵守Copernicus EMS来源条款。仅Destroyed/Damaged建筑被计入，Possibly damaged不升级解释。
- UNOSAT 4257：https://unosat.org/products/4257；行政暴露表的来源记录为CC BY-SA，不把它重新许可为MIT。

MIT仅适用于自有软件及明确开放的合成输入；下载的外部资料和Python依赖保留原许可。SOURCE_NOTICES.md保留数据提供方必要署名，这些署名不含稿件作者身份。

仓库没有个人联系方式、账户令牌、远程Git配置或发表DOI。未发布、未推送。
