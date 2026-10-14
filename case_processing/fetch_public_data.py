"""按固定公开来源下载原始数据；不上传、认证或使用私有数据。"""
from pathlib import Path
import argparse,requests,zipfile,hashlib,json
URLS={
 'CN.zip':'https://download.geonames.org/export/dump/CN.zip',
 'NP.zip':'https://download.geonames.org/export/dump/NP.zip',
 'EMSR927_AOI01.zip':'https://rapidmapping.emergency.copernicus.eu/backend/EMSR927/AOI01/GRA_PRODUCT/EMSR927_AOI01_GRA_PRODUCT_v1.zip?type=gpkg',
 'EMSR927_AOI02.zip':'https://rapidmapping.emergency.copernicus.eu/backend/EMSR927/AOI02/GRA_PRODUCT/EMSR927_AOI02_GRA_PRODUCT_v2.zip?type=gpkg',
 'UNOSAT_population.xlsx':'https://unosat.org/static/unosat_filesystem/4257/UNOSAT_PopulationExposure_20260826_20260827_%20Nepal.xlsx'}
for tile in ['N28_00_E087_00','N27_00_E085_00','N28_00_E085_00']:
 stem='Copernicus_DSM_COG_10_'+tile+'_DEM';URLS[stem+'.tif']='https://copernicus-dem-30m.s3.amazonaws.com/'+stem+'/'+stem+'.tif'
def main(folder):
 folder.mkdir(parents=True,exist_ok=True);records=[]
 for name,url in URLS.items():
  target=folder/name
  if not target.exists():
   with requests.get(url,stream=True,timeout=(30,300)) as r:
    r.raise_for_status();temp=target.with_suffix(target.suffix+'.part')
    with temp.open('wb') as f:
     for chunk in r.iter_content(1048576):f.write(chunk)
    temp.replace(target)
  h=hashlib.sha256(target.read_bytes()).hexdigest();records.append({'file':name,'url':url,'sha256':h})
  if name.startswith('EMSR927'):
   with zipfile.ZipFile(target) as z:
    for info in z.infolist():
     if info.filename.lower().endswith('.gpkg'):
      dest=folder/Path(info.filename).name
      with z.open(info) as source,dest.open('wb') as f:
       while chunk:=source.read(1048576):f.write(chunk)
 (folder/'download_manifest.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf8');print('公开源下载完成；失败会停止，不创建替代输入。')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=Path('public_data'));a=p.parse_args();main(a.output)
