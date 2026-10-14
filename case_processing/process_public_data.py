"""处理合法取得的GeoNames、Copernicus DSM和EMSR927公开源，不生成替代观测。"""
from pathlib import Path
import argparse,json,sqlite3,struct,zipfile,shutil
import numpy as np,pandas as pd,rasterio
from rasterio.merge import merge
from rasterio.windows import from_bounds,Window
from pyproj import Geod,Transformer
from shapely import wkb
ROOT=Path(__file__).resolve().parent;OUT=ROOT.parent/'source_data'
COLS=['geonameid','name','asciiname','alternatenames','lat','lon','feature_class','feature_code','country','cc2','admin1','admin2','admin3','admin4','population','elevation','dem','timezone','modification_date']
def geonames(folder,country):
 p=folder/(country+'.txt');z=folder/(country+'.zip')
 if p.exists():source=p
 elif z.exists():
  with zipfile.ZipFile(z) as f:
   with f.open(country+'.txt') as s:return pd.read_csv(s,sep='\t',names=COLS,low_memory=False)
 else:raise FileNotFoundError('缺少公开国家GeoNames文件：'+country)
 return pd.read_csv(source,sep='\t',names=COLS,low_memory=False)
def decode_geom(blob):
 assert blob[:2]==b'GP';flags=blob[3];envelope=(flags>>1)&7;length={0:0,1:32,2:48,3:48,4:64}[envelope];return wkb.loads(blob[8+length:])
def buildings(folder):
 rows=[]
 for aoi in ['aoi01','aoi02']:
  files=[p for p in folder.rglob('*.gpkg') if aoi in p.name.lower()];assert len(files)==1,'每个AOI需一个初始分级GPKG'
  with sqlite3.connect(files[0]) as c:
   layers=c.execute('SELECT table_name,srs_id FROM gpkg_contents').fetchall()
   for table,srs in layers:
    if not table.lower().startswith('builtup'):continue
    assert all(x.isalnum() or x=='_' for x in table)
    cols=[x[1] for x in c.execute('PRAGMA table_info('+table+')')];geom=c.execute('SELECT column_name FROM gpkg_geometry_columns WHERE table_name=?',(table,)).fetchone()[0];fid=cols[0];damage=next(x for x in cols if x.lower()=='damage_gra')
    if int(srs)!=4326:
     definition=c.execute('SELECT definition FROM gpkg_spatial_ref_sys WHERE srs_id=?',(srs,)).fetchone()[0];tr=Transformer.from_crs(definition,4326,always_xy=True)
    else:tr=None
    for idx,grade,blob in c.execute('SELECT '+fid+','+damage+','+geom+' FROM '+table):
     if grade not in ['Destroyed','Damaged']:continue
     g=decode_geom(blob)
     if g.geom_type!='Point':g=g.centroid
     lon,lat=g.x,g.y
     if tr:lon,lat=tr.transform(lon,lat)
     rows.append({'经度':lon,'纬度':lat,'来源文件':files[0].name,'图层':table,'来源FID':int(idx),'damage_gra':grade})
 assert len(rows)==760,'公开产品要素数与本次初始产品不符，停止而非伪造需求'
 return pd.DataFrame(rows)
def terrain(folder,name,bounds):
 direct=next(iter(folder.rglob(name+'_DEM.tif')),None);out=OUT/'地形裁切';out.mkdir(exist_ok=True)
 if direct:
  with rasterio.open(direct) as s:z=s.read(1);transform=s.transform;nodata=s.nodata
 else:
  tiles=list(folder.rglob('Copernicus_DSM*.tif'));assert tiles,'缺少公开GLO-30 TIFF'
  datasets=[rasterio.open(p) for p in tiles]
  try:
   mosaic,transform=merge(datasets,bounds=bounds);z=mosaic[0];nodata=datasets[0].nodata
  finally:
   for s in datasets:s.close()
 lon=transform.c+(np.arange(z.shape[1])+.5)*transform.a;lat=transform.f+(np.arange(z.shape[0])+.5)*transform.e
 np.savez_compressed(out/(name+'_dem.npz'),z=z,lon=lon,lat=lat,nodata=-9999. if nodata is None else nodata)
def main(folder):
 OUT.mkdir(exist_ok=True);cn=geonames(folder,'CN');npal=geonames(folder,'NP')
 # 只导出所需地理窗口，不导出无关记录。
 cn=cn[cn.lon.between(87.25,87.75)&cn.lat.between(28.30,28.75)].copy();npal=npal[npal.lon.between(85.08,85.45)&npal.lat.between(27.85,28.30)].copy()
 cn.to_csv(OUT/'CN_案例居民点.csv',index=False,encoding='utf-8-sig');npal.to_csv(OUT/'NP_案例居民点.csv',index=False,encoding='utf-8-sig')
 _,_,d=Geod(ellps='WGS84').inv(np.full(len(cn),87.45),np.full(len(cn),28.50),cn.lon.to_numpy(),cn.lat.to_numpy());dingri=cn[(d<=10000)&cn.feature_code.eq('PPL')].sort_values('geonameid');assert len(dingri)==10,'PPL点位数不符，停止'
 dingri.to_csv(OUT/'定日真实居民点.csv',index=False,encoding='utf-8-sig');buildings(folder).to_csv(OUT/'Rasuwa确认受损建筑点.csv',index=False,encoding='utf-8-sig')
 terrain(folder,'dingri',(87.25,28.30,87.75,28.75));terrain(folder,'rasuwa',(85.08,27.85,85.45,28.30));print('公开来源处理完成：定日10点、EMSR927受损建筑760点。')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--data-root',type=Path,required=True);a=p.parse_args();main(a.data_root)
