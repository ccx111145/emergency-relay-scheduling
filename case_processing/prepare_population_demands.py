"""从UNOSAT行政暴露汇总按受损建筑占比分配情景请求，不把它称为实测订单。"""
from pathlib import Path
import argparse,json,math
import pandas as pd,numpy as np
from pyproj import Transformer
def main(exposure,buildings,out):
 sheets=pd.read_excel(exposure,sheet_name=None,header=None);found=[]
 for frame in sheets.values():
  for row in frame.itertuples(index=False,name=None):
   if any(isinstance(v,str) and v.strip().lower()=='rasuwa' for v in row):found.append(row)
 assert found,'缺少Rasuwa行政汇总行，不能代填暴露量'
 # 列结构先由使用者核对：暴露值应来自行政区汇总，不能把乡镇或总人口替代。
 numbers=[float(v) for row in found for v in row if isinstance(v,(float,int,np.number)) and not pd.isna(v)]
 values=[v for v in numbers if abs(v-1915.247941872)<1e-3];assert values,'公开表字段变化或不含已核对的Rasuwa暴露量，停止'
 population=values[0];b=pd.read_csv(buildings);x,y=Transformer.from_crs(4326,32645,always_xy=True).transform(b['经度'].to_numpy(),b['纬度'].to_numpy());b['gx']=np.floor(x/2000).astype(int);b['gy']=np.floor(y/2000).astype(int);counts=b.groupby(['gx','gy'],sort=True).size().to_numpy();assert len(counts)==7 and counts.sum()==760
 rows=[]
 for rate in [.01,.05,.10,.25]:
  Q=math.ceil(rate*population);share=(Q-7)*counts/counts.sum();alloc=1+np.floor(share).astype(int);remain=Q-int(alloc.sum());order=sorted(range(7),key=lambda i:(-(share[i]-math.floor(share[i])),i))
  for i in order[:remain]:alloc[i]+=1
  rows.append({'exposed_population':population,'requests_per_person':rate,'requests':Q,'cell_requests':alloc.tolist(),'delivery_window_extensions_seconds':[300,1800,7200],'parcel_kg':1.5,'parcel_volume_m3':.002,'batch_max_requests':20})
 out.parent.mkdir(parents=True,exist_ok=True);out.write_text(json.dumps({'assumption':'行政汇总暴露按受损建筑比例分配，假定每建筑暴露均匀；不是网格人口观测或救援订单','scenarios':rows},ensure_ascii=False,indent=2),encoding='utf8');print('人口情景请求规模20/96/192/479；不是观测需求。')
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('--exposure',type=Path,required=True);p.add_argument('--buildings',type=Path,required=True);p.add_argument('--out',type=Path,default=Path('reports/population_scenarios.json'));a=p.parse_args();main(a.exposure,a.buildings,a.out)
