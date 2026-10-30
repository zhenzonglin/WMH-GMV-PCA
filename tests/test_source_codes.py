import csv
from pathlib import Path
from wmh_gmv.clinical import prepare_clinical
from wmh_gmv.config import load_config

ROOT=Path(__file__).resolve().parents[1]

def write(path,rows):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)

def test_dictionary_specific_yes_no_and_death_codes(tmp_path):
    cfg=load_config(ROOT/'config/project.yaml')
    image=tmp_path/'image.tsv';clinical=tmp_path/'clinical.tsv';out=tmp_path/'out'
    write(image,[dict(subject_key='C3:001',participant_id='001',cohort='C3',age=65,sex=1,
                      icv_ml=1400,wmh_total_ml=8,gmv_total_ml=600)])
    write(clinical,[dict(code_n='001',H_MRS=0,A_NIHSS=4,H_STROKE=1,H_HYPT=2,
        hypertension=1,H_DIAB=1,H_AF=0,TPA=1,ET_MT=2,F3_DEATH=2,F6_DEATH=2,F12_DEATH=1,
        F3_MRS=2,F6_MRS=1,F12_MRS='',IMG_SVD_CMB=99,IMG_SVD_CMB_cSS=1,y1_is=1,y1_is_dd=120)])
    prepare_clinical(image,clinical,cfg,out,'C3')
    with (out/'clinical_wide.tsv').open() as f: row=next(csv.DictReader(f,delimiter='\t'))
    assert float(row['ivt'])==0 and float(row['evt'])==1
    assert float(row['hypertension'])==1 and float(row['diabetes'])==0
    assert float(row['death_m3'])==0 and float(row['mrs_m3'])==2
    assert float(row['mrs_m12'])==6
    assert row['clinical_conflict']=='0'  # H_HYPT=2和canonical=1是同一含义。
    assert row['cmb']=='' and row['source_cmb']=='99'
    assert row['cmb_assessment_status']=='not_assessable' and float(row['css'])==0
    with (out/'events.tsv').open() as f: events=list(csv.DictReader(f,delimiter='\t'))
    event=next(e for e in events if e['event_type']=='ischemic_recurrence')
    assert float(event['relative_days'])==120 and event['time_unit']=='days'
    assert event['event_date_raw']=='' and event['time_valid']=='1'

def test_nonmatching_binary_code_remains_invalid(tmp_path):
    cfg=load_config(ROOT/'config/project.yaml')
    image=tmp_path/'image.tsv';clinical=tmp_path/'clinical.tsv';out=tmp_path/'out'
    write(image,[dict(subject_key='C3:001',participant_id='001',cohort='C3',age=65,sex=1,
                      icv_ml=1400,wmh_total_ml=8,gmv_total_ml=600)])
    write(clinical,[dict(code_n='001',TPA=7)])
    prepare_clinical(image,clinical,cfg,out,'C3')
    with (out/'variable_audit.tsv').open() as f: rows=list(csv.DictReader(f,delimiter='\t'))
    assert next(r for r in rows if r['variable']=='ivt')['invalid']=='1'
