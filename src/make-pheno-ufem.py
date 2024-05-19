#!/usr/bin/env python3

from pathlib import Path
from argparse import ArgumentParser

import io
import numpy as np
import pandas as pd

from helpers import normalize_sample_id

def make_argument_parser():
    parser = ArgumentParser()
    parser.add_argument('--samples', type=Path, required=True)
    parser.add_argument('--database', type=Path, required=True)
    parser.add_argument('--output', choices=['patients', 'unmerged'], default='patients')
    return parser

def read_edmus_file(path, prefix):
    files = list(path.glob(f'{prefix}*.txt'))
    if len(files) != 1:
        return None
    fname = files[0]
    return pd.read_csv(fname, encoding='utf-16', sep='\t')

def date_from_edmus(x):
    return pd.to_datetime(x, dayfirst=True)


parser = make_argument_parser()
args = parser.parse_args()

ufem_samples = (
    pd
        .read_excel(args.samples)
        .rename(columns={
            'id': 'sample_id',
            'fecha DNA': 'date',
            'Unnamed: 2': 'nhc',
        })
        .dropna(subset=['sample_id', 'nhc'])
        .query('nhc != "RIS"')
        .assign(
            nhc = lambda df: pd.to_numeric(df.nhc, errors='coerce'),
            sample_id = lambda df: df.sample_id.map(normalize_sample_id, na_action='ignore')
        )
)

ufem_db_personal = read_edmus_file(args.database, 'BCN4-Personal')
ufem_db_personal['MS Onset'] = date_from_edmus(ufem_db_personal['MS Onset'])
ufem_db_personal['Irreversible DSS 3'] = date_from_edmus(ufem_db_personal['Irreversible DSS 3'])
ufem_db_personal['Irreversible DSS 6'] = date_from_edmus(ufem_db_personal['Irreversible DSS 6'])

ufem_db_diagnosis = read_edmus_file(args.database, 'BCN4-Diagnosis')
ufem_db_diagnosis['MS Onset'] = date_from_edmus(ufem_db_diagnosis['MS Onset'])
ufem_db_diagnosis['Progression Onset'] = date_from_edmus(ufem_db_diagnosis['Progression Onset'])

ufem_db_clinical = read_edmus_file(args.database, 'BCN4-Clinical')
ufem_db_clinical['Date'] = date_from_edmus(ufem_db_clinical['Date'])
ufem_followups = ufem_db_clinical.sort_values(['Patient ID', 'Date'])

ufem_db_episodes = read_edmus_file(args.database, 'BCN4-Episodes')
ufem_db_episodes['Date'] = date_from_edmus(ufem_db_episodes['Date'])

ufem_patients = (
    pd
        .merge(ufem_db_personal.drop(columns=['MS Onset']), ufem_db_diagnosis, on='Patient ID')
        .assign(nhc=lambda df: pd.to_numeric(df['Other identifier'], errors='coerce'))
        .merge(
            ufem_followups
                .groupby('Patient ID')
                .head(1)
                .rename(columns={'Date': 'First Assessment'})
                [[
                    'Patient ID', 'First Assessment'
                ]],
            on='Patient ID'
        )
        .merge(
            ufem_followups
                .groupby('Patient ID')
                .tail(1)
                .rename(columns={'Date': 'Last Assessment'})
                [[
                    'Patient ID', 'Last Assessment'
                ]],
            on='Patient ID'
        )
)

ufem_relapses = ufem_db_episodes.merge(ufem_patients[['Patient ID', 'MS Onset', 'First Assessment']], on='Patient ID')
ufem_relapses_y1 = ufem_relapses[ufem_relapses['Date'] <= (ufem_relapses['MS Onset'] + pd.Timedelta(days=1*12*30))]
ufem_relapses_y1_2 = ufem_relapses[ufem_relapses['Date'] <= (ufem_relapses['MS Onset'] + pd.Timedelta(days=2*12*30))]
ufem_relapses_y1_3 = ufem_relapses[ufem_relapses['Date'] <= (ufem_relapses['MS Onset'] + pd.Timedelta(days=3*12*30))]
ufem_relapses_y1_5 = ufem_relapses[ufem_relapses['Date'] <= (ufem_relapses['MS Onset'] + pd.Timedelta(days=5*12*30))]

ufem_patients = ufem_patients.assign(
    ms_phenotype=lambda df: np.select([
        df['Disease Course'].isin([1, 2]),
        df['Disease Course'].isin([3, 4]),
        df['Disease Course'].isin([5, 6, 7]),
    ], ['RR', 'SP', 'PP']),
    disease_duration=lambda df: (df['Last Assessment'] - df['MS Onset']).dt.days / (12*30),
    timetoiedss3=lambda df: (df['Irreversible DSS 3'] - df['MS Onset']).dt.days / (12*30),
    timetoiedss6=lambda df: (df['Irreversible DSS 6'] - df['MS Onset']).dt.days / (12*30),
    timetoprogression=lambda df: (date_from_edmus(df['Progression Onset']) - df['MS Onset']).dt.days / (12*30),
    arr=lambda df: np.where(
        df['First Assessment'] <= (df['MS Onset'] + pd.Timedelta(days=12*30)),
        df[['Patient ID']].merge(
            ufem_relapses.groupby('Patient ID').count(), on='Patient ID', how='left'
        )['Episode ID'] / df.disease_duration, np.nan
    ),
    arr_y1=lambda df: np.where(
        (df.disease_duration >= 1) & (df['First Assessment'] <= (df['MS Onset'] + pd.Timedelta(days=12*30))),
        df[['Patient ID']].merge(
            ufem_relapses_y1.groupby('Patient ID').count(), on='Patient ID', how='left'
        )['Episode ID'], np.nan
    ),
    arr_y1_2=lambda df: np.where(
        (df.disease_duration >= 2) & (df['First Assessment'] <= (df['MS Onset'] + pd.Timedelta(days=12*30))),
        df[['Patient ID']].merge(
            ufem_relapses_y1_2.groupby('Patient ID').count(), on='Patient ID', how='left'
        )['Episode ID'] / 2, np.nan
    ),
    arr_y1_3=lambda df: np.where(
        (df.disease_duration >= 3) & (df['First Assessment'] <= (df['MS Onset'] + pd.Timedelta(days=12*30))),
        df[['Patient ID']].merge(
            ufem_relapses_y1_3.groupby('Patient ID').count(), on='Patient ID', how='left'
        )['Episode ID'] / 3, np.nan
    ),
    arr_y1_5=lambda df: np.where(
        (df.disease_duration >= 5) & (df['First Assessment'] <= (df['MS Onset'] + pd.Timedelta(days=12*30))),
        df[['Patient ID']].merge(
            ufem_relapses_y1_5.groupby('Patient ID').count(), on='Patient ID', how='left'
        )['Episode ID'] / 5, np.nan
    ),
)

if args.output == 'patients':
    output = pd.merge(ufem_samples, ufem_patients, on='nhc')
    output['FID'] = output['IID'] = output.sample_id
    output['RR_SP'] = np.where(output.ms_phenotype.isin(['RR', 'SP']), 2, 1)
    output['PP'] = np.where(output.ms_phenotype == 'PP', 2, 1)
    output['ageatonset'] = np.where(output['Age at onset'] > 0, output['Age at onset'], -9)
    output['duration'] = output['disease_duration'].round(1)
    output['arr'] = output['arr'].round(1)
    output['arr_y1'] = output['arr_y1'].round(1)
    output['arr_y1_2'] = output['arr_y1_2'].round(1)
    output['arr_y1_3'] = output['arr_y1_3'].round(1)
    output['arr_y1_5'] = output['arr_y1_5'].round(1)
    output['timetoiedss3'] = output['timetoiedss3'].round(1)
    output['timetoiedss6'] = output['timetoiedss6'].round(1)
    output['timetoprogression'] = output['timetoprogression'].round(1)
    output = output[[
        'FID', 'IID', 'RR_SP', 'PP', 'ageatonset', 'duration',
        'timetoiedss3', 'timetoiedss6', 'timetoprogression',
        'arr', 'arr_y1', 'arr_y1_2', 'arr_y1_3', 'arr_y1_5',
    ]].sort_values(['FID', 'IID'])

elif args.output == 'unmerged':
    output = ufem_samples.merge(
        ufem_db_personal.assign(nhc=lambda df: pd.to_numeric(df['Other identifier'], errors='coerce')),
        on='nhc', how='left'
    )
    output = output.loc[output['Patient ID'].isna(), ['sample_id', 'nhc', 'date']]

outbuf = io.StringIO()
output.to_csv(outbuf, sep='\t', index=False, na_rep='-9')
print(outbuf.getvalue(), end='')