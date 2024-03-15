#!/usr/bin/env python3

import io
import sqlite3
import numpy as np
import pandas as pd

from argparse import ArgumentParser

def ufela_parse_date(x, **kwargs):
    return pd.to_datetime(x, format='%d-%m-%Y', **kwargs)

def make_argument_parser():
    parser = ArgumentParser()
    parser.add_argument("--samples", required=True)
    parser.add_argument("--patients", required=True)
    return parser

parser = make_argument_parser()
args = parser.parse_args()

ufela_samples = (
    pd
        .read_excel("data/ufela/samples-20240201.xlsx")
        .rename(columns={
            'Unnamed: 1': 'sample_ids_and_comments',
            'Nombre': 'patient_name',
            'SAP': 'nhc',
            'Dtco': 'phenotype',
            'Fecha muestra': 'date'
        })
        .assign(
            # explode samples with more than one sample_id while preserving comments
            sample_id = lambda df: 'ELA-' + df.sample_ids_and_comments.str.extract(r'^ELA[- ]?0*(\d+\s?(?:\(\d+\))?)'),
            sample_id2 = lambda df: 'ELA-' + df.sample_ids_and_comments.str.extract(r'/\s*ELA[- ]0*?(\d+\s?(?:\(\d+\))?)'),
            comments = lambda df: df.sample_ids_and_comments.str.extract(r'\(([^0-9][^)]+)\)\s*$')
        )
)

ufela_samples = (
    pd
        .concat([
            ufela_samples[ufela_samples.sample_id2.isna()].drop(columns='sample_id2'),
            ufela_samples[ufela_samples.sample_id2.notna()].drop(columns='sample_id').rename(columns={'sample_id2': 'sample_id'})
        ])
        .assign(
            sample_id = lambda df: df.sample_id.str.replace(' ', ''),
            order = lambda df: df.sample_id.str.extract(r'ELA-(\d+)').astype(float)
        )
        .drop_duplicates(['sample_id', 'date'])
        .dropna(subset='sample_id')
        .sort_values('order')
        [[
            'order', 'sample_id', 'date', 'nhc', 'patient_name', 'phenotype', 'comments'
        ]]
)

ufela_controls = (
    pd
        .concat([
            ufela_samples[ufela_samples.patient_name.str.contains('control', case=False, na=False)],
            ufela_samples[ufela_samples.phenotype.str.contains('control', case=False, na=False)]
        ])
        .drop_duplicates('sample_id')
        .drop(columns='date')
)

with sqlite3.connect(args.patients) as ufela_db:
    ufela_db_patients = pd.read_sql("SELECT * FROM pacientes", ufela_db)
    ufela_db_clinical = (
        pd
            .read_sql("SELECT * FROM datos_clinicos", ufela_db)
            .assign(
                fecha_inicio_clinica = lambda df: ufela_parse_date(df.fecha_inicio_clinica, errors='coerce'),
                fecha_diagnostico_ELA = lambda df: ufela_parse_date(df.fecha_diagnostico_ELA, errors='coerce')
            )
    )
    ufela_db_alsfrs = (
        pd
            .read_sql("SELECT * FROM esc_val_ela", ufela_db)
            .rename(columns={ 'fecha_visita_esc_val_ela': 'fecha_visita' })
            .assign(fecha_visita = lambda df: ufela_parse_date(df.fecha_visita, errors='coerce'))
    )

ufela_patients = (
    pd
        .merge(
            ufela_db_patients.drop(columns=['id', 'created_datetime', 'updated_datetime']),
            ufela_db_clinical.drop(columns=['id', 'created_datetime', 'updated_datetime'])
        )
        .assign(nhc = lambda df: df.nhc.astype(int))
        .merge(ufela_samples, on='nhc')
        .drop_duplicates('sample_id')
        .dropna(subset='sample_id')
)

ufela_alsfrs = ufela_db_alsfrs.rename(columns={'fecha_visita_esc_val_ela': 'fecha_visita'})
ufela_alsfrs.iloc[:, 3:16] = ufela_alsfrs.iloc[:, 3:16].transform(
    lambda x: pd.to_numeric(x, errors='coerce').where(lambda x: x.between(0, 4))
)
ufela_alsfrs = ufela_alsfrs.merge(
    ufela_db_clinical[['pid', 'fecha_inicio_clinica']], on = 'pid'
)
ufela_alsfrs = ufela_alsfrs.assign(
    cortar = lambda df: np.select([
        df.cortar_sin_peg.notna() & df.cortar_con_peg.isna(),
        df.cortar_con_peg.notna() & df.cortar_sin_peg.isna(),
        df.cortar_con_peg == df.cortar_sin_peg,
    ], [df.cortar_sin_peg, df.cortar_con_peg, df.cortar_con_peg]),
    total_bulbar = lambda df: df.lenguaje + df.salivacion + df.deglucion,
    total_motor_fino = lambda df: df.escritura + df.cortar + df.vestido,
    total_motor_grosero = lambda df: df.cama + df.caminar + df.subir_escaleras,
    total_respiratorio = lambda df: df.disnea + df.ortopnea + df.insuficiencia_respiratoria,
    total = lambda df: df.total_bulbar + df.total_motor_fino + df.total_motor_grosero + df.total_respiratorio,
    delta_fs = lambda df: np.where(
        (df.fecha_visita - df.fecha_inicio_clinica).dt.days >= 6,
        (48 - df.total) / ((df.fecha_visita - df.fecha_inicio_clinica).dt.days / 30),
        np.nan
    )
)

alsfrs_baselines = (
    ufela_alsfrs
        [[
            'pid', 'fecha_visita', 'delta_fs'
        ]]
        .dropna()
        .sort_values(['pid', 'fecha_visita'])
        .groupby('pid')
        .head(1)
        .rename(columns={
            'total': 'alsfrs_total_basal',
            'fecha_visita': 'fecha_alsfrs_basal',
        })
        .merge(ufela_db_clinical[['pid', 'fecha_diagnostico_ELA']], on = 'pid')
)

# do not include deltafs for second opinions and patient transfers
ufela_patients = ufela_patients.merge(
    alsfrs_baselines.loc[
        ((alsfrs_baselines.fecha_alsfrs_basal - alsfrs_baselines.fecha_diagnostico_ELA).dt.days / 30).between(0, 6),
        ['pid', 'delta_fs']
    ], on='pid'
)

ufela_controls['FID'] = ufela_controls['IID'] = ufela_controls.sample_id
ufela_controls['ALS'] = 1
ufela_controls['DFS'] = -9

ufela_patients['FID'] = ufela_patients['IID'] = ufela_patients.sample_id
ufela_patients['ALS'] = 2
ufela_patients['DFS'] = ufela_patients.delta_fs.astype(float).round(2).fillna(-9)

output = (
    pd
        .concat([ufela_controls, ufela_patients])
        .sort_values('order')
        [[
            'FID', 'IID', 'ALS', 'DFS'
        ]]
)

outbuf = io.StringIO()
output.to_csv(outbuf, sep='\t', index=False)
print(outbuf.getvalue())