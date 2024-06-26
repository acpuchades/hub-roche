#!/usr/bin/env python3

import io
import sqlite3

from pathlib import Path
from argparse import ArgumentParser

import vcf
import numpy as np
import pandas as pd

from helpers import normalize_sample_id

def ufela_parse_date(x, **kwargs):
    return pd.to_datetime(x, format='%d-%m-%Y', **kwargs)

def make_argument_parser():
    parser = ArgumentParser()
    parser.add_argument('--vcf', type=Path, required=True)
    parser.add_argument('--samples', type=Path, required=True)
    parser.add_argument('--database', type=Path, required=True)
    parser.add_argument('--output', '-o', choices=['patients', 'controls', 'both'], required=True)
    return parser

parser = make_argument_parser()
args = parser.parse_args()

with open(args.vcf, 'r') as f:
    ufela_vcf = vcf.Reader(f)
    ufela_samples = (
        pd
            .read_excel(args.samples, skiprows=1)
            .rename(columns={
                'Código caso NorayBanks': 'id_noraybanks',
                'NHC*': 'nhc',
            })
            [[
                'id_noraybanks', 'nhc'
            ]]
            .assign(nhc=lambda df: pd.to_numeric(df.nhc, errors='coerce'))
            .drop_duplicates()
    )
    ufela_samples = ufela_samples[
        ufela_samples['id_noraybanks'].isin(ufela_vcf.samples)
    ]

with sqlite3.connect(args.database) as ufela_db:

    ufela_db_patients = (
        pd
            .read_sql("SELECT * FROM pacientes", ufela_db)
            .assign(
                nhc=lambda df: df.nhc.astype(int),
                fecha_nacimiento=lambda df: ufela_parse_date(df.fecha_nacimiento, errors='coerce'),
                fecha_exitus=lambda df: ufela_parse_date(df.fecha_exitus, errors='coerce'),
            )
            .drop(columns=['id', 'created_datetime', 'updated_datetime'])
    )

    ufela_db_clinical = (
        pd
            .read_sql("SELECT * FROM datos_clinicos", ufela_db)
            .assign(
                fecha_inicio_clinica=lambda df: ufela_parse_date(df.fecha_inicio_clinica, errors='coerce'),
                fecha_diagnostico_ELA=lambda df: ufela_parse_date(df.fecha_diagnostico_ELA, errors='coerce')
            )
            .drop(columns=['id', 'created_datetime', 'updated_datetime'])
    )

    ufela_db_alsfrs = (
        pd
            .read_sql("SELECT * FROM esc_val_ela", ufela_db)
            .rename(columns={'fecha_visita_esc_val_ela': 'fecha_visita'})
            .assign(fecha_visita = lambda df: ufela_parse_date(df.fecha_visita, errors='coerce'))
    )

ufela_patients = (
    pd
        .merge(ufela_samples, ufela_db_patients, on='nhc')
        .merge(ufela_db_clinical, on='pid')
        .assign(
            edad_inicio = lambda df: np.floor(
                (df.fecha_inicio_clinica - df.fecha_nacimiento).dt.days / 356.25
            ),
            fenotipo = lambda df: np.select([
                df.fenotipo_al_diagnostico.isin(['ELA Bulbar', 'ELA Espinal', 'ELA Respiratoria']),
                df.fenotipo_al_diagnostico == 'Esclerosis Lateral Primaria (ELP)',
                df.fenotipo_al_diagnostico == 'Atrofia Muscular Progresiva (AMP)',
                df.fenotipo_al_diagnostico == 'Parálisis bulbar progresiva',
                df.fenotipo_al_diagnostico == 'Flail-Arm',
                df.fenotipo_al_diagnostico == 'Flail-Leg',
            ], ['ALS', 'PLS', 'PMA', 'PBP', 'Flail-Arm', 'Flail-Leg']),
            forma_inicio = lambda df: np.select([
                df.fenotipo_al_diagnostico == 'ELA Bulbar',
                df.fenotipo_al_diagnostico == 'ELA Espinal',
                df.fenotipo_al_diagnostico == 'ELA Respiratoria',
            ], ['bulbar', 'spinal', 'respiratory']),
            retraso_diagnostico = lambda df: (
                df.fecha_diagnostico_ELA - df.fecha_inicio_clinica
            ).dt.days / 30,
            tiempo_supervivencia = lambda df: (
                df.fecha_exitus - df.fecha_inicio_clinica
            ).dt.days / (12*30),
        )
)

ufela_controls = (
    ufela_samples[~ufela_samples.nhc.isin(ufela_patients.nhc)]
)

if args.output == 'controls':

    output = ufela_controls.copy()
    output['FID'] = output['IID'] = output.id_noraybanks
    output = output[['FID', 'IID']].sort_values(['FID', 'IID'])

elif args.output == "patients":

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
            (df.fecha_visita - df.fecha_inicio_clinica).dt.days >= (6 * 30),
            (48 - df.total) / ((df.fecha_visita - df.fecha_inicio_clinica).dt.days / 30),
            np.nan
        ).astype(float)
    )

    alsfrs_baselines = (
        ufela_alsfrs[['pid', 'fecha_visita', 'delta_fs']]
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
        ], on='pid', how='left'
    )

    output = ufela_patients.copy()
    output['FID'] = output['IID'] = output.id_noraybanks
    output['NHC'] = output.nhc
    output['ALS'] = np.where(output.fenotipo == 'ALS', 2, 1)
    output['PLS'] = np.where(output.fenotipo == 'PLS', 2, 1)
    output['PMA'] = np.where(output.fenotipo == 'PMA', 2, 1)
    output['PBP'] = np.where(output.fenotipo == 'PBP', 2, 1)
    output['ageatonset'] = output.edad_inicio.fillna(-9).astype(int)
    output['bulbar_onset'] = np.where(output.forma_inicio == 'bulbar', 2, 1)
    output['spinal_onset'] = np.where(output.forma_inicio == 'spinal', 2, 1)
    output['respiratory_onset'] = np.where(output.forma_inicio == 'respiratory', 2, 1)
    output['dxdelay'] = output.retraso_diagnostico.round(1)
    output['deltafs'] = output.delta_fs.round(2)
    output['timetodeath'] = output.tiempo_supervivencia.round(1)
    output = output[[
        'FID', 'IID', 'NHC', 'ALS', 'PLS', 'PMA', 'PBP',
        'bulbar_onset', 'spinal_onset', 'respiratory_onset',
        'ageatonset', 'dxdelay', 'deltafs', 'timetodeath',
    ]].sort_values(['FID', 'IID'])

elif args.output == 'both':

    output = ufela_samples.copy()
    output['FID'] = output['IID'] = output.id_noraybanks
    output['ALS'] = np.where(output.nhc.isin(ufela_patients.nhc), 2, 1)
    output = output[['FID', 'IID', 'ALS']].sort_values(['FID', 'IID'])

outbuf = io.StringIO()
output.to_csv(outbuf, sep='\t', index=False, na_rep='-9')
print(outbuf.getvalue(), end='')