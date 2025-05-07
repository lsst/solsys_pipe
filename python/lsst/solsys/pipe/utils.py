import heliolinx.heliolinx as hl
import heliolinx.solarsyst_dyn_geo as solardg
import numpy as np


def df2numpy(df, dtypename):
    sa = getattr(hl, f"create_{dtypename}")(len(df))
    for col in df.columns:
        sa[col] = df[col]
    return sa


def df_to_numpy(df, dtype):
    sa = np.zeros(len(df), dtype=dtype)
    for col in df.columns:
        sa[col] = df[col]
    return sa


def make_hldet(df):
    return df_to_numpy(df, dtype=solardg.hldet)


def make_hlimage(df):
    return df_to_numpy(df, dtype=solardg.hlimage)
