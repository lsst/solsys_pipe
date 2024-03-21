import numpy as np
import heliolinx.solarsyst_dyn_geo as solardg
import heliolinx.heliolinx as hl

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

def det_to_numpy(df):
    return df_to_numpy(df, dtype=solardg.hldet)

def vis_to_numpy(df):
    return df_to_numpy(df, dtype=solardg.hlimage)


