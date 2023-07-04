""" Update canid database
This script aims to locate all the raw and processed data in ERDA and Mjolnir.
"""

import os
import paramiko
import pandas as pd
import numpy as np
from datetime import date


class ERDA:
    """ ERDA sftp connection"""

    def __init__():
        pass

    def connect() -> ERDA_session:
        """ Connect to ERDA and return a session. """

    def locate() -> list:
        """ Locate a folder with given list of id"""

    def download() -> None:
        """ Download list of files """

    def disconnetct() -> None:
        """ Disconnect session"""

    def check_alive() -> Bool:
        """ Check if session is alive."""


def bam_df(no_rows: int) -> pd.DataFrame:
    """ create a bam df
    One BAM file per row

    Args:
        no_rows: Desired number of data rows.

    Returns:
        pd.DataFrame:
          Dataframe with NUll as default values. Data columns are as follows

          ========  =========================================================
          bam       str. bam file name
          dir       str. bam file location, full directory
          id        str. Main sample id
          ref       str. Reference version mapped to
          size      int. file size
          date      datetime64
          idx       bool: bam.bai exist
          idx_file  str. index file location, full directory
          ========  =========================================================
    """
    df = pd.DataFrame({
        "bam": np.empty(no_rows, dtype='object'),
        "dir": np.empty(no_rows, dtype='object'),
        "id": np.empty(no_rows, dtype='object'),
        "ref": np.empty(no_rows, dtype='object'),
        "size": np.empty(no_rows, dtype='int'),
        "date": np.repeat(np.datetime64(date.today()), no_rows),
        "idx": np.empty(no_rows, dtype='bool'),
        "idx_file": np.empty(no_rows, dtype='object')
    })

    return df


def raw_df(no_rows: int) -> pd.DataFrame:
    """ create a raw data df
    One fastq per row

    Args:
        no_rows: Desired number of data rows.

    Returns:
        pd.DataFrame:
          Dataframe with NUll as default values. Data columns are as follows

          ========  =========================================================
          fq        str. fastq file
          gz        bool. gzipped
          paired    int. 0,1,2
          id        str. Main sample id
          size      int. file size
          date      datetime64
          batch     str. sequencing batch
          platform  str. sequencing platform
          ========  =========================================================
    """
    df = pd.DataFrame({
        "fq": np.empty(no_rows, dtype='object'),
        "gz": np.empty(no_rows, dtype='bool'),
        "paired": np.empty(no_rows, dtype='int'),
        "id": np.empty(no_rows, dtype='object'),
        "size": np.empty(no_rows, dtype='int'),
        "date": np.repeat(np.datetime64(date.today()), no_rows),
        "batch": np.empty(no_rows, dtype='object'),
        "platform": np.empty(no_rows, dtype='object')
    })

    return df


def sample_df(no_rows: int) -> pd.DataFrame:
    """ create a sample data df
    One sample per row

    Args:
        no_rows: Desired number of data rows.

    Returns:
        pd.DataFrame:
          Dataframe with NUll as default values. Data columns are as follows
          ========  =========================================================
          id        str. Main sample id
          id_alt    list of alternative ids
          dp        sequencing coverage
          batch     list of sequencing batches
          raw       list of raw data in erda
          bam       list of bam file in mjolnir
          projects  list of projects
          upload    bool. uploaded to ENA
          ========  =========================================================
    """
    df = pd.DataFrame({
        "id": np.empty(no_rows, dtype='object'),
        "id_alt": [[''] for i in range(no_rows)],
        "dp": np.empty(no_rows, dtype='float16'),
        "batch": [[''] for i in range(no_rows)],
        "raw": [[''] for i in range(no_rows)],
        "bam": [[''] for i in range(no_rows)],
        "projects": [[''] for i in range(no_rows)],
        "upload": np.empty(no_rows, dtype='bool')
    })

    return df

def read_csv(csv: list[str], sample_df: pd.DataFrame) -> None:
    """ Read csv(list) into sample_df

    Args:
        csv: list of csv files
        sample_df: sample data frame to update

    Returns:
        None
    """
    for file_csv in csv:
        tmp_df = pd.read_csv(file_csv)
        print(f'File read: {file_csv} with {tmp_df.shape[0]} rows')

        print(f'Removed {test} rows with missing Lab ID')






