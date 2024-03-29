""" Update canid database
This script aims to locate all the raw and processed data in ERDA and Mjolnir.

Remember next time:
If you only have two files, no loop
"""

import os
import paramiko
import pandas as pd
import numpy as np
import re
from datetime import date
from stat import S_ISREG, S_ISDIR
from datetime import datetime
import argparse


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


def fq_df(no_rows: int) -> pd.DataFrame:
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
          batch     list of batches, dict{'names': str; 'idx': int; 'pool': bool}
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


def clear_id(id_column: pd.Series) -> [pd.Series, pd.Series]:
    """ Clear sample id

    For main ID: AWxxx and MWxxx will be kept
    For alternative IDs, seperators and or () [], AWxxxABC MWxxxABC

    Args:
        A pd.Seriers column

    Returns:
        list(main_id, alternative_id_list)
    """
    id_main_df = id_column.str.findall(r'MW\d+|AW\d+')  # will have duplicates
    id_main_df = id_main_df.apply(lambda x: x[0] if x else np.nan)  # keep the first find

    id_alt_df = id_column.str.findall(
        r'(MW\w+)|(AW\w+)|\((\w+)\)|\[(\w*?)\]|^(\w+) | (\w+)$| (\w+) |(\w+),')
    id_alt_df = id_alt_df.apply(lambda x: [j for i in x for j in i if j != ''])

    # remember to drop len(id)<4, when use, useless
    # exclude ext(ra) amp(lification) pur(ify)
    id_alt_df = id_alt_df.apply(
        lambda x: [j for j in x if re.match(r'ext|amp|pur', j) is None])

    return pd.DataFrame({"id_main": id_main_df, "id_alt": id_alt_df})


def clean_batch(df_raw: pd.DataFrame) -> pd.Series:
    """ clean batch information from sample_df

    Args:
        df_raw: raw sample df
    
    Returns:
        pd.Series with [batch1_dict, batch2_dict, ...]
        dict{'name': str; 'idx': int; 'pool': bool}
    """

    batch_columns = np.array([
        "Contract ID",
        "Contract ID for sequencing"])
    index_columns = np.array([
        "Forward Index",
        "Reverse Index"])
    pool_columns = np.array([
        "Pool name for screening (if necessary)",
        "Pool name for Sequencing"])

    # batch nan means missing, else [batches]
    batch_pd = df_raw[batch_columns].copy()
    batch_pd = batch_pd.stack(dropna=False).str.findall(r'F[A-Z0-9]*')
    batch_pd = batch_pd.groupby(level=0).sum()

    # index 0 means missing, else [index]
    index_pd = df_raw[index_columns].copy()
    index_pd = index_pd.stack(dropna=False).str.findall(r'\d+')
    index_pd = index_pd.groupby(level=0).sum()

    # pool [] means missing, else [pools]
    pool_pd = df_raw.copy()
    pool_pd["idx"] = pool_pd.index
    pool_pd = pool_pd[np.append(pool_columns, "idx")]
    pool_pd = pool_pd.melt(id_vars=["idx"]).drop_duplicates(subset=['idx','value'])
    pool_pd = pool_pd.pivot(index="idx", columns="variable", values="value")
    pool_pd = pool_pd.applymap(lambda x: [x] if pd.notnull(x) else []).sum(axis=1)

    # return dict Series
    df_batch = pd.DataFrame({"idx":batch_pd.copy().index})
    df_batch["batch"] = df_batch["idx"].apply(
    lambda x: {
        "names" : batch_pd.loc[x],
        "idx" : index_pd.loc[x],
        "pool" : pool_pd.loc[x]}
    )
    df_batch= df_batch.drop(labels=['idx'], axis=1)

    return df_batch["batch"]


def __check_ref(x: str) -> str:
    """ check reference for bam"""
    if re.findall(r'canfam', x, re.IGNORECASE):
        return "CamFam3.1"
    elif re.findall(r'L.D', x, re.IGNORECASE):
        return  "L.Dalen_14"
    else:
        return np.nan


def __find_bam(id_main: str, df_bam: pd.DataFrame) -> dict:
    """ find bam for sample """
    tmp_df = df_bam[df_bam["id_main"].isin([id_main])].sort_values(by="size")
    tmp_df.drop_duplicates(subset="ref", keep="last", inplace=True)
    return pd.Series(tmp_df.full_path.values, index=tmp_df.ref).to_dict()


def find_fq():
    

    pass

def read_sample(csv: list[str]) -> pd.DataFrame:
    """ Read csv(list) into sample_df

    csv contain information of AW, MW samples

    Args:
        csv: list of csv files

    Returns:
        None
    """
    list_df = []
    for file_csv in csv:
        tmp_df = pd.read_csv(file_csv)
        print(f'File read: {file_csv} with {tmp_df.shape[0]} rows')

        print(f'Remove {sum(pd.isna(tmp_df["Lab ID"]))} rows with missing Lab ID')
        tmp_df.dropna(subset=["Lab ID"], inplace=True)

        print(f'Extract {tmp_df.shape[0]} ID from Lab ID')
        id_clean_pd = clear_id(tmp_df["Lab ID"])

        # plug main into the raw df
        tmp_df["id_main"] = id_clean_pd["id_main"]
        tmp_df.dropna(subset=["id_main"], inplace=True)  # drop na in main id

        # drop na and duplicates
        id_clean_pd.dropna(subset=["id_main"], inplace=True)
        id_clean_pd.drop_duplicates("id_main", inplace=True)
        n_unique_sample = id_clean_pd.shape[0]
        print(f'A total of {n_unique_sample} unique samples in {file_csv}')

        # generate a tmp sample_df
        out_df = sample_df(n_unique_sample)
        # input main id and alt
        out_df.loc[:, "id"] = id_clean_pd.loc[:, "id_main"].values  # no index
        out_df.loc[:, "id_alt"] = id_clean_pd.loc[:, "id_alt"].values

        # clean DP in raw
        DP_columns = np.array(
            [
                "# Estimated coverage from unique hits to CanFam31",
                "Estimated coverage (Canfam31)"
            ]
        )
        DP_str = DP_columns[np.isin(DP_columns, tmp_df.columns.values)][0]
        tmp_df.rename(columns={DP_str: "DP"}, inplace=True)  # rename to DP

        tmp_df.loc[
            ~tmp_df["DP"].astype(str).str.fullmatch(r'[0-9,.]*', na=False),
            "DP"
        ] = np.nan

        tmp_df["DP"] = tmp_df["DP"].apply(
            lambda x: x if type(x) == float else float(x.replace(",", "."))
        )
        # read DP
        out_df.loc[:, "dp"] = out_df["id"].apply(
            lambda x: np.nanmax(
                tmp_df.loc[
                    tmp_df["id_main"].isin([x]),
                    "DP"
                ].tolist()
            )
        )

        # TODO: read batch
        # clean batch in raw
        # one sample could have two line with two library
        tmp_df["batch_dict"] = clean_batch(tmp_df)
        out_df.loc[:, "batch"] = out_df["id"].apply(
            lambda x: 
                tmp_df.loc[
                    tmp_df["id_main"].isin([x]),
                    "batch_dict"
                ].tolist()
        )

        # add to return list
        list_df.append(out_df)

    # concate all df and return
    return pd.concat(list_df)


def walk_dir_fq(in_dir: str, list_of_files: list, sftp: paramiko.sftp) -> list:
    """ Walk through a sftp directory and search for FASTQ

    Args:
        in_dir: paramiko sftp dir
        list_of_files: a list to store file information
        sftp: sftp session

    Return:
        list with [[filename, fullpath]...]
    """
    for files_att in sftp.listdir_attr(in_dir):
        if S_ISDIR(files_att.st_mode):
            # file is directory
            walk_dir_fq(in_dir + "/" + files_att.filename, list_of_files, sftp)

        if S_ISREG(files_att.st_mode):
            # file is regular file
            # check fq.gz, fastq.gz
            file_name = files_att.filename
            if file_name.endswith(tuple([".fq.gz", ".fastq.gz"])):
                full_path = in_dir + "/" + file_name
                list_of_files.append([file_name, full_path])

    return 0


def walk_dir_bam(in_dir: str, list_of_files: list, sftp: paramiko.sftp) -> list:
    """ Walk through a sftp directory and search for BAMs

    Args:
        in_dir: paramiko sftp dir
        list_of_files: a list to store file information
        sftp: sftp session

    Return:
        list with [[filename, fullpath, date, size]...]
    """
    for file_att in sftp.listdir_attr(in_dir):
        if S_ISDIR(file_att.st_mode):
            # file is directory
            walk_dir_bam(in_dir + "/" + file_att.filename, list_of_files, sftp)

        if S_ISREG(file_att.st_mode):
            # file is regular file
            # check fq.gz, fastq.gz
            file_name = file_att.filename
            if file_name.endswith(tuple([".bam"])):
                full_path = in_dir + "/" + file_name
                file_date = datetime.fromtimestamp(file_att.st_mtime)
                file_size = file_att.st_size
                list_of_files.append(
                    [file_name, full_path, file_date, file_size]
                )

    return 0


def get_all_fq(user: str, passwd: str) -> pd.DataFrame:
    """ Scan all fastq files in ERDA session

    Args:
        user: user_name for sftp connection
        passwd: password for sftp connection

    Returns:
        bam data frame
    """
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    try:
        # Connect to the SSH server
        ssh.connect("io.erda.dk", username=user, password=passwd)

        # Create SFTP client
        sftp = ssh.open_sftp()

        # Walt through directory
        list_of_file = []
        dir_full = sftp.listdir()
        for dir_f in [d for d in dir_full if d.startswith("F")]:
            walk_dir_fq("." + "/" + dir_f, list_of_file, sftp)

        # Close the SFTP connection
        sftp.close()

    finally:
        # Close the SSH connection
        ssh.close()

    df_fq = pd.DataFrame(np.array(list_of_file),
                         columns=["file_name", "full_path"])

    # save to csv
    df_fq.to_csv("fq_erda.csv")

    return df_fq


def get_all_bam(user: str, passwd: str, dir_list: list[str]) -> pd.DataFrame:
    """ Scan all bam files in ERDA session

    Args:
        user: user_name for sftp connection
        passwd: password for sftp connection
        dir_list: list of directory to search

    Returns:
        bam data frame
    """
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    try:
        # Connect to the SSH server
        ssh.connect("mjolnirhead01fl.unicph.domain",
                    username=user,
                    password=passwd)

        # Create SFTP client
        sftp = ssh.open_sftp()

        # Walt through directory
        list_of_bam = []
        for dir_sub in dir_list:
            walk_dir_bam(dir_sub, list_of_bam, sftp)

        # Close the SFTP connection
        sftp.close()

    finally:
        # Close the SSH connection
        ssh.close()

    df_bam = pd.DataFrame(
        np.array(list_of_bam),
        columns=["file_name", "full_path", "date", "size"]
    )

    # save to csv
    df_bam.to_csv("bam_server.csv")

    return df_bam


def arg():
    """ argument management """

    parse = argparse.ArgumentParser(prog='database for canid data')

    parse.add_argument(
        '--in_sample_file',
        help='file with list of sample files, should be MW and AW list',
        type=str,
        default=None)
    parse.add_argument(
        '--in_sample_csv',
        help='csv of sample meta, preprocessed',
        type=str,
        default=None)

    parse.add_argument(
        '--in_bam_dir_file',
        help='file with list of bam directories',
        type=str,
        default=None)
    parse.add_argument(
        '--in_bam_csv',
        help='csv of bam meta, preprocessed',
        type=str,
        default=None)
    parse.add_argument(
        '--server_usr',
        help='usr_name to connect server',
        type=str,
        default=None)
    parse.add_argument(
        '--server_pw',
        help='password to connect server',
        type=str,
        default=None)

    parse.add_argument(
        '--in_fq_dir_file',
        help='file with list of fq directories',
        type=str,
        default=None)
    parse.add_argument(
        '--in_fq_csv',
        help='csv of fq meta, preprocessed',
        type=str,
        default=None)
    parse.add_argument(
        '--erda_usr',
        help='usr_name to connect erda',
        type=str,
        default=None)
    parse.add_argument(
        '--erda_pw',
        help='password to connect erda',
        type=str,
        default=None)

    return parse.parse_args()


def main():
    # read args
    args = arg()

    # read sample df
    if args.in_sample_csv is None:
        list_of_sample_files = []
        with open(args.in_sample_file) as fh:
            list_of_sample_files.append(fh.readline().strip())
        df_sample = read_sample(list_of_sample_files)
    else:
        df_sample = pd.read_csv(args.in_sample_csv)

    # read bam df
    if args.in_bam_csv is None:
        list_of_bam_dir = []
        with open(args.in_bam_dir_file) as fh:
            list_of_bam_dir.append(fh.readline().strip())
        df_bam_raw = get_all_bam(
            user=args.server_usr,
            passwd=args.server_pw,
            dir_list=list_of_bam_dir
        )
    else:
        df_bam_raw = pd.read_csv(args.in_bam_csv)

    # read raw df
    if args.in_fq_csv is None:
        df_fq_raw = get_all_fq(
            user=args.erda_usr,
            passwd=args.erda_pw
        )
    else:
        df_fq_raw = pd.read_csv(args.in_fq_csv)

    # match bams for sample
    # clean df_bam
    df_bam_raw["ref"] = df_bam_raw["file_name"].apply(__check_ref)
    df_bam_raw["is_merge"] = df_bam_raw["file_name"].apply(
        lambda x: True if re.findall(r'merge', x, re.IGNORECASE) else False
    )
    # exclude some folders, "from_pub", "temp"
    df_bam_raw = df_bam_raw[df_bam_raw.full_path.apply(
        lambda x: False if re.findall(r'(from_pub)|(temp)', x, re.IGNORECASE) else True
    ).values]
    # clean bam id
    df_bam_raw["bam_id_clean"] = df_bam_raw["file_name"].str.replace(
        r'(\.canfam.*)|(\.L\.D.*)|(\.realign.*)|(\.bam.*)|(merg[a-zA-Z]*_*)',
        "",
        flags=re.IGNORECASE)
    df_bam_raw["id_main"] = clear_id(df_bam_raw["bam_id_clean"])["id_main"]


    # match raw for sample

    # save_csv for sample


if __name__ == '__main__':
    main()
