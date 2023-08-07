""" File transfer with paramiko """

import os
import paramiko
from stat import S_ISREG, S_ISDIR

def upload(hostname, username, password, list_of_files):
    # Create SSH client
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())

    try:
        # Connect to the SSH server
        ssh.connect(hostname, username=username, password=password)

        # Create SFTP client
        sftp = ssh.open_sftp()

        for i in list_of_files:
            print(f'upload {i}')
            sftp.put(i, "./bam/" + i.split("/")[-1])
            print(f'{i} done')

        # Close the SFTP connection
        sftp.close()

    finally:
        # Close the SSH connection
        ssh.close()

def download(hostname, username, password, list_of_dir, local_dir):
    # Create SSH client
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    try:
        # Connect to the SSH server
        ssh.connect(hostname, username=username, password=password)

        # Create SFTP client
        sftp = ssh.open_sftp()


        for in_dir in list_of_dir:
            # check if local file directory exist
            id = in_dir.split("/")[-1]
            local_dir_id = local_dir + "/" + id
            if not os.path.exists(local_dir_id):
                os.mkdir(local_dir_id)

            for files_att in sftp.listdir_attr(in_dir):
                if S_ISREG(files_att.st_mode):
                    # file is regular file
                    # check fq.gz, fastq.gz
                    file_name = files_att.filename
                    if file_name.endswith(tuple([".fq.gz", ".fastq.gz"])):
                        full_path = in_dir + "/" + file_name
                        local_path = local_dir_id + "/" +file_name

                        print(f'download {local_path}')
                        sftp.get(full_path, local_path)
                        print(f'{local_path} done')
            return 0

        # Close the SFTP connection
        sftp.close()

    finally:
        # Close the SSH connection
        ssh.close()




# Set the SSH credentials and directory path
hostname = "io.erda.dk"
username = ""
password = ""

in_file_of_list_of_files = "list_czh.bam_bai"
with open(in_file_of_list_of_files, 'r') as fh:
    list_files = [l.strip() for l in fh]
upload(hostname=hostname,
       username=username,
       password=password,
       list_of_files=list_files)