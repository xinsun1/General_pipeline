""" File transfer with paramiko """

import os
import paramiko


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