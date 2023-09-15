

import pysftp
import json
from urllib.parse import urlparse
import os


class Sftp:
    def __init__(self, computer, password, port=22):
        """Constructor Method"""
        # Set connection object to None (initial value)
        with open(os.path.join("C:\\Users\\coopy\\.user_config","remote_computer_config.json")) as f:
            config = json.load(f)
        self.connection = None
        self.hostname = config["COMPUTERS"][computer]
        self.username = config["USER"]
        self.password = password
        self.port = port

    def connect(self):
        """Connects to the sftp server and returns the sftp connection object"""

        try:
            # Get the sftp connection object
            self.connection = pysftp.Connection(
                host=self.hostname,
                username=self.username,
                password=self.password,
                port=self.port,
            )
        except Exception as err:
            raise Exception(err)
        finally:
            print(f"Connected to {self.hostname} as {self.username}.")

    def disconnect(self):
        """Closes the sftp connection"""
        self.connection.close()
        print(f"Disconnected from host {self.hostname}")

    def listdir(self, remote_path):
        """lists all the files and directories in the specified path and returns them"""
        for obj in self.connection.listdir(remote_path):
            yield obj

    def listdir_attr(self, remote_path):
        """lists all the files and directories (with their attributes) in the specified path and returns them"""
        for attr in self.connection.listdir_attr(remote_path):
            yield attr

    def download(self, remote_path, target_local_path):
        """
        Downloads the file from remote sftp server to local.
        Also, by default extracts the file to the specified target_local_path
        """

        try:
            print(
                f"downloading from {self.hostname} as {self.username} [(remote path : {remote_path});(local path: {target_local_path})]"
            )

            # Create the target directory if it does not exist
            path, _ = os.path.split(target_local_path)
            if not os.path.isdir(path):
                try:
                    os.makedirs(path)
                except Exception as err:
                    raise Exception(err)

            # Download from remote sftp server to local
            self.connection.get(remote_path, target_local_path)
            print("download completed")

        except Exception as err:
            raise Exception(err)

    def upload(self, source_local_path, remote_path):
        """
        Uploads the source files from local to the sftp server.
        """

        try:
            print(
                f"uploading to {self.hostname} as {self.username} [(remote path: {remote_path});(source local path: {source_local_path})]"
            )

            # Download file from SFTP
            self.connection.put(source_local_path, remote_path)
            print("upload completed")

        except Exception as err:
            raise Exception(err)


if __name__ == "__main__":
    # sftp_url = os.environ.get("SFTPTOGO_URL")

    # if not sftp_url:
    #     print("First, please set environment variable SFTPTOGO_URL and try again.")
    #     exit(0)

    # parsed_url = urlparse(sftp_url)

    # sftp = Sftp(
    #     hostname=parsed_url.hostname,
    #     username=parsed_url.username,
    #     password=parsed_url.password,
    # )

    # Connect to SFTP
    password = input("input password")
    sftp = Sftp("PERL", password)
    sftp.connect()

    # Lists files with attributes of SFTP
    # path = "/"
    # print(f"List of files with attributes at location {path}:")
    # for file in sftp.listdir_attr(path):
    #     print(file.filename, file.st_mode, file.st_size, file.st_atime, file.st_mtime)

    # Upload files to SFTP location from local
    # local_path = "/Users/saggi/Downloads/tls2.png"
    # remote_path = "/tls2.png"
    # sftp.upload(local_path, remote_path)

    # Lists files of SFTP location after upload
    # print(f"List of files at location {path}:")
    # print([f for f in sftp.listdir(path)])

    # Download files from SFTP
    remote_path = "/pscratch/sd/c/cote3804/epoxidation/calcs/molecules/AcetAld/No_bias"
    local_path = "C:\\Users\\coopy\\Downloads"
    sftp.download(
        remote_path, local_path
    )

    # Disconnect from SFTP
    sftp.disconnect()



'''
import json
import os

linux_config_dict = {"USER": "cote3804",
                     "COMPUTERS":{
                         "PERL": "perlmutter-p1.nersc.gov",
                         "ALPINE": "login.rc.colorado.edu"
                     }}

json_object = json.dumps(linux_config_dict)

path = "C:\\Users\\coopy\\.user_config"

with open(os.path.join(path,"remote_computer_config.json"), "w") as f:
    f.write(json_object)
'''