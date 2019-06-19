"""
This file automates the installation of blast software to the
tcrdist/external/blast-x.xx.x/bin/
"""
import os

def install_blast_to_externals(download_from):
    """
    Function takes a download key "ncbi_osx", "ncbi_linux",
    "dropbox_osx", "dropbox_linux"

    curl downloads blast 2.2.16.tar.gz from ncbi or dropbox
    , uncompresses it to tcrdist/externals/ ,
    and then removes the tar.gz file.

    Parameters
    ----------
    download_from : string must be "ncbi_osx", "ncbi_linux", "dropbox_osx", "dropbox_linux"

    Returns
    -------
    installs to tcrdist/external/blast-2.2.16/*

    """
    if not isinstance(download_from, str):
        raise TypeError("the < download_from > arg must be a string")

    if download_from not in ["ncbi_osx","ncbi_linux", "dropbox_osx", "dropbox_linux"]:
        raise KeyError("the < download_from > arg passed to install_blast_to_externals \
        must be one of the following: ncbi_osx, ncbi_linux, dropbox_osx, dropbox_linux")

    address_ncbi_osx = 'https://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.16/blast-2.2.16-universal-macosx.tar.gz'
    address_ncbi_linux = 'https://ftp.ncbi.nlm.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.16/blast-2.2.16-x64-linux.tar.gz'
    address_dropbox_osx = 'https://www.dropbox.com/s/x3e8qs9pk5w6szq/blast-2.2.16-universal-macosx.tar.gz?dl=1'
    address_dropbox_linux = 'https://www.dropbox.com/s/gurbwgcys6xcttm/blast-2.2.16-x64-linux.tar.gz?dl=1'

    path_file = os.path.realpath(__file__)
    path = os.path.dirname(path_file)
    install_path = os.path.join(path, "external")

    def generate_curl(download_link):
        return('curl -o blast.tar.gz {}'.format(download_link))

    curl_url = {"ncbi_osx" : generate_curl(address_ncbi_osx),
    "ncbi_linux" : generate_curl(address_ncbi_linux) ,
    "dropbox_osx" : generate_curl(address_dropbox_osx),
    "dropbox_linux" : generate_curl(address_dropbox_linux)}

    path_file = os.path.realpath(__file__)
    path = os.path.dirname(path_file)
    install_path = os.path.join(path, "external")

    curl_cmd = curl_url[download_from]
    os.system(curl_cmd)

    untar_cmd = "tar -C {} -zxvf blast.tar.gz".format(install_path)
    os.system(untar_cmd)

    rm_cmd = "rm blast.tar.gz"
    os.system(rm_cmd)
