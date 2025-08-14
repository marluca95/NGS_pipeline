import argparse
import os
from pybis import Openbis
import sys
import getpass


try:
    # Parse arguments
    parser = argparse.ArgumentParser(description='Script to let openBIS users to download data from the new openbis ELN(April 24th 2023)')
    parser.add_argument('--path_to_output', help='Path to the output folder', required=True)
    parser.add_argument('--pool_id', help='ID of the pool in openBIS(e.g. /GFB/POOL-1)', required=True)
    parser.add_argument('--flowcell_id',help='ID of the flowcell. You can find this value in the notification of NGSvivo(e.g. Flowcell: 000000000-DM2VM). It is not mandatory. To use only in case the user wants to download sample data of a single flowcell.',
                        required=False)
    parser.add_argument('--token',help='For no ETH users to login to openBIS is mandatory to follow these instruction to get the token: ',
                        required=False)
    args = parser.parse_args()
    pool_id = args.pool_id
    path_to_output = args.path_to_output
    flowcell_id = args.flowcell_id
    token = args.token

    # Connect to openbis
    url = "https://openbis-gfb.ethz.ch"
    print("Connecting to openbis...")
    o = Openbis(url, verify_certificates=False)

    if token:
        print("Login via token")
        o.set_token(token)
    else:
        print("Login using credentials")
        # TO CHANGE before running the script
        username = "rakuhn" #input("Enter your openBIS username: ")
        password = "PASSWORD" #getpass.getpass("Enter your openBIS password: ")
        o.login(username, password, save_token=True)

    # If user is connected
    if o.is_session_active():

        # If exist the output folder
        if os.path.isdir(path_to_output):

            print("Connection to openbis with success!")

            # Get pool object
            if o.get_sample(pool_id):

                pool_data = o.get_sample(pool_id)

                # Get the samples of the pool
                for sample in pool_data.children:

                    # Get sample object
                    sample_data = o.get_sample(sample)

                    # Download the datasets
                    datasets = sample_data.get_datasets()

                    # Download all the datasets of the sample
                    for dataset in datasets:

                        os.chdir(path_to_output)

                        # In case the flowcell is set, download only the datasets that belong to the flowcell
                        if flowcell_id:
                            pattern = sample_data.code + "_" + flowcell_id
                            ds = o.get_dataset(dataset.code)
                            sss = ds.file_list

                            # Check if in the sample name there is the flowcell id
                            if any(pattern in s for s in sss):
                                print(
                                    "Downloading of the dataset " + dataset.code + " of the sample " + sample_data.code + "...")
                                if dataset.download():
                                    print("The downloading of the dataset " + dataset.code + " finished with success!")
                                else:
                                    print("ERROR downloading dataset " + dataset.code)
                                    sys.exit()

                        else:
                            # REMEMBER: An object "sample" can have multiple datasets from different flowcells and lanes.
                            if len(dataset.file_list) > 0:

                                print(
                                    "Downloading of the dataset " + dataset.code + " of the sample " + sample_data.code + "...")
                                if dataset.download():
                                    print("The downloading of the dataset " + dataset.code + " finished with success!")
                                else:
                                    print("ERROR downloading dataset " + dataset.code)
                                    sys.exit()
                            else:
                                print(
                                    "The dataset " + dataset.code + " of the sample " + sample_data.code + " doesn't have reads!")
            else:
                print(
                    "ERROR: the pool id " + pool_id + " does not exists. Please double check the ID in the openBIS UI")
                sys.exit()
        else:
            print("ERROR: the folder " + path_to_output + " doesn't exist or you don't have permission on it")
            sys.exit()
    else:
        print(
            "ERROR: Connection to openbis didn't work. If you are an ETH user, replace the variables username and password. If you are a non ETH user, pass the token as argument")
        sys.exit()

    # Close session
    if not token:
        o.logout()

except Exception as e:
    exc_type, exc_obj, exc_tb = sys.exc_info()
    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
    print(exc_type, fname, exc_tb.tb_lineno)
    print(str(e))
