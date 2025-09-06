import paramiko
import os
import re

class ssh():
    
    def __init__(self, hostname=None, port=None, username=None, password=None, mode='direct',config_file=None,target_pathname='cnspy_ssh'):
        """
        Initialize the SSH connection parameters.
        
        :param hostname: The hostname of the remote server.
        :param port: The port number of the remote server.
        :param username: The username to log in to the remote server.
        :param password: The password to log in to the remote server.
        :param mode: The mode of the SSH connection, 'direct' or 'proxy'.
        """
        # Set the default values for the connection parameters
        # if they are not provided, read them from the config file
        self.is_ssh_class = True
        self.windowspathdelimiter = '\\'
        self.unixpathdelimiter = '/'
        if config_file is None:
            self.config_file = os.path.join(os.getcwd(),'config').replace(self.windowspathdelimiter,self.unixpathdelimiter)
        else:
            self.config_file = config_file
        self.config = self.parse_ssh_config()
        self.mode = mode
        if self.mode == 'direct':
            if hostname is None:
                target_host = self.config['target_host']
                self.target_hostname = target_host['HostName']
                self.target_port = target_host['Port']
                self.target_username = target_host['User']
                self.target_password = target_host['Password']
            else:
                self.target_hostname = hostname
                self.target_port = port
                self.target_username = username
                self.target_password = password
        elif self.mode == 'proxy':
            target_host = self.config['target_host']
            self.target_hostname = target_host['HostName']
            self.target_port = target_host['Port']
            self.target_username = target_host['User']
            self.target_password = target_host['Password']
            proxy_host = self.config['proxy_host']
            self.proxy_hostname = proxy_host['HostName']
            self.proxy_port = proxy_host['Port']
            self.proxy_username = proxy_host['User']
            self.proxy_password = proxy_host['Password']
            self.proxy_ssh = paramiko.SSHClient()
            self.proxy_ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            self.proxy_command = f"ssh -q -W {self.target_hostname}:{self.target_port} {self.proxy_username}@{self.proxy_hostname}"
        else:
            raise ValueError("Invalid mode. Choose 'direct' or 'proxy'.")


        self.sftp = None
        self.target_pathname = target_pathname

    def parse_ssh_config(self):
        config = {}
        with open(self.config_file, 'r') as file:
            lines = file.readlines()
            current_host = None
            for line in lines:
                line = line.strip()
                if line.startswith('Host '):
                    current_host = line.split(' ')[1]
                    config[current_host] = {}
                elif line == '':
                    pass
                else:
                    key, value = line.split(' ')
                    config[current_host][key] = value
        return config

    def upload_directory(self, local_path, remote_path):
        """
        Recursively upload a directory to a remote server via SFTP.

        :param sftp: An active SFTP session.
        :param local_path: The local directory to upload.
        :param remote_path: The target directory on the remote server.
        """
        for root, dirs, files in os.walk(local_path):
            # Construct the corresponding remote directory path
            remote_dir = os.path.join(remote_path, os.path.relpath(root, local_path)).replace(self.windowspathdelimiter,self.unixpathdelimiter)
            
            try:
                self.sftp.mkdir(remote_dir)
            except IOError:  # Directory already exists
                pass
            
            for file in files:
                local_file = os.path.join(root, file).replace(self.windowspathdelimiter,self.unixpathdelimiter)
                remote_file = os.path.join(remote_dir, file).replace(self.windowspathdelimiter,self.unixpathdelimiter)
                print(f"Uploading {local_file} to {remote_file}")
                self.sftp.put(local_file, remote_file)
    
    def upload_file(self, local_file, remote_file):
        """
        Upload a file to a remote server via SFTP.

        :param sftp: An active SFTP session.
        :param local_file: The local file to upload.
        :param remote_file: The target file on the remote server.
        """
        print(f"Uploading {local_file} to {remote_file}")
        self.sftp.put(local_file, remote_file)
    
    def download_directory(self, remote_path, local_path):
        """
        Recursively download a directory from a remote server via SFTP.

        :param sftp: An active SFTP session.
        :param remote_path: The remote directory to download.
        :param local_path: The target directory on the local machine.
        """
        for root, dirs, files in self.sftp.walk(remote_path):
            # Construct the corresponding local directory path
            local_dir = os.path.join(local_path, os.path.relpath(root, remote_path)).replace(self.windowspathdelimiter,self.unixpathdelimiter)
            os.makedirs(local_dir, exist_ok=True)
            
            for file in files:
                remote_file = os.path.join(root, file).replace(self.windowspathdelimiter,self.unixpathdelimiter)
                local_file = os.path.join(local_dir, file).replace(self.windowspathdelimiter,self.unixpathdelimiter)
                print(f"Downloading {remote_file} to {local_file}")
                self.sftp.get(remote_file, local_file)

    def download_file(self, local_file, remote_file):
        """
        Download a file to a remote server via SFTP.

        :param sftp: An active SFTP session.
        :param local_file: The local file to upload.
        :param remote_file: The target file on the remote server.
        """
        print(f"Downloading {remote_file} to {local_file}")
        self.sftp.get(local_file, remote_file)      
                
    def connect(self):
        """
        Connect to the remote server via SSH and open an SFTP session.
        """
        if self.mode == 'direct':
            self.ssh = paramiko.SSHClient()
            self.ssh.load_system_host_keys()
            self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            self.ssh.connect(self.target_hostname, self.target_port, self.target_username, self.target_password)
            print(f"Connected to the {self.target_hostname} server")
        elif self.mode == 'proxy':
            # Connect to the proxy server
            self.proxy_ssh.connect(self.proxy_hostname, port=int(self.proxy_port), 
                            username=self.proxy_username, password=self.proxy_password)
            print("Connected to the proxy server")
            
            # Open a transport to the remote server through the proxy
            transport = self.proxy_ssh.get_transport()
            dest_addr = (self.target_hostname, int(self.target_port))
            local_addr = ('', 0)
            channel = transport.open_channel("direct-tcpip", dest_addr, local_addr)
            
            # Create an SSH client for the remote connection
            self.ssh = paramiko.SSHClient()
            self.ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            self.ssh.connect(self.target_hostname, port=int(self.target_port), 
                            username=self.target_username, password=self.target_password, 
                            sock=channel)
            print("Connected to the remote server through the proxy")
            #self.ssh.connect(self.target_hostname, self.target_port, self.target_username,  self.target_password, sock=proxy)
            #print("Connected to the remote server through the proxy")
        self.sftp = self.ssh.open_sftp()
        
    def close(self):
        """
        Close the SFTP session and the SSH connection.        
        """
        self.sftp.close()
        self.ssh.close()
        if self.mode == 'proxy':
            self.proxy_ssh.close()
        print("Connection closed")
        
    def execute(self, remote_command):
        """
        Execute a command on the remote server via SSH.
        
        :param remote_command: The command to execute on the remote server.
        """
        print(f"Executing command: {remote_command}")
        stdin, stdout, stderr = self.ssh.exec_command(remote_command)
        print("Command output:")
        print(stdout.read().decode())
        print("Command error (if any):")
        print(stderr.read().decode())
    
    def home_directory(self):
        """
        Return the home directory of the remote server.
        """
        stdin, stdout, stderr = self.ssh.exec_command('echo $HOME')
        return stdout.read().decode().strip()
    
    def working_directory(self):
        """
        Return the working directory of the remote server.
        """
        return os.path.join(self.home_directory(),self.target_pathname).replace(self.windowspathdelimiter,self.unixpathdelimiter)

    def expand_nodelist(self,nodelist):
        """
        Expands a Slurm nodelist string into a list of individual node names.
        """
        match = re.match(r'(\D+)\[(.+)\]', nodelist)
        if match:
            prefix = match.group(1)
            ranges = match.group(2).split(',')
            nodes = []
            for item in ranges:
                if '-' in item:
                    start, end = map(int, item.split('-'))
                    nodes.extend([f"{prefix}{i}" for i in range(start, end + 1)])
                else:
                    nodes.append(f"{prefix}{item}")
            return nodes
        else:
            return [nodelist]

    def find_idle_node_list(self,partition_name=None):
        """
        Return idle node list on the cluster.
        """
        if partition_name is None:
            stdin, stdout, stderr = self.ssh.exec_command('sinfo -t idle -o "%N"')
        else:
            stdin, stdout, stderr = self.ssh.exec_command(f'sinfo -p {partition_name} -t idle -o "%N"')
        return self.expand_nodelist(stdout.read().decode().splitlines()[1])
    
    def idle_node(self,n,partition_name=None):
        """
        Return n idle node string with form cu1,cu2,... .
        """
        idle_node_list = self.find_idle_node_list(partition_name=None)
        return ",".join(idle_node_list[0:n])
    
    def node_status(self,nodename):
        """
        Return the status of the specific node.
        """
        stdin, stdout, stderr = self.ssh.exec_command(f'sinfo -n {nodename} -o "%t"')
        return stdout.read().decode().splitlines()[1]
        
    

class slurm_script():
    
    def __init__(self, node=None, qos=None, np=1, jobname=None,partition_name=None):
        self.node = node
        self.qos = qos
        self.np = np
        self.jobname = jobname
        self.partition_name = partition_name
        if self.node is None:
            raise ValueError("Node name must be provided")
        if self.qos is None:
            raise ValueError("QoS name must be provided")
        if self.np is None:
            raise ValueError("Number of processors must be provided")
        if self.jobname is None:
            raise ValueError("Job name must be provided")
        self.slurmname = f"script_{self.np}_{self.qos}.slurm"
    
    def submit_command(self):
        """
        Submit a job script to the batch system.
        """
        return f"sbatch {self.slurmname}"
        
    def write_script(self,path,mlx4_0=False):
        """
        Write a job script to a file.
        
        :param script: The content of the job script.
        """
        WINDOWS_LINE_ENDING = b'\r\n'
        UNIX_LINE_ENDING = b'\n'
        
        with open(os.path.join(path,self.slurmname), 'w') as file:
            file.write("#!/bin/bash\n\n")
            file.write(f"#SBATCH --qos={self.qos} -w {self.node} -n {self.np}\n\n")
            if self.partition_name is not None:
                file.write(f"#SBATCH --partition={self.partition_name}\n\n")
            if mlx4_0:
                file.write("export OMPI_MCA_btl=^openib\n\n")
            if self.np > 1:
                file.write(f"mpirun -n $SLURM_NPROCS {self.jobname} &> {self.jobname}_out.dat\n")
            else:
                file.write(f"./{self.jobname} &> {self.jobname}_out.dat\n")
        
        with open(os.path.join(path,self.slurmname), 'rb') as open_file:
            content = open_file.read()
            
        # Windows ➡ Unix
        content = content.replace(WINDOWS_LINE_ENDING, UNIX_LINE_ENDING)

        with open(os.path.join(path,self.slurmname), 'wb') as open_file:
            open_file.write(content)

        print(f"Job script written to {path}")
    
    def clear_script(self,path):
        """
        Clear the job script file.
        """
        for item in os.listdir(path):
            if item.endswith(".slurm"):
                os.remove(os.path.join(path,item))
        print(f"All job scripts are removed")