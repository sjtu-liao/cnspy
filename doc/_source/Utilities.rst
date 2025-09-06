Utilities
=========

Connect to cluster
----------------------

CNSPy provides a simple way to connect to a cluster. There are several options that can be passed to the `ssh` function. The following code snippet shows how to connect to a cluster with a specific username and password:

.. code-block:: python

    from cnspy import cluster
    HN = cluster.ssh(username='username', password='password')
    HN.connect()
    HN.close()

if you write your username and password in a file like::

    Host proxy_host
      HostName frp-boy.top
      User Administrator
      Password 123456
      Port 123456

    Host target_host
      HostName xxx.x.xx.xxx
      User root
      Password 123456
      Port 22

Then you can connect to the cluster by:

.. code-block:: python

    from cnspy import cluster
    HN = cluster.ssh()
    HN.connect(config_file='path/to/your/config/file')
    HN.close()

The `ssh` function will automatically read the configuration file and connect to the target host.

If you want to connect through a proxy, you can use the following code snippet:

.. code-block:: python

    from cnspy import cluster
    HN = cluster.ssh(proxy=True)
    HN.connect( mode='proxy',config_file='path/to/your/config/file')
    HN.close()

The `ssh` function will automatically read the configuration file and connect to the target host through the proxy.

After connecting to the cluster, you can execute commands on the cluster. The following code snippet shows how to execute a command on the cluster:

.. code-block:: python

    HN.execute('ls')

or you can upload/download a file/folder to the cluster:

.. code-block:: python

    HN.upload_file('path/to/local/file', 'path/to/remote/file')
    HN.upload_directory('path/to/local/folder', 'path/to/remote/folder')
    HN.download_file('path/to/local/file', 'path/to/remote/file')
    HN.download_directory('path/to/local/folder', 'path/to/remote/folder')

To find the idle node in the cluster, you can use the following code snippet:

.. code-block:: python

    idle_node_list = HN.find_idle_node_list(partition_name='ifany')

or diectly use 

.. code-block:: python

    idle_node_string = HN.idle_node(2, partition_name='ifany')

which can use in slurm script.

Slurm script
---------------

CNSPy provides a simple way to generate a slurm script. The following code snippet shows how to generate a slurm script:

.. code-block:: python

    from cnspy import cluster
    HN = cluster.ssh(username='username', password='password')
    HN.connect()
    SL = cluster.slurm_script(node=HN.idle_node(2, partition_name='ifany'), qos='super', np=128, jobname = 'jobname')
    SL.write_script(lorenz.ccodepath,mlx4_0=False)

The `mlx4_0` is to solve the porblem we show in troubleshooting.