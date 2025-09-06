from setuptools import setup,find_packages
import os
import shutil

# Define the root directory of your project
ROOT_DIR = os.path.abspath(os.path.dirname(__file__))

# Define a function to delete the build and *.egg-info folders
def clean():
    build_dir = os.path.join(ROOT_DIR, 'build')
    egg_info_dir = ROOT_DIR
    for item in os.listdir(ROOT_DIR):
        if item.endswith('.egg-info'):
            egg_info_dir = os.path.join(ROOT_DIR, item)
            break
    if os.path.isdir(build_dir):
        shutil.rmtree(build_dir)
    if os.path.isdir(egg_info_dir):
        shutil.rmtree(egg_info_dir)


setup(
  name = 'cnspy',
  version = '1.0',
  description = 'Clean Numerical Simulation code generater using sympy as symbolic operation tools.',
  author = 'Bo Zhang',
  author_email = 'zilpher@sjtu.edu.cn',
  url = 'https://github.com/sjtu-liao',
  packages = find_packages(),
  include_package_data=True,
  requires = ['numpy','sympy','scipy','paramiko'],
)

clean()








