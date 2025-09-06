# Clean Numerical Simulation with Python (cnspy)

<!-- PROJECT LOGO -->
<br />

<p align="center">
  <a href="https://https://github.com/sjtu-liao/cnspy">
    <img src="icon/cnspy.png" alt="Logo" width="280">
  </a>
  <p align="center">
    <b>C</b>lean <b>N</b>umerical <b>S</b>imulation with <b>Py</b>thon.
  </p>
  <p align="center">
      <i>Bo Zhang and Shijun Liao</i>
  </p>
  <p align="center">
     <i>School of Ocean and Civil Engineering, Shanghai Jiao Tong University</i>
  </p>
  <p align="center">
     <i>Shanghai 200240, China</i>
  </p>
</p> 

## 📥 Clone or Download

### Clone the Repository

```bash
git clone https://github.com/sjtu-liao/cnspy.git
cd cnspy
```

### Or Download as a ZIP

1. Click the **Code** button on the repository page.
2. Select **Download ZIP**.
3. Extract the contents to your desired location.

---

## 🛠️ Setup Instructions

### System Requirements

Ensure the following dependencies are installed with the minimum versions:

- **GCC**: 4.8.5 or higher
- **OpenMPI**: 2.1.0 or higher
- **GMP**: 6.1.0 or higher
- **MPFR**: 4.2.0 or higher


### 1. Create and Activate a Virtual Environment

```bash
conda create -n <envname>  # Replace <envname> with your environment name
conda activate <envname>
```

### 2. Install Required Packages

```bash
conda install numpy sympy scipy paramiko setuptools ipykernel
```

### 3. Build and Install cnspy

```bash
python setup.py bdist_wheel
pip install dist/cnspy-1.0-py3-none-any.whl
```

### 4. Uninstall cnspy

To remove the package, use:

```bash
pip uninstall cnspy
```

### 5. Reinstall cnspy

To reinstall or update the package, use:

```bash
pip install dist/cnspy-1.0-py3-none-any.whl --force-reinstall
```

### 6. Run an Example

1. Navigate to the `example` folder.
2. Open a `*.ipynb` file using Jupyter Notebook.
3. Run the notebook to test the functionality.

---

## 📚 Documentation and Examples

For detailed documentation and additional examples, please refer to the [documentation page](https://sjtu-liao.github.io/cnspy/).

## 📝 License

This project is licensed under the [MPL-2.0](https://www.mozilla.org/en-US/MPL/2.0/) license. See the LICENSE file for details.

## 🤝 Contributing

Contributions, issues, and feature requests are welcome! Feel free to check the [issues page](https://github.com/sjtu-liao/cnspy/issues) or submit a pull request.

---

Happy using with `cnspy`! 🚀
