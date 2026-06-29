# Fisheye Night Sky Imager Data Processing Pipeline

*Developed and maintained by [Li-Wei Hung](mailto:li-wei_hung@nps.gov)*

> Hung, L.-W., White, J., Joyce, D., Anderson, S. J., & Banet, B. (2024). Fisheye Night Sky Imager: A Calibrated Tool to Measure Night Sky Brightness. *Publications of the Astronomical Society of the Pacific*, 136, 085002. https://doi.org/10.1088/1538-3873/ad6bc1

---

## Background

The NPS Fisheye Night Sky Imager is a camera system developed by the Night Skies Team of the U.S. National Park Service to measure and monitor night sky brightness in national parks. It comprises a Sony IMX455 CMOS sensor housed in a ZWO ASI6200MM camera, a Johnson V filter, and a Sigma 8 mm F3.5 fisheye lens — all commercially available components. The fisheye lens captures the entire sky in a single 30-second exposure. This open-source pipeline processes the resulting images through flat-field correction, astrometric plate solving, photometric calibration using Hipparcos standard stars, positional calibration, median filtering, and final projection in both fisheye and Hammer equal-area views, with a photometric calibration uncertainty of 0.12 mag.

---

## Dependencies

- Python 3.12.12
- Git
- [Astrometry.net](http://nova.astrometry.net/) account and API key

> ⚠️ Never commit your Astrometry.net API key to the repository.

---

## Getting Started

Clone the repository and install dependencies:

```bash
git clone https://github.com/liweihung/Fisheye.git
cd Fisheye/Scripts
conda create --name fisheye python=3.12.12
conda activate fisheye
pip install poetry
poetry install
```

Before each session, pull the latest changes:
```bash
git pull
```

### Contributing
Pull requests are welcome. See GitHub's [Contributing to a project](https://docs.github.com/en/get-started/exploring-projects-on-github/contributing-to-a-project) guide for details.


## Running the Process
There are two options for running the Fisheye post-process. <br>
Before you do, it is a best practice to pull the latest changes from this repository using `git pull`.

### VS Code
Step 1: Open VS code software. <br> 
Step 2: Open the *Fisheye\Scripts* repository folder. <br> 
Step 3: Open the VS Code terminal. From the menu, "View → Terminal". Or use Ctrl + backtick. <br>
Step 4: In the terminal type `conda activate fisheye` and press Enter. <br>
Step 5: In the terminal type `git pull` and press Enter. <br>
Step 6: Open `process_input.py` and enter all the deployment metadata. Save the file. <br>
Step 7: Type `ipython` to activate the IPython terminal. <br>
Step 8: `run process.py`. <br>

### Miniforge
Step 1: In the Fisheye repository scripts folder, open `process_input.py` (with a text editor) and enter all the deployment metadata. Save the file. <br>
Step 2: Open the Miniforge prompt. <br>
Step 3: Activate your Fisheye environment by typing `conda activate fisheye`. <br>
You'll see the environment in parenthesis change to (Fisheye). <br>
Step 4: Use the Windows `cd` command to change directories to *Fisheye\Scripts*. <br>
Step 5: Type `python process.py` and hit enter. <br>

## License

### Public domain

This project is in the worldwide [public domain](LICENSE.md):

> This project is in the public domain within the United States,
> and copyright and related rights in the work worldwide are waived through the
> [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
>
> All contributions to this project will be released under the CC0 dedication.
> By submitting a pull request, you are agreeing to comply with this waiver of copyright interest.
