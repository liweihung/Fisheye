# Fisheye

## Running the Process
There are two options for running the Fisheye post-process. <br>
Before you do, it is a best practice to pull the latest changes from this repository using `git pull`.

### VS Code
Step 1: Open VS code software.
Step 2: Open the *Fisheye\Scripts* repository folder. <br> 
Step 3: Open the VS Code terminal. From the menu, "View → Terminal". Or use Ctrl + backtick. <br>
Step 4: In the terminal type `conda activate fisheye` and press Enter. <br>
Step 5: In the terminal type `git pull` and press Enter. <br>
Step 6: Open `process_input.py` and enter all the deployment metadata. Save the file. <br>
Step 7: Type `ipython` to activate the IPython terminal. <br>
Step 8: Run `process.py`. <br>

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
