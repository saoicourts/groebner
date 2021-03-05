Some practice for coding related to our reading of *Ideals, Varieties, and Algorithms*.

## Suggested Tooling Setup for Contributing
There are a few people in this group that haven't done this kind of stuff before so I will be more verbose than usual here. Here are recommendations to get you using the same setup I am:
- If you are using an up-to-date version of Windows, you can install [Windows Subsystem for Linux](https://docs.microsoft.com/en-us/windows/wsl/install-win10) which is basically a lightweight virtual machine onto which you can install a copy of linux. It gives you a quick way to get access to a unix command line without having to change your operating system (many tools work more naturally with a unix command line than with Windows although this is certainly not a requirement).

**One note for what follows**: I assume you have access to a unix-style command line for the commands I give you. If you choose not to go that route, you may need to change them. I tried to include links where I could that would explain how to do that.
- Anaconda is great for managing python environments but is optional.
    - If you go this route, you can start by [installing Anaconda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
    - Once it is installed, you can set up a new environment by issuing the command
    `conda create -n <your_choice_of_name> python=3.8` to get everything ready to go.
    - Issue the command `conda activate <your_choice_of_name>` whenever you want to 
    work on the project. Anything you install will be contained within that
    environment and won't affect other things you work on.
- If you haven't installed Anaconda, make sure to install python 3.8 (using anything >= 3.6 should work fine).
    - You have to make sure you have [pip](https://pip.pypa.io/en/stable/) installed (for installing additional python packages as needed) although it is automatically packaged with python. If you are using Anaconda (below) it is also already installed.
- Once you have python and `pip` installed, you can install the following by typing (e.g.) `pip install jupyterlab numpy [...]`. Some recommended packages are:
    - `jupyterlab` (for running jupyter notebooks)
    - `numpy` (numerical package with a ton of useful functions)
    - `pytest` is the testing package we will use to build unit tests, but it should install automatically later
- [Virtual Studio Code](https://code.visualstudio.com/) is a really nice cross-platform IDE. It has a ton of extensions for every language and technology you could possibly want and it makes it a one-stop shop for development.

## Downloading This Code
We are using `git` to collaborate over the internet. If you haven't seen it before, you can read more about it in [Atlassian's wonderful documentation for it](https://www.atlassian.com/git). There is even a GUI interface for it these days, although I haven't gotten around to trying it yet.

Once you have `git` installed, you can clone this repository (copy the code to your local computer) by issuing the command `git clone https://github.com/NicoCourts/groebner` in the place you want it to appear. After it downloads, you will now have a folder named "groebner" (which you are free to rename) that contains all the files in this repo.

## Getting Ready to Run
Once you have your tools in place and the code copied to your computer (let's say it's in the folder */home/nico/groebner*), you can do the following

1. Navigate to the folder */home/nico* in your terminal
1. Install the package in "editing mode" (keeps the package up-to-date as you are coding) by typing `pip install -e groebner`
1. You should be ready to go! 
    - You can test this by typing `python` and then typing `import groebner` into the prompt. There should be no errors. Typing `groebner.QQ` should print the string `the rational numbers`. Type `exit()` to leave the REPL.

## Running Example Notebooks
I have started a few Jupyter notebooks in the [notebooks directory](https://github.com/NicoCourts/groebner/tree/master/notebooks/). To run them, navigate to */home/nico/groebner* (or any directory above it) and type `jupyter lab`. You should either see a link to follow or a browser will automatically open with the interface. From there you can choose which notebook to run and execute the code yourself.

## Running Tests
Unit tests will go in the [tests directory](https://github.com/NicoCourts/groebner/tree/master/tests/). To run them, navigate to anywhere inside the */home/nico/groebner* directory and type `pytest`. You should see a concise output with the status of all tests (there are just 2 at the moment of this writing).

## Questions/Comments
These are always welcome. Shoot me an email or pop on the discord any time to talk.