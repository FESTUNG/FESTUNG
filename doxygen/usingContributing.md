# Using and contributing to FESTUNG

The FESTUNG project is split into several repositories:

* The core: it contains grid data structures, basis functions, quadrature rules, integration and assembly routines, visualization functions etc. Every problem needs this.
* Problem solvers: it contains the steps of a solver in the generic problem framework for a certain problem and makes use of routines from the core. You will need at least one (where you implement your algorithm).
* The umbrella project: it contains the main driver routine and pulls in the core and problem solvers as submodules. Cloning this is the easiest to obtain all that is needed

## How to obtain FESTUNG

Clone the umbrella project with all submodules:

```bash
git clone --recursive https://github.com/FESTUNG/project
```

## How to update FESTUNG

If you want to obtain the latest updates on an existing copy of FESTUNG, simply update the umbrella project and pull in all changes in the submodules:

```bash
git pull --recurse-submodules
git submodule update
```

## Working on the core

If you want to modify or add files to the core, it is best practice to create your own branch on which you perform your changes.
That way you can pull in changes from the GitHub repository and selectively merge them into your version.

### Check out your own branch

To checkout your own branch, replace `<branch name>` in the following command by your desired name:
```bash
cd core
git checkout -b <branch name>
```

### Pull in changes from the GitHub repository

The following command pulls in changes from the GitHub repository and tries to apply your changes on top:
```bash
cd core
git pull --rebase origin master
```

## Working on an existing problem

If you want to modify or add files to one of the problem solvers in FESTUNG, it is best practice to create your own branch on which you perform your changes.
That way you can pull in changes from the GitHub repository and selectively merge them into your version.

### Check out your own branch

To checkout your own branch, replace `<problem name>` and `<branch name>` in the following command by your desired names:
```bash
cd <problem name>
git checkout -b <branch name>
```

### Pull in changes from the GitHub repository

The following command pulls in changes from the GitHub repository and tries to apply your changes on top:
```bash
cd <problem name>
git pull --rebase origin master
```

## Implementing a new problem solver

If you want to implement a new solver it is highly recommended to use a repository to keep track of changes.
For that, simply create a GitHub repository in your own account, possibly by forking from https://github.com/FESTUNG/template, and add it as a new submodule to your project folder:

```bash
git submodule add https://github.com/<your name>/<your repository> <problem name>
```

Then continue working on it as described above.

## Need help or found a bug?

If you need help or found a bug, please open an issue in the GitHub bug tracker at https://github.com/FESTUNG/project/issues and describe your question or problem as detailed as possible.
We are grateful for any feedback.

