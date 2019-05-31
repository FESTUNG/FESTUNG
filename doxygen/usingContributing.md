# Using and contributing to FESTUNG

## How to obtain FESTUNG

Clone the Git repository:

```bash
git clone https://github.com/FESTUNG/FESTUNG.git
```

## How to update FESTUNG

If you want to obtain the latest updates on an existing copy of FESTUNG, simply update the repository and pull in all changes:

```bash
git pull
```

## Working on FESTUNG

If you want to modify or add files in one of the problem solvers or the core, it is best practice to create a fork of our repository in your own account, on which you perform your changes.
That way you can use Git to keep track of your changes and push them to GitHub.

### Creating a fork

To create a fork, simply click "Fork" in the top right corner on the GitHub repository website.

### Clone your own fork

To clone your own fork simply run
```bash
git clone https://github.com/<your name>/FESTUNG.git
```
where you replace `<your name>` by your GitHub login.

### Pull in changes from our GitHub repository

Note that your forked repository won't receive any updates to FESTUNG. To pull in recent changes that happened in our repository, simply run:

```bash
git pull https://github.com/FESTUNG/FESTUNG.git master
```

### Implementing a new problem solver

Create a copy of the `template` folder to obtain a clean problem solver with all required functions:

```bash
cp -r template <problem name>
```

### Contributing changes

If you have developed something new and think this is useful also for others, push your changes to your forked repository and open a pull request. Please make sure that you follow the [Naming convention](namingConvention.md) and that all routines are documented using the doxygen syntax.

## Need help or found a bug?

If you need help or found a bug, please open an issue in the GitHub bug tracker at https://github.com/FESTUNG/FESTUNG/issues and describe your question or problem as detailed as possible.
We are grateful for any feedback.

