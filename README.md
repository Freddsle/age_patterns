# age_patterns




# Install and run with poetry
```console
# install poetry
# for details look for https://python-poetry.org/docs/
sudo apt-get install curl
curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python3 -

# poetry will be accessible in current session
source $HOME/.poetry/env

# prepare project
git clone https://github.com/Freddsle/age_patterns
cd ./code


# create env
poetry env use python3.10
poetry install

# Run
poetry run python code/file.py

# or for run .ipynb files
poetry run jupyter notebook
```

## Used Python and libraries:
- Python 3.10.1
- NetworkX 2.8
- pandas 1.4.2
- community 1.0.0b1
- Matplotlib 3.5.2

- Biopython 1.79
- seaborn 0.11.2
