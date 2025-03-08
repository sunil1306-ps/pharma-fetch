# Pharma-Fetch

A Python tool to fetch and filter research papers from PubMed with pharmaceutical/biotech affiliations.

---

## Table of Contents

1. [Overview](#overview)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Command-Line Options](#command-line-options)
5. [Output](#output)
6. [Tools and Libraries](#tools-and-libraries)
7. [Contributing](#contributing)
8. [License](#license)

---

## Overview

Pharma-Fetch is a command-line tool that:
- Fetches research papers from PubMed based on a user-specified query.
- Filters papers to identify those with at least one author affiliated with a pharmaceutical or biotech company.
- Generates a CSV file with the following columns:
  - `PubmedID`: Unique identifier for the paper.
  - `Title`: Title of the paper.
  - `PublicationDate`: Date the paper was published.
  - `NonAcademicAuthors`: Names of authors affiliated with non-academic institutions.
  - `CompanyAffiliations`: Names of pharmaceutical/biotech companies.
  - `CorrespondingAuthorEmail`: Email address of the corresponding author.

---

## Installation

1. **Install Poetry** (if not already installed):
   ```bash
   pip install poetry
   ```

2. **Clone the Repository**:
   ```bash
   git clone https://github.com/your-username/pharma-fetch.git
   cd pharma-fetch
   ```

3. **Install Dependencies**:
   ```bash
   poetry install
   ```

4. **Activate the Virtual Environment** (optional):
   ```bash
   poetry shell
   ```

---

## Usage

Run the program with a PubMed query:

```bash
poetry run get-papers-list "your query here"
```

### Command-Line Options

- `-h`, `--help`: Show help message and exit.
- `-f FILE`, `--file FILE`: Specify the output CSV file name (default: `output.csv`).
- `-d`, `--debug`: Enable debug output.

### Example Commands

1. Fetch papers with a query and save to the default file (`output.csv`):
   ```bash
   poetry run get-papers-list "cancer treatment"
   ```

2. Fetch papers and save to a custom file:
   ```bash
   poetry run get-papers-list "cancer treatment" -f results.csv
   ```

3. Fetch papers with debug output:
   ```bash
   poetry run get-papers-list "cancer treatment" -d
   ```

4. Display help message:
   ```bash
   poetry run get-papers-list -h
   ```

---

## Output

The program generates a CSV file with the following columns:

| Column Name             | Description                                      |
|-------------------------|--------------------------------------------------|
| `PubmedID`              | Unique identifier for the paper.                 |
| `Title`                 | Title of the paper.                              |
| `PublicationDate`       | Date the paper was published.                    |
| `NonAcademicAuthors`    | Names of authors with non-academic affiliations. |
| `CompanyAffiliations`   | Names of pharmaceutical/biotech companies.       |
| `CorrespondingAuthorEmail` | Email address of the corresponding author.    |

Example CSV output:

```csv
PubmedID,Title,PublicationDate,NonAcademicAuthors,CompanyAffiliations,CorrespondingAuthorEmail
12345678,Example Paper Title,2023-10-01,John Doe; Jane Smith,Pfizer; Novartis,john.doe@example.com
```

---

## Tools and Libraries

- **Python**: The program is written in Python 3.12.
- **Poetry**: Used for dependency management and packaging.
- **Biopython**: Used to interact with the PubMed API.
- **Pandas**: Used to generate the CSV file.
- **argparse**: Used for command-line argument parsing.

### Relevant Links

- [Biopython Documentation](https://biopython.org/)
- [PubMed API](https://pubmed.ncbi.nlm.nih.gov/)
- [Poetry Documentation](https://python-poetry.org/docs/)


---

## Contact

For questions or feedback, please contact [Your Name](mailto:saragadamsunil7@gmail.com).

