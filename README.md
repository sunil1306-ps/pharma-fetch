# Pharma-Fetch

A Python tool to fetch and filter research papers from PubMed with pharmaceutical/biotech affiliations.

---

## Table of Contents

1. [Overview](#overview)
2. [Installation](#installation)
3. [Installation from Test PyPI](#installation-from-test-pypi)
4. [Usage](#usage)
5. [Command-Line Options](#command-line-options)
6. [Output](#output)
7. [Tools and Libraries](#tools-and-libraries)
8. [Non Functional Requirements](#non-functional-requirements)
9. [Contributing](#contributing)
10. [License](#license)

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

## Installation from Test PyPI

To install the `pharma-fetch` package from Test PyPI, follow these steps:

1. **Install the package**:
   Use `pip` to install the package from Test PyPI. Since Test PyPI does not host all dependencies, you need to specify both Test PyPI and the official PyPI as package sources.

   ```bash
   pip install -i https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ pharma-fetch==0.1.1
   ```

2. **Verify the installation**:
   After installation, verify that the `get-papers-list` command is available:

   ```bash
   get-papers-list --help
   ```

   This should display the help message for the command.

3. **Run the tool**:
   Use the `get-papers-list` command to fetch and filter research papers. For example:

   ```bash
   get-papers-list "cancer treatment" -f results.csv
   ```

   This will fetch papers related to "cancer treatment" and save the results to `results.csv`.

**NOTE**:
**Ensure the Script is in Your PATH**
If the script is not found, it means the installation directory is not in your system’s `PATH`. Here’s how to fix it:

#### On Windows:
1. Find the `Scripts` directory where Python installed the script. It’s usually located in:
   ```
   C:\Python312\Scripts\
   ```
   or
   ```
   C:\Users\<YourUsername>\AppData\Local\Programs\Python\Python312\Scripts\
   ```

2. Add this directory to your system’s `PATH`:
   - Open the Start menu and search for "Environment Variables".
   - Click "Edit the system environment variables".
   - In the System Properties window, click "Environment Variables".
   - Under "System variables", find the `Path` variable and click "Edit".
   - Add the `Scripts` directory to the list (e.g., `C:\Python312\Scripts\`).

3. Restart your terminal and try running `get-papers-list` again.


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

- [DeepSeek] (https://chat.deepseek.com/)
- **Python**: The program is written in Python 3.12.
- **Poetry**: Used for dependency management and packaging.
- **Biopython**: Used to interact with the PubMed API.
- **Pandas**: Used to generate the CSV file.
- **argparse**: Used for command-line argument parsing.

---

## Non-Functional Requirements

### Typed Python
- Added type hints to all functions, arguments, and return values using Python’s `typing` module.
- Ensured type consistency across all modules.

### Performance
- Optimized API calls and data processing for efficiency.
- Added a retry mechanism for transient API failures.

### Readability
- Added docstrings, meaningful variable names, and comments for clarity.
- Formatted code consistently for better readability.

### Organization
- Modularized the code into separate files, each handling a specific task.
- Followed the Single Responsibility Principle (SRP).

### Robustness
- Added error handling for invalid queries, API failures, and missing data.
- Used `.get()` with default values to handle missing fields gracefully.

### Relevant Links

- [Biopython Documentation](https://biopython.org/)
- [PubMed API](https://pubmed.ncbi.nlm.nih.gov/)
- [Poetry Documentation](https://python-poetry.org/docs/)
- [DeepSeek] (https://chat.deepseek.com/)


---

## Contact

For questions or feedback, please contact [Saragadam Sunil](mailto:saragadamsunil7@gmail.com).

