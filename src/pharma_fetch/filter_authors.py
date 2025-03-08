from typing import List, Dict, Any

def is_non_academic(affiliation: str) -> bool:
    """
    Check if an affiliation is non-academic (pharma/biotech).

    Args:
        affiliation (str): Affiliation string to check.

    Returns:
        bool: True if the affiliation is non-academic, False otherwise.
    """
    academic_keywords = ["university", "college", "institute", "lab"]
    return not any(keyword in affiliation.lower() for keyword in academic_keywords)

def filter_authors(paper_data: Dict[str, Any], debug: bool = False) -> List[Dict[str, str]]:
    """
    Filter authors with non-academic affiliations.

    Args:
        paper_data (Dict[str, Any]): Paper data containing authors and affiliations.
        debug (bool): Enable debug output (default: False).

    Returns:
        List[Dict[str, str]]: List of non-academic authors with their affiliations.
    """
    non_academic_authors: List[Dict[str, str]] = []
    for author in paper_data.get("Authors", []):
        affiliations: List[str] = author.get("affiliations", [])
        if debug:
            print("Author:", author.get("name", ""))
            print("Affiliations:", affiliations)
        for affiliation in affiliations:
            if is_non_academic(affiliation):
                non_academic_authors.append({
                    "name": author.get("name", ""),
                    "affiliation": affiliation
                })
    if debug:
        print("Non-academic authors:", non_academic_authors)
    return non_academic_authors