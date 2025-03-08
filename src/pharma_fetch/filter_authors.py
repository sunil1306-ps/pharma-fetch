from typing import List, Dict

def is_non_academic(affiliation: str) -> bool:
    """
    Check if an affiliation is non-academic (pharma/biotech).
    """
    academic_keywords = ["university", "college", "institute", "lab"]
    return not any(keyword in affiliation.lower() for keyword in academic_keywords)

def filter_authors(paper_data: Dict, debug: bool = False) -> List[Dict]:
    """
    Filter authors with non-academic affiliations.
    """
    non_academic_authors = []
    for author in paper_data.get("Authors", []):
        affiliations = author.get("affiliations", [])
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