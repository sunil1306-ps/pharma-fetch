import pandas as pd
from typing import List, Dict, Any

def generate_csv(papers: List[Dict[str, Any]], output_file: str = "output.csv", debug: bool = False) -> None:
    """
    Generate a CSV file from the filtered papers.

    Args:
        papers (List[Dict[str, Any]]): List of papers with extracted data.
        output_file (str): Output CSV file name (default: output.csv).
        debug (bool): Enable debug output (default: False).
    """
    data: List[Dict[str, str]] = []
    for paper in papers:
        pubmed_id: str = paper.get("PubmedID", "")
        title: str = paper.get("Title", "")
        publication_date: str = paper.get("PublicationDate", "")
        non_academic_authors: List[Dict[str, str]] = paper.get("NonAcademicAuthors", [])
        company_affiliations: List[str] = list(set(author["affiliation"] for author in non_academic_authors))
        corresponding_email: str = paper.get("CorrespondingAuthorEmail", "Not available")

        if debug:
            print("Paper:", pubmed_id, title, publication_date)
            print("Non-academic authors:", non_academic_authors)
            print("Company affiliations:", company_affiliations)
            print("Corresponding email:", corresponding_email)

        data.append({
            "PubmedID": pubmed_id,
            "Title": title,
            "PublicationDate": publication_date,
            "NonAcademicAuthors": "; ".join(author["name"] for author in non_academic_authors),
            "CompanyAffiliations": "; ".join(company_affiliations),
            "CorrespondingAuthorEmail": corresponding_email
        })

    df = pd.DataFrame(data)
    df.to_csv(output_file, index=False)
    if debug:
        print(f"CSV file generated: {output_file}")