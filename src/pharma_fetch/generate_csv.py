import pandas as pd
from typing import List, Dict

def generate_csv(papers: List[Dict], output_file: str = "output.csv", debug: bool = False):
    """
    Generate a CSV file from the filtered papers.
    """
    data = []
    for paper in papers:
        pubmed_id = paper.get("PubmedID", "")
        title = paper.get("Title", "")
        publication_date = paper.get("PublicationDate", "")
        non_academic_authors = paper.get("NonAcademicAuthors", [])
        company_affiliations = list(set(author["affiliation"] for author in non_academic_authors))
        corresponding_email = paper.get("CorrespondingAuthorEmail", "")

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