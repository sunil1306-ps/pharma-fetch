from typing import List, Dict
from Bio import Entrez

def fetch_papers(query: str, max_results: int = 10, debug: bool = False) -> List[Dict]:
    """
    Fetch papers from PubMed based on a query.
    """
    if debug:
        print("Fetching papers for query:", query)

    Entrez.email = "saragadamsunil7@gmail.com"
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()

    if debug:
        print("PubMed IDs found:", record["IdList"])

    paper_ids = record["IdList"]
    papers = []

    for paper_id in paper_ids:
        if debug:
            print("Fetching details for paper ID:", paper_id)

        handle = Entrez.efetch(db="pubmed", id=paper_id, retmode="xml")
        paper_data = Entrez.read(handle)
        handle.close()

        if debug:
            print("Raw paper data for ID", paper_id, ":", paper_data) 

        # Extract relevant fields from the paper data
        try:
            # Navigate through the nested structure
            pubmed_article = paper_data["PubmedArticle"][0]
            medline_citation = pubmed_article["MedlineCitation"]
            article = medline_citation["Article"]

            # Extract title
            title = article.get("ArticleTitle", "No title available")

            # Extract publication date
            publication_date = ""
            article_dates = article.get("ArticleDate", [])
            if article_dates: 
                pub_date = article_dates[0]
                publication_date = f"{pub_date.get('Year', '')}-{pub_date.get('Month', '')}-{pub_date.get('Day', '')}"

            # Extract authors and affiliations
            authors = []
            for author in article.get("AuthorList", []):
                author_name = f"{author.get('LastName', '')}, {author.get('ForeName', '')}"
                affiliations = [affil.get("Affiliation", "") for affil in author.get("AffiliationInfo", [])]
                authors.append({
                    "name": author_name,
                    "affiliations": affiliations
                })

            # Extract corresponding author email (if available)
            corresponding_email = ""
            for author in article.get("AuthorList", []):
                if "Identifier" in author:
                    for identifier in author["Identifier"]:
                        # Handle StringElement type
                        if hasattr(identifier, "attributes") and identifier.attributes.get("Source", "") == "Email":
                            corresponding_email = str(identifier)
                            break

            # Append the extracted data to the papers list
            papers.append({
                "PubmedID": paper_id,
                "Title": title,
                "PublicationDate": publication_date,
                "Authors": authors,
                "CorrespondingAuthorEmail": corresponding_email
            })

        except KeyError as e:
            print(f"Error parsing paper {paper_id}: {e}")
            continue

    return papers