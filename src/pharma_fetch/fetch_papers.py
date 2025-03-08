from typing import List, Dict, Optional
from Bio import Entrez

def fetch_papers(query: str, max_results: int = 10, debug: bool = False) -> List[Dict[str, str]]:
    """
    Fetch papers from PubMed based on a query.

    Args:
        query (str): PubMed search query.
        max_results (int): Maximum number of results to fetch (default: 10).
        debug (bool): Enable debug output (default: False).

    Returns:
        List[Dict[str, str]]: List of papers with extracted data.
    """

    Entrez.email = "saragadamsunil7@gmail.com"
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
    except:
        print(f"Error searching for keyword {query}: {e}")
        continue
    finally:
        handle.close()

    if debug:
        print("PubMed IDs found:", record["IdList"])

    paper_ids: List[str] = record["IdList"]
    papers: List[Dict[str, str]] = []

    for paper_id in paper_ids:
        if debug:
            print("Fetching details for paper ID:", paper_id)

        try:
            handle = Entrez.efetch(db="pubmed", id=','.join(paper_ids), retmode="xml")
            paper_data = Entrez.read(handle)
        except:
            print(f"Error fetching paper {paper_id}: {e}")
            continue
        finally:
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
            title: str = article.get("ArticleTitle", "No title available")

            # Extract publication date
            publication_date: str = ""
            article_dates = article.get("ArticleDate", [])
            if article_dates: 
                pub_date = article_dates[0]
                publication_date = f"{pub_date.get('Year', '')}-{pub_date.get('Month', '')}-{pub_date.get('Day', '')}"

            # Extract authors and affiliations
            authors: List[Dict[str, str]] = []
            for author in article.get("AuthorList", []):
                author_name: str = f"{author.get('LastName', '')}, {author.get('ForeName', '')}"
                affiliations: List[str] = [affil.get("Affiliation", "") for affil in author.get("AffiliationInfo", [])]
                authors.append({
                    "name": author_name,
                    "affiliations": affiliations
                })

                # Debug: Print author details
                if debug:
                    print("Author:", author_name)
                    print("Affiliations:", affiliations)
                    print("Identifiers:", author.get("Identifier", []))

            # Extract corresponding author email (if available)
            corresponding_email: str = ""
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
                "CorrespondingAuthorEmail": corresponding_email if corresponding_email else "Not available"
            })

        except KeyError as e:
            print(f"Error parsing paper {paper_id}: {e}")
            continue

    return papers