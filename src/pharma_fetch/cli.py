import argparse
from typing import List, Dict, Any
from pharma_fetch.fetch_papers import fetch_papers
from pharma_fetch.filter_authors import filter_authors
from pharma_fetch.generate_csv import generate_csv

def main() -> None:
    """
    Command-line interface for fetching and filtering research papers from PubMed.
    """
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Fetch and filter research papers from PubMed based on a query."
    )
    
    # Required argument: PubMed query
    parser.add_argument(
        "query",
        type=str,
        help="PubMed search query (e.g., 'cancer treatment')"
    )
    
    # Optional argument: Output file name
    parser.add_argument(
        "-f", "--file",
        type=str,
        default="output.csv",
        help="Output CSV file name (default: output.csv)"
    )
    
    # Optional argument: Enable debug output
    parser.add_argument(
        "-d", "--debug",
        action="store_true",
        help="Enable debug output"
    )
    
    # Parse arguments
    args = parser.parse_args()

    # Fetch papers
    papers: List[Dict[str, Any]] = fetch_papers(args.query, debug=args.debug)
    if args.debug:
        print("Fetched papers:", papers)

    # Filter authors
    filtered_papers: List[Dict[str, Any]] = []
    for paper in papers:
        if args.debug:
            print("Processing paper:", paper.get("PubmedID", ""))
        non_academic_authors: List[Dict[str, str]] = filter_authors(paper, debug=args.debug)
        paper["NonAcademicAuthors"] = non_academic_authors
        filtered_papers.append(paper)

    # Generate CSV
    generate_csv(filtered_papers, args.file, debug=args.debug)
    print(f"Results saved to {args.file}")

if __name__ == "__main__":
    main()