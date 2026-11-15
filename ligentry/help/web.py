"""Where the Help window can send you on the web, and how it reads what you type (no Qt).

ligentry never contacts these sites by itself. They are only used when you search
from the Help window's search box or click a link, and the page then opens in the
Help window's Web tab (or in your normal browser if the Web tab is not available).

The search addresses below are written from the sites' public search pages. Sites
change them from time to time; if one stops working the others are unaffected, and
the plain web search always works.

`resolve_input` turns what you type into an address:

    https://example.org/page      opened as it is
    pubmed.ncbi.nlm.nih.gov/123   "https://" is added
    10.1000/xyz123  or  doi:...   the DOI resolver
    PMID 12345678                 the PubMed record
    PMC1234567                    the PubMed Central article
    NCT01234567                   the ClinicalTrials.gov study
    PDB 1HSG                      the RCSB PDB entry
    anything else                 a search with the source you picked

`is_safe_url` is the gate every address passes before the Web tab loads it: only
http and https, a real host name, no embedded passwords, and nothing on your own
machine or network.
"""
from __future__ import annotations

import ipaddress
import re
from dataclasses import dataclass
from urllib.parse import quote, quote_plus, urlsplit


@dataclass(frozen=True)
class Source:
    key: str
    name: str
    kind: str  # papers, clinical, chemistry, protein, video, reference, general
    description: str
    template: str  # contains {q}, replaced by the URL-encoded search text


SOURCES: tuple[Source, ...] = (
    Source("pubmed", "PubMed", "papers",
           "Citations and abstracts of biomedical papers (US National Library of Medicine).",
           "https://pubmed.ncbi.nlm.nih.gov/?term={q}"),
    Source("europepmc", "Europe PMC", "papers",
           "Papers, many with free full text, plus preprints and patents.",
           "https://europepmc.org/search?query={q}"),
    Source("pmc", "PubMed Central", "papers",
           "Free full-text archive of biomedical papers.",
           "https://pmc.ncbi.nlm.nih.gov/search/?term={q}"),
    Source("scholar", "Google Scholar", "papers",
           "Broad search over scholarly literature, including books and theses.",
           "https://scholar.google.com/scholar?q={q}"),
    Source("clinicaltrials", "ClinicalTrials.gov", "clinical",
           "Registry of clinical studies and their results.",
           "https://clinicaltrials.gov/search?term={q}"),
    Source("dailymed", "DailyMed", "clinical",
           "Current US drug labels (prescribing information).",
           "https://dailymed.nlm.nih.gov/dailymed/search.cfm?labeltype=all&query={q}"),
    Source("medlineplus", "MedlinePlus", "clinical",
           "Plain-language health information from the US National Library of Medicine.",
           "https://duckduckgo.com/?q={q}+site%3Amedlineplus.gov"),
    Source("pubchem", "PubChem", "chemistry",
           "Chemical structures, names and properties.",
           "https://pubchem.ncbi.nlm.nih.gov/#query={q}"),
    Source("chembl", "ChEMBL", "chemistry",
           "Bioactivity data for drug-like molecules and their targets.",
           "https://www.ebi.ac.uk/chembl/g/#search_results/all/query={q}"),
    Source("uniprot", "UniProt", "protein",
           "Protein sequences and annotation.",
           "https://www.uniprot.org/uniprotkb?query={q}"),
    Source("rcsb", "RCSB PDB", "protein",
           "Experimental 3D structures of proteins (searched through a web search).",
           "https://duckduckgo.com/?q={q}+site%3Arcsb.org"),
    Source("wikipedia", "Wikipedia", "reference",
           "A quick overview of a term. Check its references before relying on it.",
           "https://en.wikipedia.org/w/index.php?search={q}"),
    Source("youtube", "YouTube", "video",
           "Explanatory videos and lectures. Quality varies a great deal.",
           "https://www.youtube.com/results?search_query={q}"),
    Source("web", "Web search", "general",
           "A general web search (DuckDuckGo).",
           "https://duckduckgo.com/?q={q}"),
)

_BY_KEY = {s.key: s for s in SOURCES}
DEFAULT_SOURCE = "web"

# The sources offered next to a manual page's "look it up" links, in this order.
QUICK_SOURCES = ("pubmed", "europepmc", "scholar", "wikipedia", "youtube")


def get_source(key: str) -> Source:
    try:
        return _BY_KEY[key]
    except KeyError:
        raise ValueError(f"unknown web source '{key}'") from None


def build_search_url(source_key: str, query: str) -> str:
    """The address of a search for `query` at the given source."""
    query = " ".join(query.split())
    if not query:
        raise ValueError("Type something to search for.")
    if len(query) > 500:
        raise ValueError("The search text is too long (500 characters at most).")
    source = get_source(source_key)
    return source.template.format(q=quote_plus(query))


def is_safe_url(url: str) -> bool:
    """True for an ordinary public http(s) address; False for anything else."""
    try:
        parts = urlsplit(url.strip())
        host = parts.hostname
    except ValueError:
        return False
    if parts.scheme not in ("http", "https") or not host:
        return False
    if parts.username is not None or parts.password is not None:
        return False
    host = host.lower().rstrip(".")
    if host == "localhost" or host.endswith((".localhost", ".local", ".internal")):
        return False
    try:
        ip = ipaddress.ip_address(host)
    except ValueError:
        return "." in host  # a name must have a dot (no bare machine names)
    return not (ip.is_private or ip.is_loopback or ip.is_link_local or ip.is_reserved
                or ip.is_multicast or ip.is_unspecified)


_DOI = re.compile(r"^(?:doi:\s*|https?://(?:dx\.)?doi\.org/)?(10\.\d{4,9}/\S+)$", re.IGNORECASE)
_PMID = re.compile(r"^pmid[:\s#]*([1-9]\d{0,8})$", re.IGNORECASE)
_PMC = re.compile(r"^(PMC\d{4,9})$", re.IGNORECASE)
_NCT = re.compile(r"^(NCT\d{8})$", re.IGNORECASE)
_PDB = re.compile(r"^pdb[:\s#]*([0-9][A-Za-z0-9]{3})$", re.IGNORECASE)
_DOMAIN = re.compile(r"^(?:[a-z0-9-]+\.)+[a-z]{2,}(?:[/?#]\S*)?$", re.IGNORECASE)


@dataclass(frozen=True)
class Target:
    url: str
    label: str  # what to show in the status line, e.g. "PubMed record 123"
    is_search: bool


def resolve_input(text: str, source_key: str = DEFAULT_SOURCE) -> Target:
    """Work out where `text` should go. Raises ValueError for empty or unsafe input."""
    text = " ".join(text.split())
    if not text:
        raise ValueError("Type something to search for.")

    m = _DOI.match(text)
    if m:
        doi = m.group(1).rstrip(".,;)")
        return Target(f"https://doi.org/{quote(doi, safe='/:()')}", f"DOI {doi}", False)
    m = _PMID.match(text)
    if m:
        return Target(f"https://pubmed.ncbi.nlm.nih.gov/{m.group(1)}/", f"PubMed record {m.group(1)}", False)
    m = _PMC.match(text)
    if m:
        pmc = m.group(1).upper()
        return Target(f"https://pmc.ncbi.nlm.nih.gov/articles/{pmc}/", f"PubMed Central {pmc}", False)
    m = _NCT.match(text)
    if m:
        nct = m.group(1).upper()
        return Target(f"https://clinicaltrials.gov/study/{nct}", f"ClinicalTrials.gov {nct}", False)
    m = _PDB.match(text)
    if m:
        pdb = m.group(1).upper()
        return Target(f"https://www.rcsb.org/structure/{pdb}", f"PDB entry {pdb}", False)

    if " " not in text and (text.lower().startswith(("http://", "https://")) or _DOMAIN.match(text)):
        url = text if text.lower().startswith(("http://", "https://")) else "https://" + text
        if not is_safe_url(url):
            raise ValueError("That address cannot be opened here (only public http and https pages).")
        return Target(url, url, False)

    url = build_search_url(source_key, text)
    return Target(url, f"{get_source(source_key).name}: {text}", True)
