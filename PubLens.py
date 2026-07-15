
import streamlit as st
import pandas as pd
from Bio import Entrez, Medline
import requests
import io
from datetime import datetime
import calendar
import urllib3
import time
import re
import http.client
from io import BytesIO

try:
    from PyPDF2 import PdfReader
except Exception:
    PdfReader = None

# --- SSL FIX ---
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# --- CONFIGURATION ---
st.set_page_config(page_title="PubLens", page_icon="🔬", layout="wide")

# NOTE: This key was explicitly provided by the user in this chat.
NCBI_API_KEY = "58bebfaec2c401d129ce9f3f48996b66c109"

REQUEST_HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
        "(KHTML, like Gecko) Chrome/124.0 Safari/537.36"
    )
}

# --- INTERNAL DATABASE (2023 JIF) ---
INTERNAL_METRICS = {
    "nature": 50.5, "science": 44.7, "cell": 45.5,
    "pnas": 9.6, "proc natl acad sci u s a": 9.6,
    "nature communications": 14.7, "nat commun": 14.7,
    "scientific reports": 3.8, "sci rep": 3.8,
    "plos one": 2.9, "elife": 6.4,
    "new england journal of medicine": 96.2, "n engl j med": 96.2, "nejm": 96.2,
    "the lancet": 98.4, "lancet": 98.4,
    "jama": 63.1, "j am med assoc": 63.1,
    "bmj": 93.6, "british medical journal": 93.6,
    "nature medicine": 58.7, "nat med": 58.7,
    "cancer cell": 48.8,
    "cancer discovery": 28.2, "cancer discov": 28.2,
    "journal of clinical oncology": 42.1, "j clin oncol": 42.1,
    "clinical cancer research": 10.0, "clin cancer res": 10.0,
    "annals of oncology": 34.4, "ann oncol": 34.4,
    "cancer research": 12.5, "cancer res": 12.5,
    "molecular cancer": 27.7, "mol cancer": 27.7,
    "immunity": 25.5,
    "nature immunology": 27.7, "nat immunol": 27.7,
    "journal of experimental medicine": 12.6, "j exp med": 12.6,
    "science immunology": 17.6, "sci immunol": 17.6,
    "nature methods": 36.1, "nat methods": 36.1,
    "nature biotechnology": 33.1, "nat biotechnol": 33.1,
    "nature genetics": 31.7, "nat genet": 31.7,
    "cell stem cell": 19.8,
    "molecular cell": 14.5, "mol cell": 14.5,
    "neuron": 14.7,
    "journal of biological chemistry": 4.8, "j biol chem": 4.8, "jbc": 4.8,
    "current biology": 8.1, "curr biol": 8.1,
    "nature materials": 37.2, "nat mater": 37.2,
    "advanced materials": 27.4, "adv mater": 27.4,
    "journal of the american chemical society": 14.4, "j am chem soc": 14.4, "jacs": 14.4,
    "angewandte chemie": 16.1, "angew chem int ed engl": 16.1,
    "acs nano": 15.8,
    "nano letters": 9.6, "nano lett": 9.6,
    "advanced functional materials": 18.5, "adv funct mater": 18.5,
    "small": 13.0,
    "nature nanotechnology": 38.1, "nat nanotechnol": 38.1,
    "nature physics": 17.6, "nat phys": 17.6,
    "physical review letters": 8.1, "phys rev lett": 8.1,
    "nature reviews cancer": 72.5, "nat rev cancer": 72.5,
    "nature reviews immunology": 88.1, "nat rev immunol": 88.1,
    "nature reviews genetics": 39.1, "nat rev genet": 39.1,
    "nature reviews drug discovery": 110.5, "nat rev drug discov": 110.5,
    "cell research": 28.1, "cell res": 28.1
}

# =============================================================================
# HELPERS
# =============================================================================

def normalize_journal_name(name):
    if not name:
        return ""
    clean = str(name).lower().replace(".", "")
    return " ".join(clean.split())


def get_impact_factor(journal_name):
    if not journal_name:
        return 0.0
    clean_name = normalize_journal_name(journal_name)
    if clean_name in INTERNAL_METRICS:
        return INTERNAL_METRICS[clean_name]
    for key in INTERNAL_METRICS:
        if key in clean_name and len(key) > 4:
            return INTERNAL_METRICS[key]
    return 0.0


def determine_article_type(type_string_or_list):
    if not type_string_or_list:
        return "Primary Research"
    if isinstance(type_string_or_list, str):
        type_list = [type_string_or_list]
    else:
        type_list = type_string_or_list
    for t in type_list:
        t_str = str(t).lower()
        if "preprint" in t_str:
            return "Preprint"
        if "review" in t_str:
            return "Review"
        if any(x in t_str for x in ["conference", "meeting", "abstract", "proceeding", "poster"]):
            return "Conference/Poster"
    return "Primary Research"


def parse_pub_date(date_str):
    if not date_str or pd.isna(date_str):
        return pd.NaT
    date_str = str(date_str).strip()
    formats = ["%Y-%m-%d", "%Y %b %d", "%Y %B %d", "%Y %b", "%Y %B", "%Y"]
    for fmt in formats:
        try:
            return pd.to_datetime(datetime.strptime(date_str, fmt))
        except ValueError:
            continue
    try:
        return pd.to_datetime(date_str, errors="coerce")
    except Exception:
        return pd.NaT


def generate_citation(row):
    try:
        authors = str(row.get("Authors", "")).split(",")[0]
        if "," in str(row.get("Authors", "")):
            authors += ", et al"
        journal = row.get("Abbr", "")
        if not journal or journal == "nan":
            journal = row.get("Journal", "")
        date = str(row.get("Date", ""))
        doi = str(row.get("DOI", ""))
        details = date
        vol = str(row.get("Vol", ""))
        if vol and vol != "nan":
            details += f";{vol}"
            issue = str(row.get("Issue", ""))
            if issue and issue != "nan":
                details += f"({issue})"
        pages = str(row.get("Pages", ""))
        if pages and pages != "nan":
            details += f":{pages}"
        citation = f"{authors}. {journal}. {details}."
        if doi and doi != "nan":
            citation += f" DOI: {doi}."
        return citation
    except Exception:
        return "Error generating citation."


def normalize_phrase_text(text):
    if not text:
        return ""
    text = str(text).lower()
    text = re.sub(r"\s+", " ", text)
    return text.strip()


def phrase_pattern(phrase):
    words = [re.escape(w) for w in normalize_phrase_text(phrase).split() if w]
    if not words:
        return None
    return re.compile(r"(?<!\w)" + r"\s+".join(words) + r"(?!\w)", re.IGNORECASE)


def contains_exact_phrase(text, phrase):
    pat = phrase_pattern(phrase)
    if pat is None:
        return False
    return bool(pat.search(normalize_phrase_text(text)))


def extract_visible_text_from_html(html):
    if not html:
        return ""
    html = re.sub(r"(?is)<script.*?>.*?</script>", " ", html)
    html = re.sub(r"(?is)<style.*?>.*?</style>", " ", html)
    html = re.sub(r"(?s)<[^>]+>", " ", html)
    html = re.sub(r"&nbsp;|&#160;", " ", html)
    html = re.sub(r"&amp;", "&", html)
    html = re.sub(r"\s+", " ", html)
    return html.strip()


def fetch_url_text(url, timeout=18, max_pdf_pages=25):
    try:
        r = requests.get(url, headers=REQUEST_HEADERS, timeout=timeout, verify=False, allow_redirects=True)
        r.raise_for_status()
        content_type = r.headers.get("Content-Type", "").lower()
        final_url = r.url.lower()
        if ("pdf" in content_type or final_url.endswith(".pdf")) and PdfReader is not None:
            reader = PdfReader(BytesIO(r.content))
            text_chunks = []
            for page in reader.pages[:max_pdf_pages]:
                try:
                    text_chunks.append(page.extract_text() or "")
                except Exception:
                    continue
            return normalize_phrase_text(" ".join(text_chunks))
        return normalize_phrase_text(extract_visible_text_from_html(r.text))
    except Exception:
        return ""


# =============================================================================
# SOURCE 1: PMC FULL TEXT
# =============================================================================

def search_pmc_fulltext(query, max_results, email, start_date, end_date, search_mode):
    Entrez.email = email
    Entrez.api_key = NCBI_API_KEY
    mindate = start_date.strftime("%Y/%m/%d")
    maxdate = end_date.strftime("%Y/%m/%d")
    # Build query based on search_mode
    if search_mode == "Exact Phrase":
        final_query = f'"{query.strip()}"'
    elif search_mode == "All Terms (AND)":
        terms = [t.strip() for t in query.split(",") if t.strip()]
        final_query = " AND ".join([f'"{t}"' if " " in t else t for t in terms])
    elif search_mode == "Any Terms (OR)":
        terms = [t.strip() for t in query.split(",") if t.strip()]
        final_query = " OR ".join([f'"{t}"' if " " in t else t for t in terms])
    else:
        final_query = f'"{query.strip()}"'

    try:
        search_handle = Entrez.esearch(
            db="pmc",
            term=final_query,
            retmax=max_results,
            mindate=mindate,
            maxdate=maxdate,
            datetype="pdat"
        )
        pmc_results = Entrez.read(search_handle, validate=False)
        search_handle.close()
        pmc_ids = pmc_results.get("IdList", [])
        if not pmc_ids:
            return pd.DataFrame(), {"pmc_candidates": 0, "pmc_link_loss": 0}

        linked_pubmed_ids = []
        batch_size = 100
        for i in range(0, len(pmc_ids), batch_size):
            chunk = pmc_ids[i:i + batch_size]
            try:
                link_handle = Entrez.elink(dbfrom="pmc", db="pubmed", id=chunk)
                link_results = Entrez.read(link_handle, validate=False)
                link_handle.close()
                for linkset in link_results:
                    if "LinkSetDb" in linkset and len(linkset["LinkSetDb"]) > 0:
                        for link in linkset["LinkSetDb"][0]["Link"]:
                            linked_pubmed_ids.append(link["Id"])
            except Exception:
                pass
            time.sleep(0.1)

        link_loss = max(0, len(pmc_ids) - len(linked_pubmed_ids))
        if not linked_pubmed_ids:
            return pd.DataFrame(), {"pmc_candidates": len(pmc_ids), "pmc_link_loss": link_loss}

        all_records = []
        fetch_batch = 50
        for i in range(0, len(linked_pubmed_ids), fetch_batch):
            chunk = linked_pubmed_ids[i:i + fetch_batch]
            success = False
            for attempt in range(5):
                try:
                    fetch_handle = Entrez.efetch(db="pubmed", id=chunk, rettype="medline", retmode="text")
                    chunk_records = list(Medline.parse(fetch_handle))
                    all_records.extend(chunk_records)
                    success = True
                    break
                except (http.client.IncompleteRead, Exception):
                    time.sleep(1 * (attempt + 1))
            if not success:
                st.warning(f"PMC batch {i}-{i + fetch_batch} skipped after retries.")
            time.sleep(0.1)

        data = []
        for record in all_records:
            j_name = record.get("JT", "")
            j_abbr = record.get("TA", "")
            if_score = get_impact_factor(j_name)
            if if_score == 0.0:
                if_score = get_impact_factor(j_abbr)
            doi = ""
            for aid in record.get("AID", []):
                if "[doi]" in aid:
                    doi = aid.replace(" [doi]", "")
                    break
            if not doi:
                for lid in record.get("LID", []):
                    if "[doi]" in lid:
                        doi = lid.replace(" [doi]", "")
                        break
            data.append({
                "Select": False,
                "Source": "🇺🇸 PMC",
                "Type": determine_article_type(record.get("PT", [])),
                "Journal": j_name,
                "Abbr": j_abbr,
                "2023 IF": if_score,
                "Title": record.get("TI", ""),
                "Authors": ", ".join(record.get("AU", [])),
                "Date": record.get("DP", ""),
                "Vol": record.get("VI", ""),
                "Issue": record.get("IP", ""),
                "Pages": record.get("PG", ""),
                "DOI": doi,
                "PMID": record.get("PMID", ""),
                "Link": "",
                "CitedBy": 0,
            })
        return pd.DataFrame(data), {"pmc_candidates": len(pmc_ids), "pmc_link_loss": link_loss}
    except Exception as e:
        st.error(f"PMC Search Error: {e}")
        return pd.DataFrame(), {"pmc_candidates": 0, "pmc_link_loss": 0}


# =============================================================================
# SOURCE 2: EUROPE PMC FULL TEXT
# =============================================================================

def search_europe_pmc(query, max_results, start_date, end_date, email, include_preprints, search_mode):
    base_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    headers = {"User-Agent": f"PubLens ({email})"}
    s_str = start_date.strftime("%Y-%m-%d")
    e_str = end_date.strftime("%Y-%m-%d")
    # Build query based on search_mode
    if search_mode == "Exact Phrase":
        term = f'"{query.strip()}"'
    elif search_mode == "All Terms (AND)":
        terms = [t.strip() for t in query.split(",") if t.strip()]
        term = " AND ".join([f'"{t}"' if " " in t else t for t in terms])
    elif search_mode == "Any Terms (OR)":
        terms = [t.strip() for t in query.split(",") if t.strip()]
        term = " OR ".join([f'"{t}"' if " " in t else t for t in terms])
    else:
        term = f'"{query.strip()}"'
    date_part = f" AND FIRST_PDATE:[{s_str} TO {e_str}]"
    source_filter = " AND (SRC:MED OR SRC:PMC OR SRC:PPR)" if include_preprints else " AND (SRC:MED OR SRC:PMC)"
    full_query = f"{term}{date_part}{source_filter}"

    all_results = []
    cursor = "*"
    page_size = 25
    while len(all_results) < max_results:
        params = {
            "query": full_query,
            "format": "json",
            "pageSize": min(page_size, max_results - len(all_results)),
            "resultType": "core",
            "cursorMark": cursor,
        }
        success = False
        for attempt in range(3):
            try:
                response = requests.get(base_url, params=params, headers=headers, verify=False, timeout=15)
                response.raise_for_status()
                data = response.json()
                result_list = data.get("resultList", {}).get("result", [])
                if not result_list:
                    return pd.DataFrame(all_results)
                for item in result_list:
                    j_name = item.get("journalTitle", "") or item.get("journalInfo", {}).get("journal", {}).get("title", "") or ""
                    j_abbr = item.get("journalInfo", {}).get("journal", {}).get("medlineAbbreviation", "") or ""
                    pub_types = item.get("pubTypeList", {}).get("pubType", [])
                    all_results.append({
                        "Select": False,
                        "Source": "🇪🇺 Europe PMC",
                        "Type": determine_article_type(pub_types),
                        "Journal": j_name,
                        "Abbr": j_abbr,
                        "2023 IF": get_impact_factor(j_name) or get_impact_factor(j_abbr),
                        "Title": item.get("title", ""),
                        "Authors": item.get("authorString", ""),
                        "Date": item.get("firstPublicationDate", ""),
                        "Vol": item.get("journalInfo", {}).get("volume", ""),
                        "Issue": item.get("journalInfo", {}).get("issue", ""),
                        "Pages": item.get("pageInfo", ""),
                        "DOI": item.get("doi", "") or "",
                        "PMID": item.get("pmid", "") or "",
                        "Link": "",
                        "CitedBy": int(item.get("citedByCount", 0) or 0),
                    })
                    if len(all_results) >= max_results:
                        break
                next_cursor = data.get("nextCursorMark", "")
                if not next_cursor or next_cursor == cursor:
                    return pd.DataFrame(all_results)
                cursor = next_cursor
                success = True
                break
            except requests.exceptions.HTTPError as err:
                if getattr(err.response, "status_code", None) == 503:
                    time.sleep(2 * (attempt + 1))
                    continue
                return pd.DataFrame(all_results)
            except Exception:
                time.sleep(1)
        if not success:
            return pd.DataFrame(all_results)
        time.sleep(0.3)
    return pd.DataFrame(all_results)


# =============================================================================
# SOURCE 3: GOOGLE SCHOLAR + LANDING PAGE/PDF VALIDATION
# =============================================================================

def search_google_scholar(query, api_key, max_results, start_year, end_year, search_mode, validate_links=True):
    base_url = "https://serpapi.com/search"
    raw_count = 0
    kept_results = []
    link_validated_kept = 0
    # Build query and validation terms
    if search_mode == "Exact Phrase":
        query_norm = normalize_phrase_text(query)
        serp_query = f'"{query.strip()}"'
        serp_as_epq = query.strip()
        and_terms = []
        or_terms = []
    elif search_mode == "All Terms (AND)":
        terms = [t.strip() for t in query.split(",") if t.strip()]
        query_norm = terms
        serp_query = " ".join([f'"{t}"' if " " in t else t for t in terms])
        serp_as_epq = ""
        and_terms = terms
        or_terms = []
    elif search_mode == "Any Terms (OR)":
        terms = [t.strip() for t in query.split(",") if t.strip()]
        query_norm = terms
        serp_query = " OR ".join([f'"{t}"' if " " in t else t for t in terms])
        serp_as_epq = ""
        and_terms = []
        or_terms = terms
    else:
        query_norm = normalize_phrase_text(query)
        serp_query = f'"{query.strip()}"'
        serp_as_epq = query.strip()
        and_terms = []
        or_terms = []

    # Search more candidates than needed because validation will discard some.
    candidate_budget = max(max_results * 3, max_results + 30)
    max_pages = (candidate_budget + 9) // 10

    for page in range(max_pages):
        start_index = page * 10
        if len(kept_results) >= max_results:
            break

        params = {
            "engine": "google_scholar",
            "q": serp_query,
            "as_epq": serp_as_epq,
            "api_key": api_key,
            "start": start_index,
            "as_ylo": start_year,
            "as_yhi": end_year,
        }

        try:
            response = requests.get(base_url, params=params, headers=REQUEST_HEADERS, verify=False, timeout=20)
            response.raise_for_status()
            data = response.json()
            organic = data.get("organic_results", [])
            if not organic:
                break

            for item in organic:
                raw_count += 1
                title = item.get("title", "") or ""
                snippet = item.get("snippet", "") or ""
                pub_info = item.get("publication_info", {}).get("summary", "") or ""
                link = item.get("link", "") or ""
                cited_by = item.get("inline_links", {}).get("cited_by", {}).get("total", 0)

                local_text = normalize_phrase_text(f"{title} {snippet} {pub_info}")
                keep = False
                matched_by = ""

                # Validation logic for each mode
                if search_mode == "Exact Phrase":
                    if contains_exact_phrase(local_text, query_norm):
                        keep = True
                        matched_by = "Snippet/Title"
                    elif validate_links and link:
                        page_text = fetch_url_text(link)
                        if contains_exact_phrase(page_text, query_norm):
                            keep = True
                            matched_by = "Landing Page/PDF"
                            link_validated_kept += 1
                elif search_mode == "All Terms (AND)":
                    if all(term.lower() in local_text for term in query_norm):
                        keep = True
                        matched_by = "Snippet/Title"
                    elif validate_links and link:
                        page_text = fetch_url_text(link)
                        if all(term.lower() in page_text for term in query_norm):
                            keep = True
                            matched_by = "Landing Page/PDF"
                            link_validated_kept += 1
                elif search_mode == "Any Terms (OR)":
                    if any(term.lower() in local_text for term in query_norm):
                        keep = True
                        matched_by = "Snippet/Title"
                    elif validate_links and link:
                        page_text = fetch_url_text(link)
                        if any(term.lower() in page_text for term in query_norm):
                            keep = True
                            matched_by = "Landing Page/PDF"
                            link_validated_kept += 1
                else:
                    if contains_exact_phrase(local_text, query_norm):
                        keep = True
                        matched_by = "Snippet/Title"
                    elif validate_links and link:
                        page_text = fetch_url_text(link)
                        if contains_exact_phrase(page_text, query_norm):
                            keep = True
                            matched_by = "Landing Page/PDF"
                            link_validated_kept += 1

                if not keep:
                    continue

                doi = ""
                if "doi.org/" in link:
                    doi = link.split("doi.org/")[-1]

                authors_raw = ""
                journal_raw = ""
                year_str = ""
                parts = pub_info.split(" - ")
                if len(parts) >= 2:
                    authors_raw = parts[0].strip()
                    journal_year = parts[1].strip()
                    jy_parts = journal_year.rsplit(",", 1)
                    if len(jy_parts) == 2:
                        journal_raw = jy_parts[0].strip()
                        year_match = re.search(r"(\d{4})", jy_parts[1])
                        year_str = year_match.group(1) if year_match else ""
                    else:
                        journal_raw = journal_year
                elif len(parts) == 1:
                    authors_raw = parts[0].strip()

                check_str = (pub_info + " " + title).lower()
                article_type = "Unknown (Scholar)"
                if "[citation]" in check_str:
                    article_type = "Citation Only"
                elif any(x in check_str for x in ["conference", "meeting", "abstracts of", "proceedings", "poster"]):
                    article_type = "Conference/Poster"
                elif any(x in check_str for x in ["biorxiv", "medrxiv", "preprint"]):
                    article_type = "Preprint"
                elif "review" in check_str:
                    article_type = "Review"
                else:
                    article_type = "Primary Research"

                kept_results.append({
                    "Select": False,
                    "Source": "🎓 Scholar",
                    "Type": article_type,
                    "Journal": journal_raw,
                    "Abbr": "",
                    "2023 IF": get_impact_factor(journal_raw),
                    "Title": title,
                    "Authors": authors_raw,
                    "Date": year_str,
                    "Vol": "",
                    "Issue": "",
                    "Pages": "",
                    "DOI": doi,
                    "PMID": "",
                    "Link": link,
                    "CitedBy": cited_by,
                    "PhraseMatch": matched_by,
                })
                if len(kept_results) >= max_results:
                    break
        except Exception as e:
            st.warning(f"Scholar error: {e}")
            break
        time.sleep(0.5)

    stats = {
        "scholar_raw": raw_count,
        "scholar_kept": len(kept_results),
        "scholar_filtered_out": max(0, raw_count - len(kept_results)),
        "scholar_link_validated_kept": link_validated_kept,
    }
    return pd.DataFrame(kept_results), stats


# =============================================================================
# ENRICHMENT
# =============================================================================

def enrich_with_pubmed(title, doi, email):
    Entrez.email = email
    Entrez.api_key = NCBI_API_KEY
    search_terms = []
    if doi and len(str(doi)) > 5:
        search_terms.append(f"{doi}[DOI]")
    if title and len(str(title)) > 10:
        clean_title = re.sub(r"[^a-zA-Z0-9 ]", "", str(title))
        search_terms.append(f"{clean_title}[Title]")

    for search_term in search_terms:
        try:
            handle = Entrez.esearch(db="pubmed", term=search_term, retmax=1)
            result = Entrez.read(handle, validate=False)
            handle.close()
            ids = result.get("IdList", [])
            if ids:
                fetch_handle = Entrez.efetch(db="pubmed", id=[ids[0]], rettype="medline", retmode="text")
                records = list(Medline.parse(fetch_handle))
                if records:
                    rec = records[0]
                    rec_doi = ""
                    for aid in rec.get("AID", []):
                        if "[doi]" in aid:
                            rec_doi = aid.replace(" [doi]", "")
                            break
                    if not rec_doi:
                        for lid in rec.get("LID", []):
                            if "[doi]" in lid:
                                rec_doi = lid.replace(" [doi]", "")
                                break
                    return {
                        "Journal": rec.get("JT", ""),
                        "Abbr": rec.get("TA", ""),
                        "Type": determine_article_type(rec.get("PT", [])),
                        "Authors": ", ".join(rec.get("AU", [])),
                        "Date": rec.get("DP", ""),
                        "Vol": rec.get("VI", ""),
                        "Issue": rec.get("IP", ""),
                        "Pages": rec.get("PG", ""),
                        "DOI": rec_doi,
                        "PMID": rec.get("PMID", ""),
                    }
        except Exception:
            pass
        time.sleep(0.1)
    return None


def enrich_missing(df, email, progress_bar):
    total = len(df)
    enriched_count = 0
    for i, idx in enumerate(df.index):
        row = df.loc[idx]
        has_pmid = str(row.get("PMID", "")).strip() not in ["", "nan"]
        has_doi = str(row.get("DOI", "")).strip() not in ["", "nan"] and len(str(row.get("DOI", ""))) > 5
        has_full_date = len(str(row.get("Date", "")).strip()) > 4
        if has_pmid and has_doi and has_full_date:
            progress_bar.progress((i + 1) / total, text=f"Enriching {i+1}/{total} (already complete)")
            continue
        metadata = enrich_with_pubmed(str(row.get("Title", "")), str(row.get("DOI", "")), email)
        if metadata:
            enriched_count += 1
            for field in ["Journal", "Abbr", "Type", "Authors", "Date", "Vol", "Issue", "Pages", "DOI", "PMID"]:
                new_val = str(metadata.get(field, "")).strip()
                old_val = str(df.at[idx, field]).strip() if field in df.columns else ""
                if new_val and new_val != "nan" and (not old_val or old_val == "nan" or len(new_val) > len(old_val)):
                    df.at[idx, field] = new_val
            j_name = str(df.at[idx, "Journal"])
            j_abbr = str(df.at[idx, "Abbr"])
            if_score = get_impact_factor(j_name)
            if if_score == 0.0:
                if_score = get_impact_factor(j_abbr)
            df.at[idx, "2023 IF"] = if_score
        progress_bar.progress((i + 1) / total, text=f"Enriching {i+1}/{total}...")
    return df, enriched_count


# =============================================================================
# DEDUPLICATION
# =============================================================================

def dedup_and_merge(frames):
    df = pd.concat(frames, ignore_index=True)
    if df.empty:
        return df
    df["DOI_clean"] = df["DOI"].fillna("").astype(str).str.strip().str.lower()
    df["PMID_clean"] = df["PMID"].fillna("").astype(str).str.strip()
    df["_clean_title"] = df["Title"].fillna("").astype(str).str.lower().str.replace(r"[^a-z0-9]", "", regex=True)
    df["_dedup_key"] = df["DOI_clean"]
    mask = df["_dedup_key"] == ""
    df.loc[mask, "_dedup_key"] = df.loc[mask, "PMID_clean"]
    mask = df["_dedup_key"] == ""
    df.loc[mask, "_dedup_key"] = df.loc[mask, "_clean_title"]

    df["_score"] = 0
    df.loc[df["DOI_clean"].str.len() > 5, "_score"] += 20
    df.loc[df["PMID_clean"].str.len() > 0, "_score"] += 15
    df.loc[df["Abbr"].fillna("").astype(str).str.len() > 1, "_score"] += 10
    df.loc[df["Date"].fillna("").astype(str).str.len() > 4, "_score"] += 10
    df.loc[df["Source"].astype(str).str.contains("PMC"), "_score"] += 5
    df.loc[df["Source"].astype(str).str.contains("Europe"), "_score"] += 4

    source_map = df.groupby("_dedup_key")["Source"].apply(lambda x: " | ".join(sorted(set(x)))).reset_index().rename(columns={"Source": "_merged_source"})
    cited_map = df.groupby("_dedup_key")["CitedBy"].max().reset_index().rename(columns={"CitedBy": "_max_cited"})

    df = df.sort_values("_score", ascending=False)
    df = df.drop_duplicates(subset=["_dedup_key"], keep="first")
    df = df.drop(columns=["Source", "CitedBy"], errors="ignore")
    df = df.merge(source_map, on="_dedup_key", how="left")
    df = df.merge(cited_map, on="_dedup_key", how="left")
    df = df.rename(columns={"_merged_source": "Source", "_max_cited": "CitedBy"})
    df = df.drop(columns=["_score", "DOI_clean", "PMID_clean"], errors="ignore")
    return df


def to_excel(df):
    df_export = df.drop(columns=["Select", "_clean_title", "_dedup_key", "PhraseMatch"], errors="ignore")
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        df_export.to_excel(writer, index=False, sheet_name='Results')
    return output.getvalue()


# =============================================================================
# UI
# =============================================================================

st.title("🔬 PubLens")
st.markdown(
    "Find research papers and related records that mention a Leica technology term. "
    "The app can search **PMC full text**, **Europe PMC**, and **Google Scholar** using exact-phrase, all-terms, or any-terms matching, then validates matches to reduce false positives."
)

with st.sidebar:
    st.header("Search")
    search_mode = st.selectbox(
        "Search mode",
        ["Exact Phrase", "All Terms (AND)", "Any Terms (OR)"],
        help="Choose how to match your search terms."
    )
    if search_mode == "Exact Phrase":
        query = st.text_input("Technology name", value="Cell DIVE", help="Exact phrase only")
    else:
        query = st.text_input(
            "Search terms (comma-separated)",
            value="Cell DIVE, Leica",
            help="Enter terms separated by commas."
        )
    email = st.text_input("Email (required for NCBI)", placeholder="researcher@lab.edu")
    serp_key = st.text_input("SerpApi Key (optional)", type="password", help="Needed only for Google Scholar")
    st.caption("🔑 NCBI API key configured (10 req/sec)")

    st.subheader("Sources")
    use_pmc = st.checkbox("PMC full text", value=True)
    use_eu = st.checkbox("Europe PMC", value=True)
    use_gs = st.checkbox("Google Scholar (validated)", value=bool(serp_key), disabled=not bool(serp_key))

    st.subheader("Options")
    include_preprints = st.checkbox("Include BioRxiv/MedRxiv", value=False)
    do_enrich = st.checkbox("Enrich incomplete records from PubMed", value=True)
    limit = st.slider("Max results per source", 10, 500, 200)

    st.subheader("Date range")
    today = datetime.today()
    years = list(range(today.year, 1980, -1))
    months = list(calendar.month_name)[1:]
    c1, c2 = st.columns(2)
    with c1:
        s_month = st.selectbox("Start month", months, index=0)
    with c2:
        s_year = st.selectbox("Start year", years, index=min(3, len(years) - 1))
    c3, c4 = st.columns(2)
    with c3:
        e_month = st.selectbox("End month", months, index=today.month - 1)
    with c4:
        e_year = st.selectbox("End year", years, index=0)

    search_btn = st.button("Search", type="primary", use_container_width=True)

# =============================================================================
# MAIN
# =============================================================================

if search_btn:
    if not email:
        st.warning("Please enter your email.")
        st.stop()
    if not query.strip():
        st.warning("Please enter a technology name.")
        st.stop()
    if not (use_pmc or use_eu or use_gs):
        st.warning("Please select at least one source.")
        st.stop()

    month_map = {m: i + 1 for i, m in enumerate(months)}
    start_dt = datetime(s_year, month_map[s_month], 1)
    last_day = calendar.monthrange(e_year, month_map[e_month])[1]
    end_dt = datetime(e_year, month_map[e_month], last_day)

    frames = []
    stats = {
        "pmc_candidates": 0,
        "pmc_link_loss": 0,
        "scholar_raw": 0,
        "scholar_kept": 0,
        "scholar_filtered_out": 0,
        "scholar_link_validated_kept": 0,
        "enriched": 0,
    }

    if use_pmc:
        with st.status("Searching PMC full text...", expanded=False):
            df_pmc, s = search_pmc_fulltext(query, limit, email, start_dt, end_dt, search_mode)
            stats.update({**stats, **s})
            if not df_pmc.empty:
                frames.append(df_pmc)

    if use_eu:
        with st.status("Searching Europe PMC...", expanded=False):
            df_eu = search_europe_pmc(query, limit, start_dt, end_dt, email, include_preprints, search_mode)
            if not df_eu.empty:
                frames.append(df_eu)

    if use_gs and serp_key:
        with st.status("Searching Google Scholar and validating links...", expanded=False):
            df_gs, s = search_google_scholar(query, serp_key, limit, s_year, e_year, search_mode, validate_links=True)
            stats.update({**stats, **s})
            if not df_gs.empty:
                frames.append(df_gs)

    if not frames:
        st.warning("No results found.")
        st.stop()

    raw_total = sum(len(f) for f in frames)
    df = dedup_and_merge(frames)
    duplicates_removed = raw_total - len(df)

    if do_enrich and not df.empty:
        with st.status("Enriching incomplete records...", expanded=False):
            progress_bar = st.progress(0, text="Starting enrichment...")
            df, enriched_count = enrich_missing(df, email, progress_bar)
            progress_bar.empty()
            stats["enriched"] = enriched_count

    df["ParsedDate"] = df["Date"].apply(parse_pub_date)
    df["2023 IF"] = pd.to_numeric(df["2023 IF"], errors='coerce').fillna(0.0)
    df["CitedBy"] = pd.to_numeric(df["CitedBy"], errors='coerce').fillna(0).astype(int)

    # Prioritize primary research, then IF, then citations.
    type_rank = {
        "Primary Research": 0,
        "Unknown (Scholar)": 1,
        "Preprint": 2,
        "Review": 3,
        "Conference/Poster": 4,
        "Citation Only": 5,
    }
    df["_type_rank"] = df["Type"].map(type_rank).fillna(9)
    df = df.sort_values(by=["_type_rank", "2023 IF", "CitedBy"], ascending=[True, False, False]).reset_index(drop=True)

    st.session_state["results_df"] = df
    st.session_state["stats"] = stats
    st.session_state["duplicates_removed"] = duplicates_removed
    st.session_state["raw_total"] = raw_total

if 'results_df' in st.session_state:
    df = st.session_state['results_df']
    stats = st.session_state.get('stats', {})
    raw_total = st.session_state.get('raw_total', len(df))
    duplicates_removed = st.session_state.get('duplicates_removed', 0)

    primary_count = int((df["Type"] == "Primary Research").sum()) if "Type" in df.columns else 0
    dated = df[df["ParsedDate"].notna()] if "ParsedDate" in df.columns else pd.DataFrame()

    m1, m2, m3, m4, m5, m6 = st.columns(6)
    m1.metric("Unique papers", len(df))
    m2.metric("Primary research", primary_count)
    m3.metric("Duplicates removed", duplicates_removed)
    m4.metric("Scholar kept", stats.get("scholar_kept", 0))
    m5.metric("Scholar discarded", stats.get("scholar_filtered_out", 0))
    if not dated.empty:
        newest = dated.sort_values("ParsedDate", ascending=False).iloc[0]
        m6.metric("Most recent", newest["ParsedDate"].strftime("%b %Y"))
        st.info(f"Most recent: **{newest['Title'][:140]}** — *{newest.get('Journal', '')}* ({newest['ParsedDate'].strftime('%B %Y')})")

    expl = []
    if stats.get("pmc_link_loss", 0) > 0:
        expl.append(f"PMC link loss: {stats['pmc_link_loss']} PMC records had no PubMed link")
    if stats.get("scholar_link_validated_kept", 0) > 0:
        expl.append(f"Scholar link/PDF validations kept: {stats['scholar_link_validated_kept']}")
    if stats.get("enriched", 0) > 0:
        expl.append(f"Enriched from PubMed: {stats['enriched']}")
    if expl:
        st.caption(" • ".join(expl))

    display_cols = [
        "Select", "Source", "Type", "Journal", "2023 IF", "CitedBy",
        "Title", "Date", "DOI", "PMID", "Abbr", "Vol", "Issue", "Pages", "Link"
    ]
    final_cols = [c for c in display_cols if c in df.columns]

    edited_df = st.data_editor(
        df[final_cols],
        use_container_width=True,
        column_config={
            "Select": st.column_config.CheckboxColumn("Select"),
            "Type": st.column_config.TextColumn("Type", width="small"),
            "2023 IF": st.column_config.NumberColumn("2023 IF", format="%.1f"),
            "CitedBy": st.column_config.NumberColumn("Cited", format="%d"),
            "DOI": st.column_config.LinkColumn("DOI"),
            "Link": st.column_config.LinkColumn("Link"),
        },
        hide_index=True,
    )

    c_a, c_b = st.columns([1, 4])
    with c_a:
        st.download_button("📥 Download Excel", data=to_excel(df), file_name="PubLens_Results.xlsx")
    with c_b:
        if st.button("📝 Generate Citations"):
            selected = edited_df[edited_df["Select"] == True]
            if not selected.empty:
                txt = ""
                for idx, row in selected.iterrows():
                    try:
                        txt += generate_citation(df.loc[idx]) + "\n\n"
                    except Exception:
                        pass
                st.text_area("Citations", value=txt.strip(), height=220)
            else:
                st.warning("Select papers first.")
