import streamlit as st
import pandas as pd
from Bio import Entrez, Medline
import requests
import io
from datetime import datetime
import calendar
import urllib3

# --- SSL FIX ---
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# --- CONFIGURATION ---
st.set_page_config(page_title="PubLens", page_icon="🔬", layout="wide")

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

# --- HELPER FUNCTIONS ---

def normalize_journal_name(name):
    if not name: return ""
    clean = name.lower().replace(".", "")
    clean = " ".join(clean.split())
    return clean

def get_impact_factor(journal_name):
    if not journal_name: return 0.0
    clean_name = normalize_journal_name(journal_name)
    if clean_name in INTERNAL_METRICS:
        return INTERNAL_METRICS[clean_name]
    return 0.0

def determine_article_type(type_string_or_list):
    if not type_string_or_list: return "Primary Research"
    if isinstance(type_string_or_list, str): type_list = [type_string_or_list]
    else: type_list = type_string_or_list
    for t in type_list:
        if "Review" in str(t): return "Review"
    return "Primary Research"

def generate_citation(row):
    """
    Format: Pachynski RK, et al. Clin Cancer Res. 2021 Jun 15;27(12):3478-3490. DOI: 10.1158/...
    """
    try:
        # Authors: Take first one, add et al.
        authors = str(row.get("Authors", "")).split(",")[0]
        if "," in str(row.get("Authors", "")):
            authors += ", et al"
        
        # Journal: Prefer Abbreviation (TA) -> Full Name (Journal)
        journal = row.get("Abbr", "")
        if not journal or journal == "nan":
            journal = row.get("Journal", "")
            
        date = str(row.get("Date", ""))
        vol = str(row.get("Vol", ""))
        issue = str(row.get("Issue", ""))
        pg = str(row.get("Pages", ""))
        doi = str(row.get("DOI", ""))
        
        # Construct Citation parts
        # Logic: Date;Vol(Issue):Pages.
        details = date
        
        if vol and vol != "nan":
            details += f";{vol}"
            if issue and issue != "nan":
                details += f"({issue})"
        
        if pg and pg != "nan":
            details += f":{pg}"
            
        citation = f"{authors}. {journal}. {details}. DOI: {doi}."
        return citation
    except:
        return "Error generating citation."

# --- SEARCH ENGINE 1: US PUBMED ---

def search_us_pubmed(query, max_results, email, start_date, end_date, exact_phrase, use_full_text):
    Entrez.email = email
    mindate = start_date.strftime("%Y/%m/%d")
    maxdate = end_date.strftime("%Y/%m/%d")
    
    final_query = f'"{query}"' if exact_phrase else query
    
    try:
        if use_full_text:
            search_handle = Entrez.esearch(db="pmc", term=final_query, retmax=max_results, mindate=mindate, maxdate=maxdate, datetype="pdat")
            pmc_results = Entrez.read(search_handle)
            search_handle.close()
            ids = pmc_results["IdList"]
            if not ids: return pd.DataFrame()
            
            link_handle = Entrez.elink(dbfrom="pmc", db="pubmed", id=ids)
            link_results = Entrez.read(link_handle)
            link_handle.close()
            final_ids = []
            for linkset in link_results:
                if "LinkSetDb" in linkset and len(linkset["LinkSetDb"]) > 0:
                    for link in linkset["LinkSetDb"][0]["Link"]:
                        final_ids.append(link["Id"])
        else:
            search_handle = Entrez.esearch(db="pubmed", term=final_query, retmax=max_results, mindate=mindate, maxdate=maxdate, datetype="pdat")
            search_results = Entrez.read(search_handle)
            search_handle.close()
            final_ids = search_results["IdList"]

        if not final_ids: return pd.DataFrame()

        fetch_handle = Entrez.efetch(db="pubmed", id=final_ids, rettype="medline", retmode="text")
        records = Medline.parse(fetch_handle)
        
        data = []
        for record in records:
            j_name = record.get("JT", "")
            j_abbr = record.get("TA", "") # Abbreviation (e.g. Clin Cancer Res)
            
            if_score = get_impact_factor(j_name)
            if if_score == 0.0:
                if_score = get_impact_factor(j_abbr)

            # ROBUST DOI EXTRACTION
            doi = ""
            # Check AID field (Article ID)
            for aid in record.get("AID", []):
                if "[doi]" in aid:
                    doi = aid.replace(" [doi]", "")
                    break
            # Fallback to LID (Location ID)
            if not doi:
                for lid in record.get("LID", []):
                    if "[doi]" in lid:
                        doi = lid.replace(" [doi]", "")
                        break

            data.append({
                "Select": False, # For Checkbox
                "Source": "🇺🇸 US PubMed",
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
                "PMID": record.get("PMID", "")
            })
            
        return pd.DataFrame(data)

    except Exception as e:
        st.error(f"US PubMed Error: {e}")
        return pd.DataFrame()

# --- SEARCH ENGINE 2: EUROPE PMC ---

def search_europe_pmc(query, max_results, start_date, end_date, exact_phrase):
    base_url = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
    s_str = start_date.strftime("%Y-%m-%d")
    e_str = end_date.strftime("%Y-%m-%d")
    date_query = f" AND FIRST_PDATE:[{s_str} TO {e_str}]"
    term = f'"{query}"' if exact_phrase else query
    full_query = f"{term}{date_query}"
    
    params = {"query": full_query, "format": "json", "pageSize": max_results, "resultType": "core"}
    
    try:
        response = requests.get(base_url, params=params, verify=False)
        response.raise_for_status()
        data = response.json()
        result_list = data.get("resultList", {}).get("result", [])
        
        parsed_data = []
        for item in result_list:
            j_name = item.get("journalTitle", "")
            if not j_name: j_name = item.get("journalInfo", {}).get("journal", {}).get("title", "")
            if j_name is None: j_name = ""
            
            # Abbreviation
            j_abbr = item.get("journalInfo", {}).get("journal", {}).get("medlineAbbreviation", "")
            
            # Metadata Extraction
            j_vol = item.get("journalInfo", {}).get("volume", "")
            j_issue = item.get("journalInfo", {}).get("issue", "")
            j_pages = item.get("pageInfo", "")
            
            parsed_data.append({
                "Select": False,
                "Source": "🇪🇺 Europe PMC",
                "Type": determine_article_type(item.get("pubTypeList", {}).get("pubType", [])),
                "Journal": j_name,
                "Abbr": j_abbr,
                "2023 IF": get_impact_factor(j_name),
                "Title": item.get("title", ""),
                "Authors": item.get("authorString", ""),
                "Date": item.get("firstPublicationDate", ""),
                "Vol": j_vol,
                "Issue": j_issue,
                "Pages": j_pages,
                "DOI": item.get("doi", ""),
                "PMID": item.get("pmid", "")
            })
            
        return pd.DataFrame(parsed_data)

    except Exception as e:
        st.error(f"Europe PMC Error: {e}")
        return pd.DataFrame()

def to_excel(df):
    # Drop 'Select' column for Excel export
    df_export = df.drop(columns=["Select"], errors="ignore")
    output = io.BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        df_export.to_excel(writer, index=False, sheet_name='Results')
    return output.getvalue()

# --- UI LAYOUT ---

st.title("🔬 PubLens")
st.markdown("Search **US PubMed** and **Europe PMC**.")

with st.sidebar:
    st.header("Settings")
    email = st.text_input("Email (Required)", placeholder="researcher@lab.edu")
    query = st.text_input("Query", value="Cell DIVE")
    
    st.subheader("Sources")
    use_us = st.checkbox("🇺🇸 US PubMed", value=True)
    use_eu = st.checkbox("🇪🇺 Europe PMC", value=True)
    
    st.subheader("Filters")
    exact = st.checkbox("Exact Phrase Match", value=True)
    full_text_us = st.checkbox("Search Full Text (US PMC)", value=True)
    limit = st.slider("Max Results (Per Source)", 10, 200, 50)
    
    st.divider()
    st.subheader("Date Range")
    today = datetime.today()
    current_year = today.year
    years = list(range(current_year, 1980, -1))
    months = list(calendar.month_name)[1:]

    st.markdown("**Start:**")
    c1, c2 = st.columns(2)
    with c1: s_month = st.selectbox("Month", months, index=0, key="sm")
    with c2: s_year = st.selectbox("Year", years, index=3, key="sy")

    st.markdown("**End:**")
    c3, c4 = st.columns(2)
    with c3: e_month = st.selectbox("Month", months, index=today.month-1, key="em")
    with c4: e_year = st.selectbox("Year", years, index=0, key="ey")

    search_btn = st.button("Search", type="primary")

# --- MAIN LOGIC ---

if search_btn:
    if not email:
        st.warning("⚠️ Please enter your email.")
    elif not use_us and not use_eu:
        st.warning("⚠️ Please select at least one source.")
    else:
        month_map = {m: i+1 for i, m in enumerate(months)}
        start_dt = datetime(s_year, month_map[s_month], 1)
        last_day = calendar.monthrange(e_year, month_map[e_month])[1]
        end_dt = datetime(e_year, month_map[e_month], last_day)

        with st.spinner("Focusing PubLens..."):
            frames = []
            if use_us:
                df_us = search_us_pubmed(query, limit, email, start_dt, end_dt, exact, full_text_us)
                frames.append(df_us)
            if use_eu:
                df_eu = search_europe_pmc(query, limit, start_dt, end_dt, exact)
                frames.append(df_eu)
            
            if frames:
                df_combined = pd.concat(frames, ignore_index=True)
            else:
                df_combined = pd.DataFrame()

            if not df_combined.empty:
                # SCORING & DEDUPING
                df_combined["_score"] = 0
                # Prefer rows with DOIs
                df_combined.loc[df_combined["DOI"].fillna("").astype(str).str.len() > 5, "_score"] += 20
                # Prefer rows with Journal Abbrev (for citations)
                df_combined.loc[df_combined["Abbr"].fillna("").astype(str).str.len() > 1, "_score"] += 10
                # Prefer US source
                df_combined.loc[df_combined["Source"].astype(str).str.contains("US"), "_score"] += 5
                
                df_combined["_clean_doi"] = df_combined["DOI"].fillna("").str.lower().str.strip()
                df_combined = df_combined.sort_values(by="_score", ascending=False)
                
                # Dedup by DOI then Title
                df_dedup = df_combined.drop_duplicates(subset=["_clean_doi"])
                # Remove rows with empty DOIs from result to dedup by title
                # Actually simplest is to just Dedup by Title separately
                df_dedup = df_dedup.drop_duplicates(subset=["Title"])
                
                st.success(f"✅ Found {len(df_dedup)} publications")
                
                # Format IF
                df_dedup["2023 IF"] = pd.to_numeric(df_dedup["2023 IF"], errors='coerce').fillna(0.0)
                df_dedup = df_dedup.sort_values(by="2023 IF", ascending=False)
                
                # Store in session state so selection persists
                st.session_state['results_df'] = df_dedup

# --- RESULTS DISPLAY & CITATION GENERATION ---

if 'results_df' in st.session_state:
    df = st.session_state['results_df']
    
    # Reorder columns for display (Hide metadata cols used for citation but maybe clutter view?)
    # We keep them for now so you can verify data
    display_cols = ["Select", "Source", "Type", "Journal", "2023 IF", "Title", "Date", "DOI", "Abbr", "Vol", "Issue", "Pages"]
    final_cols = [c for c in display_cols if c in df.columns]
    
    # INTERACTIVE EDITOR
    edited_df = st.data_editor(
        df[final_cols],
        use_container_width=True,
        column_config={
            "Select": st.column_config.CheckboxColumn("Select", help="Check to generate citation"),
            "Type": st.column_config.TextColumn("Type", width="small"),
            "2023 IF": st.column_config.NumberColumn("2023 IF", format="%.1f"),
            "DOI": st.column_config.LinkColumn("DOI")
        },
        hide_index=True
    )
    
    col_a, col_b = st.columns([1, 4])
    
    with col_a:
        st.download_button(
            "📥 Download Excel",
            data=to_excel(df),
            file_name=f"PubLens_Results.xlsx"
        )
        
    with col_b:
        # CITATION GENERATOR
        if st.button("📝 Generate Citations for Selected"):
            # Filter for selected rows
            selected_rows = edited_df[edited_df["Select"] == True]
            
            if not selected_rows.empty:
                st.markdown("### Formatted Citations")
                
                # Generate text block
                citation_text = ""
                for index, row in selected_rows.iterrows():
                    # We need to pull the full data from the original DF to get Authors (which might be hidden/truncated)
                    # But here we included Authors in display? No I removed it from 'display_cols' list above to save space?
                    # Wait, I need Authors for citation.
                    # Let's look up the row in the original DF using Title or DOI matches, 
                    # OR just ensure Authors is in the edited_df. 
                    # Authors is NOT in 'final_cols' list above!
                    # Correction: Authors IS needed. I will grab it from session state using Index.
                    
                    # Get original row data using the index from the edited dataframe
                    # (Indices should align if we didn't sort the edited view differently, 
                    # but data_editor preserves index).
                    original_row = df.loc[index]
                    citation_text += generate_citation(original_row) + "\n\n"
                
                st.text_area("Copy Citations", value=citation_text, height=200)
            else:
                st.warning("Please select at least one paper using the checkboxes.")