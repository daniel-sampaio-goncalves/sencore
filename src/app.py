import re
import pandas as pd
import numpy as np
import sqlite3
import streamlit as st


def parse_dataset_ids(dataset_id_str: list) -> list:
    """
    def to parse dataset ids with ranges and single numbers
    """
    dataset_list = []
    if dataset_id_str:
        parts = dataset_id_str.split(',')
        for part in parts:
            part = part.strip()
            if re.match(r'^\d+$', part):  # Single numbers
                dataset_list.append(int(part))
            elif re.match(r'^\d+-\d+$', part):  # Range like 1-5
                start, end = map(int, part.split('-'))
                dataset_list.extend(range(start, end + 1))
    return sorted(set(dataset_list))  # Remove duplicates and sort

def coerce_to_float(df_obj: pd.DataFrame) -> pd.DataFrame:
    """
    Coerce numeric columns to proper float dtype (invalid values -> NaN)
    """
    for col_ in ['padj', 'log2FoldChange', 'NC_CTRL', 'NC_Sen']:
        if col_ in df_obj.columns:
            df_obj[col_] = pd.to_numeric(df_obj[col_], errors='coerce')
    return df_obj

def get_Multiple_Genes_dev_extend_log_padj(sql_pathway_:str,dataset_id=None , gene_name=None, extended_search=False, searchmode=None, 
                                           padj_filtering_variable=None,pvalue_threashold=None,
                                           log2_filtering_variable=None, log2_filter_side=None, log2_value=None,
                                           SF_TF_filter=None) ->pd.DataFrame:
    """
    main sql retrieval function with several options
    """
    path_to_db=sql_pathway_
    # Connect to the SQLite database
    conn = sqlite3.connect(path_to_db)

    # Define the view-like query
    query = '''
        SELECT DISTINCT Genes, log2FoldChange, padj, NC_CTRL, NC_Sen, Organism, Cells, Dataset
        FROM Unified_Gene_Table
        WHERE 1=1
    '''
    # Dynamic filters: where 1=1 is a trick so we can add dynamic filters
    params = []
    if extended_search: #begining of genes:
        if dataset_id:
            #this allows for ranges like 1-5 and multiple like 1,5,6
            dataset_list = parse_dataset_ids(dataset_id)
            placeholders = ','.join(['?'] * len(dataset_list))
            query += f' AND Dataset IN ({placeholders})'
            params.extend(dataset_list)
        
        if gene_name and searchmode =="Starts with": 

            gene_list = [name.strip().lower() for name in gene_name.split(',') if name.strip()]
            for gene in gene_list:
                query += f' AND LOWER(Genes) LIKE ?'
                params.append(f'{gene}%')  # Matches gene names starting with the input
        
        if gene_name and searchmode =="Contains": 
            gene_list = [name.strip().lower() for name in gene_name.split(',') if name.strip()]
            for gene in gene_list:
                query += f' AND LOWER(Genes) LIKE ?'
                params.append(f'%{gene}%')  # Matches gene names containing the input

    else:    
        # Handle multiple dataset IDs
        if dataset_id:
            #this chunck allows for ranges like 1-5 and multiple like 1,5,6
            dataset_list = parse_dataset_ids(dataset_id)
            placeholders = ','.join(['?'] * len(dataset_list))
            query += f' AND Dataset IN ({placeholders})'
            params.extend(dataset_list)
        
        # Add gene name filter (multiple supported)
        if gene_name:
            gene_list = [name.strip().lower() for name in gene_name.split(',') if name.strip()]
            placeholders = ','.join(['?'] * len(gene_list))
            query += f' AND LOWER(Genes) IN ({placeholders})'
            params.extend(gene_list)
    #pvalue:
    if padj_filtering_variable:
        query += f" AND padj < {pvalue_threashold}"
    #log2 fc value:
    if log2_filtering_variable:
        if log2_filter_side=="Log2< (Lower than)":
            query += f" AND log2FoldChange < {log2_value}"
        else :
            query += f" AND log2FoldChange > {log2_value}"
    #TF and SF filtering
    if SF_TF_filter=="Transcription Factors":
        query += f" AND LOWER(Genes) IN (SELECT LOWER(Symbol) FROM Transcription_Factors_Table)"
    elif SF_TF_filter=="Secreted Proteins":
        query += f''' AND LOWER(Genes) IN (SELECT LOWER(Symbol_Mouse) FROM Mouse_Secreted_Factors_Table UNION 
                    SELECT LOWER(Symbol_Human) FROM Human_Secreted_Factors_Table
                    )
    '''
    # Execute query and load into pandas DataFrame
    df = pd.read_sql_query(query, conn, params=params)
    # Add a datasetcount table:
    # Count occurrences of each value in column A
    count_series = df['Genes'].str.lower().value_counts()

    # Map the counts back to the DataFrame as column B
    df['Datasets_Counted'] = df['Genes'].str.lower().map(count_series)
    # Close the connection
    conn.close()
    # Replace any invalide numbers with valid ones (eg ''/na/NULL -> NaN)
    df=coerce_to_float(df)
    return df

def retrieve_data_set_table(sql_pathway_:str)-> pd.DataFrame:
    """
    get the reference table names:dataset_ids to use to retrieve genes to specific dataset
    """
    path_to_db=sql_pathway_
    conn = sqlite3.connect(path_to_db)
    query = '''
        SELECT *
        FROM Dataset_Table
    '''
    df_dataset = pd.read_sql_query(query, conn)
    conn.close()
    return df_dataset

def main() -> None:
    sql_pathway="./SENCORE_08_11_2025.db"
    ### Streamlit UI###
    #page layout:
    st.set_page_config(
        page_title="SenCore",
        page_icon="🧬",
        layout="wide",
    )
    st.title("SENCORE Viewer")
    st.divider()
    # Create a form so Enter key can trigger submission
    
    #checkbox for extended search
    extended_search = st.checkbox("Use Extended Genes Search (Pattern Search)", value=False)

    # Checkbox for enabling log filtering
    log_filtering = st.checkbox("Log filtering")

    # Checkbox for padj filtering
    padj_filtering = st.checkbox("padj filtering")

    # Add Radiobutton for TF and SF:
    Radio_TForSF=st.radio("Filter genes for:", ["All (default)", "Transcription Factors", "Secreted Proteins"], index=0)

    # Checkbox for different sets
    st.divider()
    st.text("Optional: Quickly select datasets to include in the search:")
    Set_Number =""
    col1, col2, col3, col4, col5, col6 = st.columns(6)
    # Place a checkbox in each column for different predefined dataset id ranges
    with col1:
        chk1 = st.checkbox("MetaSet")
        if chk1:
            Set_Number="1-18"

    with col2:
        chk2 = st.checkbox("IMR90/Huh7 DOX & IMR90 siEDA2R")
        if chk2:
            if Set_Number!="":
                Set_Number+=",19-21"
            else:
                Set_Number="19-21"
    with col3:
        chk3 = st.checkbox("Muscle cells in vivo")
        if chk3:
            if Set_Number!="":
                Set_Number+=",22-33"
            else:
                Set_Number="22-33"
    with col4:
        chk4 = st.checkbox("AER vs ectoderm")
        if chk4:
            if Set_Number!="":
                Set_Number+=",34"
            else:
                Set_Number="34"
    with col5:
        chk5 = st.checkbox("Muscle fibers p21 high vs low")
        if chk5:
            if Set_Number!="":
                Set_Number+=",35-37"
            else:
                Set_Number="35-37"
    with col6:
        chk6 = st.checkbox("SASP treated cells")
        if chk6:
            if Set_Number!="":
                Set_Number+=",38-44"
            else:
                Set_Number="38-44"



    with st.form(key="gene_query_form"):
        search_mode=None
        padj_value=None
        log_comparison_operator=None
        log_threshold_value=None
        dataset_id = st.text_input("Optional: Enter a dataset ID (OR  multiple (eg. 5,10,18) OR a range (eg. 1-10) OR combined (eg. 1-5,8,10-15,4)", value=Set_Number)
        gene_name = st.text_input("Optional: Search a specific gene (Or multiple: sox2, sox4,...)/ A gene pattern if extended search is enabeld ( eg. sox)")
        

        # Radio button for search mode (only shown if extended search is enabled)
        if extended_search:
            search_mode = st.radio(
                "Extended Search Mode",
                options=["Starts with", "Contains"],
                index=0,
                help="Choose how gene name patterns should be matched"
            )
        if log_filtering:
            log_comparison_operator = st.radio("Select comparison operator", options=["Log2\\> (Higher than)","Log2< (Lower than)"], key="log_threshold_filtering")
            log_threshold_value = st.number_input("Enter log2 value (defaults to 1.5 FOLD)", value=np.log2(1.5), key="log_threshold")

        
        if padj_filtering:
            padj_value = st.number_input("Enter padj value", value=0.05, key="padj_threshold")      
        
        submit = st.form_submit_button("Show Genes")



    with st.sidebar:
        st.sidebar.markdown("### Dataset IDs for reference:")
        data_set_table_sidebar = retrieve_data_set_table(sql_pathway)
        st.dataframe(data_set_table_sidebar, height=1575, hide_index=True)

    # Trigger query when form is submitted
    if submit:
        with st.spinner("Loading... (Large queries may exceed memory)", show_time=True):
            df = get_Multiple_Genes_dev_extend_log_padj(sql_pathway,dataset_id, gene_name, extended_search, search_mode,
                                                        padj_filtering,padj_value,
                                                        log_filtering, log_comparison_operator,log_threshold_value, 
                                                        Radio_TForSF)
            
            # Display filtering criteria
            # adjust column widths based on filtering options
            if log_filtering:
                log_screen_value=0.12
                Genes_screen_value=0.08
            
            else:
                log_screen_value=0.2
                Genes_screen_value=0.10

            col_res1, col_res2, col_res3 = st.columns([Genes_screen_value,log_screen_value,0.2])
            with col_res1:
                    st.write(f"Filtering for {Radio_TForSF.lower()} genes")
            if log_filtering:
                col_log=col_res2
            if padj_filtering:
                if log_filtering:
                    col_padj=col_res3
                else:
                    col_padj=col_res2

            if log_filtering:
                with col_log:
                        st.write(f"Filtering logs where value {log_comparison_operator} {log_threshold_value}")
            if padj_filtering:
                with col_padj:
                        st.write(f"Filtering padj<{padj_value}")

            dataset_id_display = "1-44" if dataset_id == '' or dataset_id == '1-18,19-21,22-33,34,35-37,38-44' else dataset_id
            Genes_display = "all possible genes" if gene_name == '' else gene_name
            styled_df=df.sort_values(by="Dataset", ascending=True).reset_index(drop=True)
            st.write(f"Results for datasets '{dataset_id_display}' and '{Genes_display}' with a total of {len(df)} genes found.")
            st.dataframe(styled_df, height=800)

if __name__=="__main__":
    main()