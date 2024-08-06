import streamlit as st 
from st_aggrid import AgGrid, GridOptionsBuilder
import duckdb
from stmol import showmol
import py3Dmol
import toyplot
import toytree
import toyplot.svg
import os
import pathlib
import boto3
import io
st.set_page_config(layout="wide")


s3 = boto3.resource('s3',
  endpoint_url = os.environ['R2_ENDPOINT_URL'],
  aws_access_key_id = os.environ['R2_KEY_ID'],
  aws_secret_access_key = os.environ['R2_KEY']
)

st.session_state['tree'] = None

def download_protein_s3(id):
    f = io.BytesIO()
    # print(f'{id}.pdb')
    # s3.Bucket('nsdb').download_fileobj('nsdb-000001.pdb', f)
    s3.Bucket('nsdb').download_fileobj(f'{id}.pdb', f)
    # path_to_pdb = pathlib.Path(os.getcwd()).parent / 'structures/pdb/{0}.pdb'.format(id)
    f.seek(0)
    return f.read().decode('utf-8')

@st.cache_resource
def open_database():
    duckdb.execute(        
        "CREATE TABLE ref AS FROM read_csv('reference.csv', AUTO_DETECT=TRUE, header=True)"
    )
    duckdb.execute(
        "CREATE TABLE chainref AS FROM read_csv('chain-reference.csv', AUTO_DETECT=TRUE, header=True)"
    )
    duckdb.execute(
        "CREATE TABLE phylo AS FROM read_csv('phylogenetic-relationships.csv', AUTO_DETECT=TRUE, header=True)"
    )
    # duckdb.read_csv('../data/tree/phylogenetic-relationships.csv', header=True)


@st.cache_resource
def open_tree():
    
    tre = toytree.tree(open('AGNifAlign103.asr.tre').read(), tree_format=1)
    rtre = tre.root(wildcard='BchChl')
    rtre = rtre.drop_tips(wildcard="BchChl")
    return rtre


def find_phylogenetic_relatives(selection):
    if selection['nitrogenase_type'] != 'Anc':
        df = query("SELECT * FROM phylo WHERE x = '{0}'".format(
            selection['nitrogenase_type'] + '_' + selection['scientific_name'].replace(' ', '_')
        ))
    else:
        df = query("SELECT * FROM phylo WHERE x = '{0}'".format(
            selection['scientific_name'].replace('_map', '').replace('_altall', '')
        ))
    return df[['type', 'y']].rename(columns=dict(y='relative'))



open_database()
st.session_state['tree'] = open_tree()

def query(x):
    return duckdb.sql(x).df()

def convert_query_terms(x):
    """
    This function converts a fuzzy query into a structured SQL search term
    """
    if x[:6] == 'taxid:':
        return "SELECT * FROM ref WHERE taxond_id == '{0}'".format(x[6:])
    elif x[:5] == 'nsdb-':
        return "SELECT * FROM ref WHERE id == '{0}'".format(x)
    else:
        return "SELECT * FROM ref WHERE (scientific_name ILIKE '%{0}%' or lineage ILIKE '%{0}%')".format(x)
    
    

def query_components(x):
    return duckdb.sql("SELECT * FROM chainref WHERE id = '{0}'".format(x)).df()



def plot_tree(reference_tips, reference_nodes, query_tip, query_node):
    tip_labels = []
    tip_colors = []
    node_labels = []
    node_colors = []
    for node in st.session_state['tree'].get_node_values('name', False, False):
        if node == query_node:
            node_labels.append(node)
            node_colors.append('red')
        elif node in reference_nodes:
            node_labels.append(node)
            node_colors.append('white')
        else:
            node_labels.append('')
            node_colors.append('gray')
        

    for node in st.session_state['tree'].get_tip_labels():
        if node == query_tip:
            tip_labels.append(node)
            tip_colors.append('red')
        elif node in reference_tips:
            tip_labels.append(node)
            tip_colors.append('gray')
        else:
            tip_labels.append('')
            tip_colors.append('')

    return st.session_state['tree'].draw(
        tip_labels_align=False, tree_style='d', 
        tip_labels=tip_labels, tip_labels_colors=tip_colors,
        node_labels=node_labels, node_colors=node_colors,
        node_sizes=10
    )



col1,  col2 = st.columns(2)
col1.title('Nitrogenase Structure DB')
col1.text("""

Nitrogenases are fundamental for our biosphere, as they 
carry  out the only biological mechanism to fix molecular nitrogen 
into bioavailable nitrogen. However, despite their 
importance, there is little coverage of their structural 
diversity. In this database, we attempt to explore 
such diversity by predicting the structures of many 
extant and ancestral nitrogenases, including variants.

""")
col2.image('nitrospace-pet.png', caption='nitrospace-pet', width=300)

st.divider()

with st.container():

    st.header('Search')
    col1, col2 = st.columns(2)
    query_terms = col1.text_input('vinelandii, taxid:1076, bacterium, Anc_794, nsdb-000006')
    with col2:
        st.write('Nitrogenase type')
        subcol1, subcol2 = st.columns(2)
        include_nif = subcol1.checkbox('Nif', value=True)
        include_vnf = subcol1.checkbox('Vnf', value=True)
        include_anf = subcol2.checkbox('Anf', value=True)
        include_anc = subcol2.checkbox('Anc', value=True)        



        
include_gold = True
include_silver = True
type_includes = [
    ('Nif', include_nif), ('Vnf', include_vnf), 
    ('Anf', include_anf), ('Anc', include_anc)
]

dataset_includes = [
    ('gold', include_gold), ('silver', include_silver), 
]

# st.divider()

if query_terms != '':
    st.header('Results')

    st.write(
        "Gold entries had their MSA built from scratch."
        "Silver entries were built by realigning the main variant (map)"
        "MSA against another variant (e.g. alt2)."
    )

    sql_query = convert_query_terms(query_terms)
    results_df = query(sql_query)[['id', 'nitrogenase_type', 'scientific_name', 'variant', 'status',  'stochiometry', 'lineage']]

    for include in filter(lambda x: x[1] is False, type_includes):
        results_df = results_df.query('nitrogenase_type != "{0}"'.format(include[0]))

    for include in filter(lambda x: x[1] is False, dataset_includes):
        results_df = results_df.query('status != "{0}"'.format(include[0]))

    options_builder = GridOptionsBuilder.from_dataframe(results_df)
    options_builder.configure_selection("single")
    grid_options = options_builder.build()
    ag = AgGrid(results_df, grid_options)


selected = False 
try:
    selection = ag['selected_rows'][0]
    selected = True
except NameError:
    pass
except IndexError:
    pass

if selected:
    st.header('Selection')
    st.subheader('Subunits')

    st.dataframe(query_components(selection['id'])[['chain', 'pLDDT', 'subunit', 'sequence']], width=1200)

    pdb_text = download_protein_s3(selection['id'])

    st.download_button(label='download pdb', data=pdb_text, file_name=selection['id'] + '.pdb')

    col1, col2 = st.columns(2)
    with col1:
        st.subheader('Structure')
        view = py3Dmol.view(height = 600,width=600)
        view.addModelsAsFrames(pdb_text)
        view.setStyle({'cartoon':{'color':'spectrum'}})
        view.zoomTo()
        showmol(view, height = 600,width=600, )
        phylogenetic_relationships = find_phylogenetic_relatives(selection)
        st.dataframe(phylogenetic_relationships, width=600)
    with col2:
        st.subheader('Phylogenetics')
        st.write('_'.join(selection['scientific_name'].split('_')[:2]))
        canvas, axes, mark = plot_tree(
            reference_tips=[
                'Nif_Rhodopseudomonas_palustris',
                'Nif_Azotobacter_vinelandii', 
                'Anf_Azotobacter_vinelandii',
                'Vnf_Azotobacter_vinelandii', 
                'Nif_Clostridium_pasteurianum', 
                'Nif_Methanotorris_igneus',
                'Nif_Geoalkalibacter_ferrihydriticus',
                'Nif_Lucifera_butyrica',
                'Nif_Orenia_metallireducens',
                'Nif_Desulfobacca_acetoxidans'
            ],
            reference_nodes=[
                'anc_1206', 'anc_821', 'anc_1345', 'anc_808', 'anc_1352', 'anc_821'
            ],
            query_tip=selection['nitrogenase_type'] + '_' + selection['scientific_name'].replace(' ', '_'),
            query_node='_'.join(selection['scientific_name'].split('_')[:2]).lower()
        )    
        toyplot.svg.render(canvas, "/tmp/tree-plot.svg")
        
        st.image('/tmp/tree-plot.svg', use_column_width=True)


st.divider()
with st.container():

    st.header('About us')

    st.markdown(
"""
If you like this work, cite us:

    [TBA]

"""

    )