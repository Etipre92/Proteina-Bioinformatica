import streamlit as st
from stmol import showmol
import py3Dmol
import requests
import biotite.structure.io as bsio
import pandas as pd
import altair as alt
import matplotlib.pyplot as plt
import numpy as np
        
st.sidebar.title('Visualizador de estructuras proteicas en base a una cadena de aminoacidos "SNMFold"')
st.sidebar.write(' El [*ESMFold*](https://esmatlas.com/about) es un software que predice la estructura de una proteína utilizando una cadena de aminoácidos, y esta basado en el modelo de lenguaje ESM-2. Esto con el objetivo de visualizar como sería la estructura esperada de una proteína cualquiera, incluso una inventada. ')

def render(pdb):
    pdbview = py3Dmol.view()
    pdbview.addModel(pdb, 'pdb')
    pdbview.setStyle({'cartoon': {'color': 'spectrum'}})
    pdbview.setBackgroundColor('white')
    pdbview.zoomTo()
    pdbview.zoom(2, 800)
    pdbview.spin(True)
    showmol(pdbview, height=500, width=800)

DEFAULT_SEQ = "MGSSHHHHHHSSGLVPRGSHMRGPNPTAASLEASAGPFTVRSFTVSRPSGYGAGTVYYPTNAGGTVGAIAIVPGYTARQSSIKWWGPRLASHGFVVITIDTNSTLDQPSSRSSQQMAALRQVASLNGTSSSPIYGKVDTARMGVMGWSMGGGGSLISAANNPSLKAAAPQAPWDSSTNFSSVTVPTLIFACENDSIAPVNSSALPIYDSMSRNAKQFLEINGGSHSCANSGNSNQALIGKKGVAWMKRFMDNDTRYSTFACENPNSTRVSDFRTANCSLEDPAANKARKEAELAAATAEQ"
txt = st.sidebar.text_area('Secuencia de entrada:', DEFAULT_SEQ, height=275)

def update(sequence=txt):
    headers = {
        'Content-Type': 'application/x-www-form-urlencoded',
    }

    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', headers=headers, data=sequence, verify=False)
    
    if response.status_code == 200:
        pdb_string = response.content.decode('utf-8')
        
        if "ATOM" not in pdb_string:
            st.error("La respuesta no contenia los datos adecuados para formar una proteina.")
            return
        
        with open('predicted.pdb', 'w') as f:
            f.write(pdb_string)

        try:
            struct = bsio.load_structure('predicted.pdb', extra_fields=["b_factor"])
            b_value = round(struct.b_factor.mean(), 4)

            st.subheader('1.Visualización de la estructura proteica ingresada')
            render(pdb_string)

            st.subheader('Bioinformatica tercer semestre Unison, Ing. Biomédica')
            st.write('Este programa, realizado en la clase de bioinformatica, tiene como objetivo que el estudiante o persona interesada en el tema de la formación y estructura de las proteinas pueda crear o utilizar proteinas existentes y visualizar su estructura.')
            st.info(f'plDDT (Estimador de confianza en la predicción de la proteina del 0-100): {b_value}')

            st.download_button(
                label="Descargar proteina creada",
                data=pdb_string,
                file_name='proteina_creada.pdb',
                mime='text/plain',
            )
            st.write('Programa realizado por los alumnos:')
            st.write('Etienne Preciado Tapia')
            st.write('Victor Beltran Valenzuela')
            st.write('Juan Ruiz Valenzuela')

            def DNA_nucleotide_count(seq):
                d = dict([
                    ('A', seq.count('A')),
                    ('R', seq.count('R')),
                    ('N', seq.count('N')),
                    ('D', seq.count('D')),
                    ('C', seq.count('C')),  
                    ('Q', seq.count('Q')),
                    ('E', seq.count('E')),
                    ('G', seq.count('G')),
                    ('H', seq.count('H')),
                    ('I', seq.count('I')),
                    ('L', seq.count('L')),
                    ('K', seq.count('K')),
                    ('M', seq.count('M')),
                    ('F', seq.count('F')),
                    ('P', seq.count('P')),
                    ('S', seq.count('S')),
                    ('T', seq.count('T')),
                    ('W', seq.count('W')),
                    ('Y', seq.count('Y')),
                    ('V', seq.count('V'))
                ])
                return d

            X = DNA_nucleotide_count(sequence)

            st.subheader('3. Display DataFrame')
            df = pd.DataFrame.from_dict(X, orient='index')
            df = df.rename({0: 'count'}, axis='columns')
            df.reset_index(inplace=True)
            df = df.rename(columns={'index': 'aminoacid'})
            st.write(df)

            st.subheader('4. Display Bar chart')
            p = alt.Chart(df).mark_bar().encode(
                x='aminoacid',
                y='count'
            )

            p = p.properties(
                height=alt.Step(120)  # controls width of bar.
            )
            st.write(p)

        except ValueError as e:
            st.error(f"Error al cargar la estructura de la proteina: {e}")
            st.text_area("PDB File Content", pdb_string)
    else:
        st.error(f"API request failed with status code: {response.status_code}")
        st.text_area("API Response", response.text)

predict = st.sidebar.button('Visualizar', on_click=update)

inventory = {
    "A": "Alanina",
    "R": "Arginina",
    "N": "Asparagina",
    "D": "A. Aspartico",
    "C": "Cisteina",
    "Q": "Glutamina",
    "E": "A. Glutamico",
    "G": "Glicina",
    "H": "Histidina",
    "I": "Isoleucina",
    "L": "Leucina",
    "K": "Lisina",
    "M": "Metionina",
    "F": "Fenilalanina",
    "P": "Prolina",
    "S": "Serina",
    "T": "Treonina",
    "W": "Triptofano",
    "Y": "Tirosina",
    "V": "Valina"
}
st.sidebar.write('Aminoacidos utilizados en la creación de proteinas:')
st.sidebar.dataframe(inventory)

if not predict:
    st.warning('👈 Escribe una secuencia de aminoacidos y oprime visualizar para visualizar tu propia proteina!')

st.header('Ingrese su secuencia de ADN')

sequence_input = "-Cadena de ADN aqui:\nGAACACGTGGAGGCAAACAGGAAGGTGAAGAAGAACTTATCCTATCAGGACGGAAGGTCCTGTGCTCGGG\nATCTTCCAGACGTCGCGACTCTAAATTGCCCCCTCTGAGGTCAAGGAACACAAGATGGTTTTGGAAATGC\nTGAACCCGATACATTATAACATCACCAGCATCGTGCCTGAAGCCATGCCTGCTGCCACCATGCCAGTCCT"

#sequence = st.sidebar.text_area("Sequence input", sequence_input, height=250)
sequence = st.text_area("Secuencia de entrada:", sequence_input, height=250)
sequence = sequence.splitlines()
sequence = sequence[1:] # Skips the sequence name (first line)
sequence = ''.join(sequence) # Concatenates list to string

st.write("""
***
""")

## Prints the input DNA sequence
#st.header('Secuencia de ADN')
#sequence

## DNA nucleotide count
st.header('Composiciòn de los nucleotidos del ADN')


# 1. Print dictionary
col1, col2 = st.columns(2)
def DNA_nucleotide_count(seq):
    d = dict([
        ('A', seq.count('A')),
        ('T', seq.count('T')),
        ('G', seq.count('G')),
        ('C', seq.count('C'))
    ])
    return d

X = DNA_nucleotide_count(sequence)

with col1:
    # 2. Print text
    st.subheader('1. Descripciòn de nucleotidos ')
    st.write('Existen ' + str(X['A']) + ' de adenina (A) en la secuencia de ADN')
    st.write('Existen ' + str(X['T']) + ' de timina (T) en la secuencia de ADN')
    st.write('Existen ' + str(X['G']) + ' de guanina (G) en la secuencia de ADN')
    st.write('Existen ' + str(X['C']) + ' de citocina (C) en la secuencia de ADN')
with col2:
    # 3. Display DataFrame
    st.subheader('2. Tabla de nucleotidos')
    df = pd.DataFrame.from_dict(X, orient='index')
    df = df.rename({0: 'cantidad'}, axis='columns')
    df.reset_index(inplace=True)
    df = df.rename(columns={'index': 'nucleotido'})
    st.write(df)

# Add CSS styling for subheaders
st.markdown(
    """
    <style>
    .stHeader > .deco-btn-container > div {
        display: inline-block;
        margin-right: 20px;
    }
    </style>
    """,
    unsafe_allow_html=True
)

st.subheader('3.Grafica cantidad de Nucletidos')
p = alt.Chart(df).mark_bar().encode(
    x='nucleotido',
    y='cantidad'
).configure_mark(
    opacity=0.8,
    color='violet'
)

p = p.properties(
    width=alt.Step(80)  # controls width of bar.
)
st.write(p)

### 5. Display Pie Chart using Altair
st.subheader('4.Grafica de convergencia')

import altair as alt

# Reshape the data for animated pie chart
df_pivot = df.melt('nucleotido', var_name='metric', value_name='value')

# Create animated pie chart
animated_pie_chart = alt.Chart(df_pivot).mark_arc().encode(
    alt.X('value:Q', stack='zero'),
    color='nucleotido:N',
    tooltip=['nucleotido', 'metric', 'value']
).properties(
    width=500,
    height=400
).transform_joinaggregate(
    total='sum(value)',
    groupby=['nucleotido']
).transform_calculate(
    percentage='datum.value / datum.total'
).encode(
    text=alt.Text('percentage:Q', format='.1%')
).configure_mark(
    opacity=0.8
)

# Adjust properties of the animated pie chart
animated_pie_chart = animated_pie_chart.properties(
    width=300,
    height=300
)

# Display the animated pie chart
st.altair_chart(animated_pie_chart, use_container_width=True)

p = alt.Chart(df).mark_bar().encode(
    x='nucleotido',
    y='cantidad',
    column='nucleotido'
)

p = p.properties(
    width=alt.Step(80),  # controls width of bar
    height=alt.Step(40),  # controls height of bar
    column=alt.Column(
        spacing=10  # controls spacing between grouped bars
    )
)

st.header('Creadores')
st.markdown('- **Etienne Preciado Tapia**')
st.markdown('- **Victor Beltran Valenzuela**')
st.markdown('- **Juan Ruiz Valenzuela**')
