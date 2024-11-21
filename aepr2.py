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

            st.subheader('2. Print dictionary')
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
            st.write(X)

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
            
            for k, v in X.items():
                X[k] = float(v)
            
            fig1, ax1 = plt.subplots()
            ax1.pie(X.values() , labels=X.keys(), autopct='%1.1f%%', radius=10)
            ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

            st.pyplot(fig1)

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
    st.warning('👈 Escribe la secuencia de aminoacidos para visualizar la proteina!')
