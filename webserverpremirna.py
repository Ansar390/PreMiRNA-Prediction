from io import StringIO
from Bio import SeqIO
import pandas as pd
import streamlit as st
from PIL import Image
# import extractFeatures as fe
import fet_ext as fe
import numpy as np
import pickle
import base64
from sklearn.preprocessing import StandardScaler


icon = Image.open('fav.png')
st.set_page_config(page_title='MiRNA', page_icon = icon)
# with open("./scaler.pkl", 'rb') as file:
#     scaler = pickle.load(file)


# with open("./model.pkl", 'rb') as file:
#     model = pickle.load(file)
filename = 'ET_Pre_Mirna.sav'
model = pickle.load(open(filename, 'rb'))

preMirna=pd.read_csv('premrna_testX.csv')

# print(preMirna.head(2))


def seqValidator(seq):
    allowed_chars = set('ACGU')
    if set(seq).issubset(allowed_chars):
        return True
    return False

final_df = pd.DataFrame(columns=['Sequence ID', 'Sequence', 'Label'])
seq = ""
st.header("""MiRNA Webserver""")

file_ = open("./WebPic.gif", "rb")
contents = file_.read()
data_url = base64.b64encode(contents).decode("utf-8")
file_.close()
st.markdown(f'<img src="data:image/gif;base64,{data_url}" alt="Header Image">', unsafe_allow_html=True)

st.sidebar.subheader(("Input Sequence(s) (FASTA FORMAT ONLY)"))
fasta_string  = st.sidebar.text_area("Sequence Input", height=200)
            
st.subheader("Click the Example Button for Sample Data or Datset Button to download complete dataset")
with open("dataset.rar", "rb") as fp:
    btn = st.download_button(
        label="Download Dataset",
        data=fp,
        file_name="dataset.rar",
        mime="application/rar"
    )
    



if st.button('Example'):
    st.code(">osa-MIR2919 MI0013261 Oryza sativa miR2919 stem-loop\nGCGGGGAGCCGAUGGGCCGCGGGAGAGGAGAGAGAGAGGGAGAGGGGAGUGGGCCGAGCCGGCCCAAGAGGGAAGGGGGGGGGGGGAAAGAGGACUUUUUAGGAUUUUUCUUUUUAUAAACCUUUUUCCACUUUGUUUAUUUCUUUUCAAUUAUUAUUUG", language="markdown")
    st.code(">mtr-MIR5276 MI0018417 Medicago truncatula miR5276 stem-loop\nAAGCAUAUAUAGGGGGAGCACCUUGCUGGGGCAUGGCUGGAGCUGGUUGGGGGUCGCCCAGACACUCCCCAACAGGCUGCUCCCCCCAAAUCUGCCCCA", language="markdown")

# st.selectbox()

# def mirna()    
if st.sidebar.button("SUBMIT"):
    if(fasta_string==""):
        st.info("Please input the sequence first.")
    fasta_io = StringIO(fasta_string) 
    records = SeqIO.parse(fasta_io, "fasta") 
    for rec in records:
        seq_id = str(rec.id)
        seq=str(rec.seq)
        if(seqValidator(seq)):
            df_temp = pd.DataFrame([[seq_id, seq,'None']], columns=['Sequence ID', 'Sequence','Label'] )
            final_df = pd.concat([final_df,df_temp], ignore_index=True)
            print(final_df)
        else:
            st.info("Sequence with Sequence ID: " + str(seq_id) + " is invalid, containing letters other than standard amino acids")
    fasta_io.close()
    if(final_df.shape[0]!=0):
        for iter in range(final_df.shape[0]):
            temp_seq =  final_df.iloc[iter, 1]
            print(temp_seq,"temp seq")


            
            # vec=fe.calcFV(temp_seq.lower())
            # vec= np.array(vec)
            # print("In array")

            testx = preMirna.drop(['class'],axis=1)

            X1 = testx.to_numpy()
            scaler = StandardScaler().fit(X1)
            # vec = scaler.transform(X1)
            # print("In scaler",vec.shape) .reshape(1, -1)

            fv_array = scaler.transform(np.array(fe.calcFV(temp_seq.lower())).reshape(1, 522))
            # fv_array = np.nan_to_num(fv_array.astype('float64'))
            # fv_array = np.nan_to_num(fv_array)

            fv_array = np.nan_to_num(fv_array.astype('float32'))
            # print("fv_array",fv_array)
            # fv_array=fv_array.reshape(1, -1)

            # print (fv_array,"fv array to model")
            score = model.predict(fv_array)
            pred_label = np.round_(score, decimals=0, out=None)
            if(pred_label==1):
                pred_label="Stress Pre-MIRNA"
            else:
                pred_label="No Stress Pre-MIRNA"
            final_df.iloc[iter, 2] = str(pred_label)

    st.dataframe(final_df)