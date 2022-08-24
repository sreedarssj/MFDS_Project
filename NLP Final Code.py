
# Imports
from datasets import load_dataset
import pandas as pd
import numpy as np
from tqdm import tqdm
tqdm.pandas()

# Load the English STSB dataset
stsb_dataset = load_dataset('stsb_multi_mt', 'en')
stsb_train = pd.DataFrame(stsb_dataset['train'])
#stsb_test=pd.DataFrame(stsb_dataset['test'])
stsb_test = pd.read_csv("MFDS_testdataset.csv")



from sklearn.metrics.pairwise import cosine_similarity
import spacy
nlp = spacy.load("en_core_web_sm")

def text_processing(sentence):
    """
    Lemmatize, lowercase, remove numbers and stop words
    
    """
    sentence = [token.lemma_.lower()          
                for token in nlp(sentence) 
                    if token.is_alpha and not token.is_stop]
    
    return sentence

def cos_sim(sentence1_emb, sentence2_emb):
    """
    Returns:
      The row-wise cosine similarity between the two columns.
      For instance is sentence1_emb=[a,b,c] and sentence2_emb=[x,y,z]
      Then the result is [cosine_similarity(a,x), cosine_similarity(b,y), cosine_similarity(c,z)]
    """
    cos_sim = cosine_similarity(sentence1_emb, sentence2_emb)
    return np.diag(cos_sim)*100

#The jaccard similarity algorithm

import textdistance

def jaccard_sim(row):
    # Text Processing
    sentence1 = text_processing(row['sentence1'])
    sentence2 = text_processing(row['sentence2'])
    print(sentence1)
    
    # Jaccard similarity
    return textdistance.jaccard.normalized_similarity(sentence1, sentence2)*100


stsb_test['Jaccard_score'] = stsb_test.progress_apply(jaccard_sim, axis=1)


#The USE ALgorithm

import tensorflow as tf
import tensorflow_hub as hub

# Load the pre-trained model
gpus = tf.config.list_physical_devices('GPU')
for gpu in gpus:
    # Control GPU memory usage
    tf.config.experimental.set_memory_growth(gpu, True)

module_url = 'https://tfhub.dev/google/universal-sentence-encoder/4'
model = hub.load(module_url)

# Generate Embeddings
sentence1_emb = model(stsb_test['sentence1']).numpy()
sentence2_emb = model(stsb_test['sentence2']).numpy()

# Cosine Similarity
stsb_test['USE_cosine_score'] = cos_sim(sentence1_emb, sentence2_emb)


#The SBERT Algorithm

from sentence_transformers import CrossEncoder

# Load the pre-trained model
model = CrossEncoder('cross-encoder/stsb-roberta-base')

sentence_pairs = []
for sentence1, sentence2 in zip(stsb_test['sentence1'], stsb_test['sentence2']):
    sentence_pairs.append([sentence1, sentence2])
    
stsb_test['SBERT CrossEncoder_score'] = model.predict(sentence_pairs, show_progress_bar=True)*100


score_cols = [col for col in stsb_test.columns if '_score' in col]

#Prints the final result
print(stsb_test.head(11).to_string())

#Stores the final result as a csv file called output
stsb_test.head(11).to_csv('output.csv')




