# NOMBRE: Bazantes, Ronald

# Carga de librerías necesarias
from Bio import Entrez
import re
from IPython.core.display import Image
from Bio.Seq import Seq
from Bio.SeqUtils import GC
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
import csv
import re
import pandas as pd
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import math
import seaborn as sns

def download_pubmed(keyword):
    """ Función que extrae listado de artículos desde pubmed a traves de un keyword entre comillas """
    Entrez.email = "ronald.bazantes@est.ikiam.edu.ec"
    Entr=Entrez.read(Entrez.esearch(db="pubmed",
                        term="Cancer",
                        usehistory="y"))
    
    webenv=Entr["WebEnv"]
    query_key=Entr["QueryKey"]
    hand1=Entrez.efetch(db="pubmed",
                        rettype='medline',
                        retmode="text",
                        retstart=0,
                        retmax=543, webenv=webenv, query_key=query_key)
    out_hand1 = open("data/Cancer_pubs.txt", "w")
    m=hand1.read()
    out_hand1.write(m)
    out_hand1.close()
    hand1.close()
    return m

   
def sciense_plot(keyword):
    """Esta funcion nos ayuda a extraer datos según se indique el keyword, en relacion a las mismas palabras usadas en el ejercicio 1, pues  de esta manera se realizó una función en general."""
    
     with open("data/Cancer_pubs.txt", errors="ignore") as l: 
        texto = l.read()
    texto = re.sub(r"\n\s{6}", " ", texto)
    countries_1 = re.findall (r"AD\s{2}-\s[A-Za-z].*,\s([A-Za-z]*)\.\s", texto)
    unique_countries = list(set(countries_1))
    conteo=Counter(countries_1)
    resultado={}
    for clave in conteo:  
        valor=conteo[clave]
        if valor > 1:
            resultado[clave] = valor
    ordenar = (sorted(resultado.values()))
    ordenar.sort(reverse=True)
    import operator
    pais = [] 
    contador = []
    
    reverse = sorted(resultado.items(), key=operator.itemgetter(1), reverse=True)   
    for name in enumerate(reverse):
        pais.append(name[1][0])
        contador.append(resultado[name[1][0]])
    paises_top = pais[0:5]
    frecuencia_cinco = contador [0:5]
    fig = plt.figure(figsize =(10, 7))
    plt.pie(frecuencia_cinco, labels = paises_top)
    (plt.savefig("img/Cancer_pubs.png", dpi=100, bbox_inches='tight'))
    plt.show()
    
    return plt.show()
 
    
    
    
    
    
    
    
    
    



def download_pubmed21(keyword):
    """
    Se descarcará información de PubMed para la obtención de información
    sobre opiniones, metodologías y técnicas empleadas sobre Escherichia coli.   
    """
    Entrez.email = "ronald.bazantes@est.ikiam.edu.ec"
    Entr=Entrez.read(Entrez.esearch(db="pubmed",
                        term="Escherichia coli",
                        usehistory="y"))
    
    webenv=Entr["WebEnv"]
    query_key=Entr["QueryKey"]
    hand1=Entrez.efetch(db="pubmed",
                        rettype='medline',
                        retmode="text",
                        retstart=0,
                        retmax=543, webenv=webenv, query_key=query_key)
    out_hand1 = open("data/E_coli_pubs.txt", "w")
    m=hand1.read()
    out_hand1.write(m)
    out_hand1.close()
    hand1.close()
    return m

def sciense_plot21(keyword):
    """Esta funcion nos ayuda a extraer datos según se indique el keyword, en relacion a las mismas palabras usadas en el ejercicio 1, pues  de esta manera se realizó una función en general."""
    
     with open("data/Escherichia coli_pubs.txt", errors="ignore") as l: 
        texto = l.read()
    texto = re.sub(r"\n\s{6}", " ", texto)
    countries_1 = re.findall (r"AD\s{2}-\s[A-Za-z].*,\s([A-Za-z]*)\.\s", texto)
    unique_countries = list(set(countries_1))
    conteo=Counter(countries_1)
    resultado={}
    for clave in conteo:  
        valor=conteo[clave]
        if valor > 1:
            resultado[clave] = valor
    ordenar = (sorted(resultado.values()))
    ordenar.sort(reverse=True)
    import operator
    pais = [] 
    contador = []
    
    reverse = sorted(resultado.items(), key=operator.itemgetter(1), reverse=True)   
    for name in enumerate(reverse):
        pais.append(name[1][0])
        contador.append(resultado[name[1][0]])
    paises_top = pais[0:5]
    frecuencia_cinco = contador [0:5]
    fig = plt.figure(figsize =(10, 7))
    plt.pie(frecuencia_cinco, labels = paises_top)
    (plt.savefig("img/Escherichia coli_pubs.png", dpi=100, bbox_inches='tight'))
    plt.show()
    
    return plt.show()


def download_pubmed22(keyword):
    """
    Se descarcará información de PubMed para la obtención de información
    sobre opiniones, metodologías y técnicas empleadas sobre el veneno de arañas.   
    """
    Entrez.email = "ronald.bazantes@est.ikiam.edu.ec"
    Entr=Entrez.read(Entrez.esearch(db="pubmed",
                        term="Spider venom",
                        usehistory="y"))
    
    webenv=Entr["WebEnv"]
    query_key=Entr["QueryKey"]
    hand1=Entrez.efetch(db="pubmed",
                        rettype='medline',
                        retmode="text",
                        retstart=0,
                        retmax=543, webenv=webenv, query_key=query_key)
    out_hand1 = open("data/Sp_ven_pubs.txt", "w")
    m=hand1.read()
    out_hand1.write(m)
    out_hand1.close()
    hand1.close()
    return m


def sciense_plot22(keyword):
    """Esta funcion nos ayuda a extraer datos según se indique el keyword, en relacion a las mismas palabras usadas en el ejercicio 1, pues  de esta manera se realizó una función en general."""
    
     with open("data/Spider venom_pubs.txt", errors="ignore") as l: 
        texto = l.read()
    texto = re.sub(r"\n\s{6}", " ", texto)
    countries_1 = re.findall (r"AD\s{2}-\s[A-Za-z].*,\s([A-Za-z]*)\.\s", texto)
    unique_countries = list(set(countries_1))
    conteo=Counter(countries_1)
    resultado={}
    for clave in conteo:  
        valor=conteo[clave]
        if valor > 1:
            resultado[clave] = valor
    ordenar = (sorted(resultado.values()))
    ordenar.sort(reverse=True)
    import operator
    pais = [] 
    contador = []
    
    reverse = sorted(resultado.items(), key=operator.itemgetter(1), reverse=True)   
    for name in enumerate(reverse):
        pais.append(name[1][0])
        contador.append(resultado[name[1][0]])
    paises_top = pais[0:5]
    frecuencia_cinco = contador [0:5]
    fig = plt.figure(figsize =(10, 7))
    plt.pie(frecuencia_cinco, labels = paises_top)
    (plt.savefig("img/Spider venom_pubs.png", dpi=100, bbox_inches='tight'))
    plt.show()
    
    return plt.show()